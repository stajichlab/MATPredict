"""SLURM-sized batch orchestration for the genome-scale detection rollout.

Task 1/2 (`MATPredict.detect.genome_acquisition.acquire_genomes`) resolve each
pilot taxid to an `AcquiredGenome`. This module does two things with that list:

1. `plan_batches` groups genomes into per-job batches sized toward a target
   wall-clock runtime per SLURM job (default 5400s / 1.5h -- see this
   project's global HPCC job-sizing rule: size scatter-gather units toward
   ~1-1.5h of real runtime, not the smallest granularity that still
   parallelizes).
2. `run_batch` actually runs the detection pipeline (`run_pipeline`) for every
   genome in one batch, writing each genome's own report under
   `out_dir/<taxid>_<accession>/`.

**Why `run_pipeline` is called directly, not via the `matpredict detect`
CLI**: `MATPredict.detect.cli._cmd_detect` (read in full before this module was
written) does two things that would be pure per-genome waste in a batch loop
that already knows its `db_root`: it re-resolves `MatpredictConfig.from_env`
and, more importantly, calls `build_reference_fasta(db_root, ...)` -- which
walks every `db/*/*/*/proteins.faa` and rewrites/concatenates all of them --
on *every single invocation*, even though the reference FASTA is identical for
every genome in a batch (it depends only on `db_root`, never on the genome).
Shelling out to the CLI once per genome would also pay `argparse`/subprocess
startup cost N times for no benefit. `run_batch` instead takes `reference_fasta`
as an already-built `Path` (built once by the caller, e.g. one
`build_reference_fasta` call before batches are dispatched) and calls
`run_pipeline` in-process for each genome in the batch.

**Scope (Task 4)**: that build-once design originally forced the batch path to
be unrestricted -- one reference FASTA had to serve genomes with different
taxids, so it held every phylum's curated proteins, and a Mucoromycota-only
rollout still aligned all 181 of them per genome instead of the 19 that can
match. `run_batch`'s `phylum` parameter resolves that without giving up
build-once: when the whole batch is one phylum, the routing decision is the
same for every genome, so it is made once and the restricted query set is
built once from it. With no `phylum` the batch stays mixed-phylum and
self-routing, exactly as before.

**Gzipped genomes**: every `AcquiredGenome.fasta_path` from the local-BFD-library
acquisition path (Task 1) points directly at a read-only `.fa.gz` /
`.masked.fasta.gz` file on shared storage -- no copy was made. `run_pipeline`
shells out to `makeblastdb`/`tblastn`/`exonerate`/`miniprot` (none gzip-aware in
this codebase's usage) and calls `Bio.SeqIO.parse` directly on the path (no
gzip handling), so `run_batch` decompresses each genome into `$SCRATCH`
(node-local, finite, per-job scratch -- never a hardcoded `/scratch/$USER/...`
path, and never the read-only BFD tree or a directory under version control)
immediately before that genome's `run_pipeline` call, and deletes the
decompressed copy immediately after (success or failure) so a batch of several
large genomes never accumulates more than one decompressed copy on scratch at
a time.
"""
from __future__ import annotations

import gzip
import logging
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from MATPredict.detect.family_registry import available_phyla, load_all_families, route
from MATPredict.detect.genome_acquisition import AcquiredGenome
from MATPredict.detect.pipeline import EvidenceFloor
from MATPredict.detect.pipeline import run_pipeline as _default_run_pipeline
from MATPredict.detect.reference_fasta import build_reference_fasta
from MATPredict.detect.report import write_detection_gff3, write_detection_report

logger = logging.getLogger(__name__)

# Seed timing evidence: docs/superpowers/plans/2026-09-17-mat-detection-search-
# localization-benchmark-notes.md, Step 2. A real `matpredict detect` run
# against GCF_000143185.2 (Schizophyllum commune H4-8; ~39 Mb uncompressed
# FASTA, 12 MB gzipped from NCBI RefSeq FTP) completed in 1m36.4s wall clock
# (96.4s) with real tblastn/exonerate/miniprot binaries. This is the single
# real data point this estimator has; it is a seed, not a calibrated model.
_REFERENCE_GZIPPED_BYTES = 12 * 1_000_000  # ~12 MB gzipped, per the benchmark notes
_REFERENCE_SECONDS = 96.4


def estimate_seconds_from_gzipped_size(fasta_path: Path) -> float:
    """Estimate `run_pipeline` wall-clock seconds for one genome from its
    on-disk (gzipped) file size.

    Scaling assumption (v1, documented, not yet recalibrated): linear scaling
    from the one real timing data point on record -- a ~12 MB gzipped genome
    (~39 Mb uncompressed) took 96.4s. Gzipped size is used (rather than
    decompressed size) because it is what `Path.stat()` can read for free,
    before touching the file at all; computing a genome's true decompressed
    size would require decompressing it first, which defeats the purpose of a
    cheap pre-batch estimate. Gzipped size is a reasonable proxy for genome
    size/complexity (assembly size roughly tracks compressed size for
    unannotated fungal genome FASTA, which compresses fairly uniformly), not
    an exact one -- repeat/GC-content differences between assemblies will make
    real per-genome runtime deviate from this estimate. This should be
    recalibrated with real multi-genome batch-run timing data once it exists
    (see Task 3 report).
    """
    gzipped_bytes = fasta_path.stat().st_size
    return _REFERENCE_SECONDS * (gzipped_bytes / _REFERENCE_GZIPPED_BYTES)


def plan_batches(
    genomes: list[AcquiredGenome],
    target_seconds_per_job: int = 5400,
    estimate_seconds: Callable[[Path], float] = estimate_seconds_from_gzipped_size,
) -> list[list[AcquiredGenome]]:
    """Group `genomes` into per-SLURM-job batches sized toward
    `target_seconds_per_job` (default 5400s = 1.5h, this project's HPCC
    scatter-gather sizing target).

    Greedy bin-packing in the input order: genomes are added to the current
    batch until the next genome would push the batch's estimated total past
    the target, at which point a new batch starts. A single genome whose own
    estimate already exceeds the target gets its own one-genome batch rather
    than being merged with anything else or split further (this function
    groups genomes; it never splits a genome's own detection run into
    pieces). Order of genomes within and across batches is preserved from the
    input list -- this is not a bin-packing optimizer, just a cheap forward
    pass that keeps batches close to the target.
    """
    batches: list[list[AcquiredGenome]] = []
    current: list[AcquiredGenome] = []
    current_seconds = 0.0
    for genome in genomes:
        estimated = estimate_seconds(genome.fasta_path)
        if current and current_seconds + estimated > target_seconds_per_job:
            batches.append(current)
            current = []
            current_seconds = 0.0
        current.append(genome)
        current_seconds += estimated
    if current:
        batches.append(current)
    return batches


def _decompressed_dest_path(fasta_path: Path, scratch_dir: Path, tag: str) -> Path:
    """The deterministic scratch destination for `fasta_path`'s decompressed
    copy. `tag` (e.g. `"<taxid>_<accession>"`) is prefixed onto the destination
    filename so two genomes processed in the same batch can never collide on
    scratch, even if their source basenames happen to match. Pure path
    arithmetic, no I/O -- safe to compute before a `try` block so the same
    path is known for cleanup even if the decompression that would write it
    never starts or fails partway through.
    """
    stripped_name = fasta_path.with_suffix("").name  # "foo.fa.gz" -> "foo.fa"
    return scratch_dir / f"{tag}_{stripped_name}"


def _decompress_genome(fasta_path: Path, dest: Path) -> None:
    """Decompress `fasta_path` (a `.gz` file) to the plain-FASTA path `dest`.

    If this raises partway through (truncated/corrupt source, `$SCRATCH` full
    or otherwise unwritable -- both realistic on shared HPCC storage), `dest`
    may exist as a partial file; the caller cleans it up unconditionally in a
    `finally` block using the same deterministic path, regardless of how far
    this function got.
    """
    with gzip.open(fasta_path, "rb") as src, dest.open("wb") as dst:
        shutil.copyfileobj(src, dst)


@dataclass(frozen=True)
class GenomeRunFailure:
    """One genome in a batch whose `run_pipeline` call raised, and why."""

    taxid: int
    accession: str
    reason: str


def run_batch(
    genomes: list[AcquiredGenome],
    db_root: Path,
    reference_fasta: Path,
    out_dir: Path,
    run_pipeline: Callable = _default_run_pipeline,
    failures: list[GenomeRunFailure] | None = None,
    *,
    emit_cds_fasta: bool = False,
    phylum: str | None = None,
    evidence_floor: EvidenceFloor | None = None,
) -> None:
    """Run the detection pipeline for every genome in one batch.

    Each genome's `.gz` FASTA (Task 1's local-library acquisition path never
    decompresses or copies the source file) is decompressed into `$SCRATCH`
    immediately before its `run_pipeline` call and deleted immediately after
    (success or failure), so scratch never holds more than one batch member's
    decompressed genome at a time. `$SCRATCH` is read via `os.environ["SCRATCH"]`
    (not `.get()`) so a missing/unset `$SCRATCH` fails loudly and immediately,
    per this project's global HPCC rule, rather than silently falling back to a
    hardcoded or relative path.

    `emit_cds_fasta` (keyword-only, default False so every existing caller's
    output is unchanged) forwards the ALREADY-DECOMPRESSED scratch copy to
    `write_detection_gff3` as its `genome_fasta`, which makes it additionally
    emit CDS features with `translation=` attributes and a companion
    `detected_loci.fasta`. The decompressed copy is reused deliberately: it
    exists on scratch for the duration of this genome's run anyway, so
    nothing is re-decompressed or refetched, and the `.gz` source is not what
    `run_pipeline` itself was given. The companion FASTA holds full contig
    sequences and is therefore large, which is why this is opt-in.

    Each genome's report is written under `out_dir/<taxid>_<accession>/`
    (`detected_loci.gff3`, `detection_report.yaml`); `accession` values seen in
    practice (e.g. `GCA_004115165.2`) contain `.` but no path separators, and
    `Path` treats `.` as an ordinary filename character on Linux, so no
    additional sanitization is applied here -- this was confirmed against this
    repo's actual acquired accessions before writing this function, not
    assumed.

    A single genome's `run_pipeline` failure does not abort the rest of the
    batch (consistent with `acquire_genomes`'s per-taxid resilience): the
    failure is logged, recorded in `failures` (if a list is passed), and the
    loop continues to the next genome, since one bad genome in a multi-genome
    SLURM job should not lose every other genome's already-computed result.

    **`phylum` (keyword-only, default None) gives the batch an explicit
    scope.** Without it this function self-routes per genome inside
    `run_pipeline` and searches whatever the caller's unrestricted reference
    FASTA holds -- which is why, before this parameter existed, a
    Mucoromycota-only rollout still aligned all 181 curated proteins against
    every one of its genomes even though only 19 of them belong to a
    Mucoromycota family (measured against the live `db/`; see
    `build_reference_fasta`'s docstring and
    `tests/detect/test_reference_fasta.py::test_per_phylum_builds_partition_the_unrestricted_build`).
    When `phylum` IS given, this function routes ONCE with
    `family_registry.route(..., phylum=...)` -- the same helper and the same
    `explicit_phylum` rule `cli._cmd_detect` uses, never a reimplementation --
    and hands that one `RoutingDecision` to every genome's `run_pipeline`
    call, so no genome's taxid can re-route it somewhere else. It also builds
    the tblastn query set ONCE for exactly those families, writing it to the
    `reference_fasta` path it was given: under `phylum`, that argument is the
    DESTINATION for the restricted query set rather than an already-built
    file, because the routed family set is only known here, and letting the
    caller build a second one is how the query set and the routing come to
    disagree. The build-once design is unchanged -- one build for the batch,
    not one per genome.

    An explicit `phylum` matching no curated family raises BEFORE the
    per-genome loop starts. It cannot be raised inside the loop: the loop's
    try/except exists so one bad genome never sinks the batch, and it would
    swallow this error once per genome, letting an 813-genome job run to
    completion having searched nothing and exit as if it had looked. This is
    the same policy `cli._cmd_detect` applies to `--phylum`.

    **`evidence_floor` (keyword-only, default None)** is forwarded to
    `run_pipeline` so a rollout can tighten or loosen the admission bar
    (before this parameter existed, batches silently inherited whatever
    `run_pipeline`'s defaults happened to be, with no opt-out and no way to
    run a comparison). `None` means "use `run_pipeline`'s own default" and is
    implemented by NOT passing the keyword at all, never by constructing an
    `EvidenceFloor()` here: a second copy of the default in this module could
    silently drift from the real one, which is exactly the bug the CLI's
    argparse defaults hit.
    """
    if "SCRATCH" not in os.environ:
        raise RuntimeError(
            "SCRATCH is not set. run_batch must run inside a SLURM job (or with "
            "SCRATCH exported to a real node-local scratch directory) -- it "
            "refuses to guess a scratch path."
        )
    scratch_dir = Path(os.environ["SCRATCH"])

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    routing = None
    if phylum is not None:
        # Routed once for the whole batch, with the same helper the
        # single-genome CLI path uses. `route`'s `explicit_phylum` rule
        # ignores the taxid entirely, so this one decision is correct for
        # every genome in the batch regardless of how each genome's own taxid
        # would have routed.
        routing = route(phylum=phylum, taxid=None, families=load_all_families(db_root))
        if not routing.families:
            raise ValueError(
                f"phylum {phylum} matches no curated family in {db_root}; "
                f"available phyla: {', '.join(available_phyla(db_root)) or 'none'}"
            )
        # The destination directory is created above, so the restricted
        # query set can be written next to the batch's own output.
        Path(reference_fasta).parent.mkdir(parents=True, exist_ok=True)
        reference_fasta = build_reference_fasta(
            db_root, reference_fasta, family_keys={f.key for f in routing.families},
        )

    # Both are passed ONLY when set, never as an explicit `None`/default
    # value. `run_pipeline` treats `routing=None` as "self-route", which is
    # the unchanged mixed-phylum behaviour, and its `evidence_floor` default
    # is the single source of truth for the admission bar -- restating either
    # here would create a second default in this module that can silently
    # drift from the real one.
    scope_kwargs: dict = {}
    if routing is not None:
        scope_kwargs["routing"] = routing
    if evidence_floor is not None:
        scope_kwargs["evidence_floor"] = evidence_floor

    for genome in genomes:
        tag = f"{genome.taxid}_{genome.accession}"
        genome_out_dir = out_dir / tag
        genome_out_dir.mkdir(parents=True, exist_ok=True)

        decompressed_path = _decompressed_dest_path(genome.fasta_path, scratch_dir, tag)
        try:
            # Decompression is inside the try/except, not before it: a
            # corrupt/truncated .gz source or an I/O error on $SCRATCH (both
            # realistic on shared HPCC storage) must be treated exactly like a
            # run_pipeline failure for this genome -- logged, recorded in
            # `failures`, and NOT allowed to propagate out of run_batch and
            # abort every remaining genome in the batch.
            _decompress_genome(genome.fasta_path, decompressed_path)
            outcome = run_pipeline(
                genome_fasta=decompressed_path,
                proteome_fasta=None,
                taxid=genome.taxid,
                db_root=db_root,
                reference_fasta=reference_fasta,
                evidence_diagnostics_path=genome_out_dir / "evidence_diagnostics.jsonl",
                **scope_kwargs,
            )
            # The genome is already decompressed on local scratch for
            # `run_pipeline`, so the CDS/companion-FASTA path costs one more
            # read of a file that is already there -- no second decompression
            # and no refetch. This stays inside the per-genome try/except so a
            # failure here is handled like any other single-genome failure.
            gff3_kwargs = {"genome_fasta": decompressed_path} if emit_cds_fasta else {}
            write_detection_gff3(
                outcome, genome_out_dir / "detected_loci.gff3", **gff3_kwargs,
            )
            write_detection_report(outcome, genome_out_dir / "detection_report.yaml")
        except Exception as exc:  # noqa: BLE001 - one bad genome must not sink the batch
            logger.warning("detection failed for %s: %s", tag, exc)
            if failures is not None:
                failures.append(
                    GenomeRunFailure(taxid=genome.taxid, accession=genome.accession, reason=str(exc))
                )
        finally:
            # `decompressed_path` is a deterministic path computed before the
            # decompression attempt, so this cleans up both a fully-written
            # copy and a partial one left behind by a mid-copy failure; it is
            # a no-op (missing_ok=True) if decompression never got far enough
            # to create the file at all.
            decompressed_path.unlink(missing_ok=True)
