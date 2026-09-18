"""Genome acquisition for the genome-scale detection rollout's pilot target list.

Acquisition tries two sources, in order, per taxid:

1. **Local BFD genome library** (`local_manifest_path` / `local_library_root`):
   a shared-storage collection of already-downloaded, unannotated fungal genome
   assemblies this lab already maintains for a different project (BFD). Preferred
   because it has no network dependency, is already-vetted lab data, and already
   matches this rollout's "unannotated genome" requirement. Both paths are
   *parameters with sensible on-HPCC defaults*, not hardcoded requirements: if
   either path does not exist (e.g. running outside this HPCC environment), the
   local-library check is skipped entirely and every taxid falls through to (2)
   -- this is not an error.
2. **NCBI `datasets` CLI** (fallback): `datasets` is not a pixi dependency of
   this project; on UCR HPCC it is provided by the environment module
   `ncbi_datasets/18.30.1` (`module load ncbi_datasets/18.30.1`), which must be
   loaded in the shell/job that calls `acquire_genomes` before this module's
   `_default_runner` can find `datasets` on PATH. `datasets` handles assembly
   discovery, download, and unzip-ready packaging in one tool, so no
   `NcbiClient` (esummary/FTP) fallback path is implemented here.
"""
from __future__ import annotations

import csv
import json
import logging
import re
import subprocess
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from MATPredict.db.ncbi_client import NcbiClient

logger = logging.getLogger(__name__)

# Real, already-downloaded, unannotated fungal genome assemblies from a sibling
# lab project (BFD), read-only from this project's point of view -- never write
# into this tree. Overridable per-call for environments where it is not mounted.
DEFAULT_LOCAL_LIBRARY_ROOT = Path(
    "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
)
DEFAULT_LOCAL_MANIFEST_PATH = Path(
    "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/backup_samples.csv"
)

_ASMID_ACCESSION_RE = re.compile(r"^(GC[AF]_\d+\.\d+)_")


@dataclass(frozen=True)
class AcquiredGenome:
    """One successfully-acquired genome assembly, ready for downstream detection.

    This exact shape (`taxid`, `accession`, `fasta_path`, `species`) is consumed
    unchanged by Task 3's batch orchestrator -- do not add, rename, or reorder
    fields without updating that consumer.
    """

    taxid: int
    accession: str
    fasta_path: Path
    species: str


@dataclass(frozen=True)
class AcquisitionFailure:
    """A taxid from the input list that could not be acquired, and why."""

    taxid: int
    reason: str


class AcquisitionError(RuntimeError):
    """Raised internally when a single taxid's genome cannot be acquired.

    Never escapes `acquire_genomes` -- it is caught per-taxid and converted into
    an `AcquisitionFailure` plus a logged warning, so one bad taxid in a pilot
    list does not abort acquisition of the rest.
    """


DatasetsRunner = Callable[[list[str]], "subprocess.CompletedProcess[str]"]


def _default_runner(cmd: list[str]) -> "subprocess.CompletedProcess[str]":
    """Run `datasets` as a real subprocess. Requires `datasets` on PATH (see module
    docstring: `module load ncbi_datasets/18.30.1` on UCR HPCC)."""
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def acquire_genomes(
    taxids: list[int],
    out_dir: Path,
    ncbi: NcbiClient | None = None,
    runner: DatasetsRunner = _default_runner,
    failures: list[AcquisitionFailure] | None = None,
    local_library_root: Path = DEFAULT_LOCAL_LIBRARY_ROOT,
    local_manifest_path: Path = DEFAULT_LOCAL_MANIFEST_PATH,
) -> list[AcquiredGenome]:
    """Acquire one representative genome assembly per taxid, preferring an
    already-downloaded local copy over a fresh NCBI download.

    For each taxid:

    1. Look it up in the local BFD genome library (`local_manifest_path`, a CSV
       joined on its `NCBI_TAXONID` column; `local_library_root`, the directory
       holding the actual `<ASMID>.fa.gz` / `<ASMID>.masked.fasta.gz` files). If
       found, that file's existing path is returned directly -- **no copy or
       download is performed**, this genome may be multi-GB. If either the
       manifest or the library root does not exist on this filesystem, the local
       check is skipped for every taxid (not an error) and NCBI is used instead.
    2. Otherwise, fall back to `datasets download genome taxon <taxid> --reference
       --include genome`, unzip the resulting package under `out_dir`, and read
       the package's `assembly_data_report.jsonl` to recover the resolved
       accession and organism name, then locate that assembly's `*_genomic.fna`.

    `out_dir` is caller-supplied (e.g. a `$SCRATCH`-based path from Task 3's batch
    orchestrator) and is created if it does not already exist; it is only used by
    the NCBI fallback path (the local path never writes anywhere). This function
    never hardcodes a scratch path itself.

    `ncbi` is accepted for interface parity with other MATPredict lookup call
    sites and is reserved for a future `NcbiClient`-based (esummary/FTP) fallback
    path if `datasets` ever becomes unavailable; neither acquisition path
    implemented here uses it.

    A taxid that cannot be resolved by either source is skipped, not raised: it
    is recorded in `failures` (if a list is passed) and logged as a warning, so a
    caller iterating a pilot list still gets every genome that *could* be
    acquired, plus visibility into which ones did not come through.
    """
    out_dir = Path(out_dir)
    genomes: list[AcquiredGenome] = []
    for taxid in taxids:
        try:
            genome = _acquire_local(taxid, local_library_root, local_manifest_path)
            if genome is None:
                out_dir.mkdir(parents=True, exist_ok=True)
                genome = _acquire_one(taxid, out_dir, runner)
            genomes.append(genome)
        except AcquisitionError as exc:
            logger.warning("genome acquisition failed for taxid=%s: %s", taxid, exc)
            if failures is not None:
                failures.append(AcquisitionFailure(taxid=taxid, reason=str(exc)))
    return genomes


def _acquire_local(
    taxid: int, library_root: Path, manifest_path: Path
) -> AcquiredGenome | None:
    """Resolve a taxid against the local BFD genome library. Returns `None` (a
    local miss, not a failure) when the manifest/library aren't mounted, when no
    manifest row matches the taxid, or when the matched row's files aren't
    actually present on disk -- in every `None` case the caller falls through to
    the NCBI path instead.

    Tie-break when multiple manifest rows match a taxid (common: a single
    species can have hundreds of local strain assemblies): prefer a row whose
    `ASMID` starts with `GCF_` (RefSeq/reference-quality) over `GCA_`, then the
    alphabetically-first `ASMID` as a final, deterministic tiebreak. This is an
    arbitrary-but-documented and reproducible choice, not "whatever the CSV
    iteration order gives."

    File preference: when both `<ASMID>.fa.gz` (unmasked) and
    `<ASMID>.masked.fasta.gz` (soft-masked) exist for the chosen row, the
    unmasked file is preferred -- `detect/search.py`'s tblastn/exonerate/miniprot
    calls have no soft-mask-awareness (no masking-related handling found there),
    so the full, unmasked sequence is the safer default for search sensitivity.

    A malformed manifest row (missing an expected column, e.g. `ASMID`) or an
    I/O error reading the manifest (plausible on shared storage: permissions, a
    transient NFS hiccup) is treated the same as a local miss -- logged and
    `None` returned, so the caller falls through to the NCBI path -- rather than
    letting `KeyError`/`OSError`/`csv.Error` propagate out of this function and
    abort acquisition for every remaining taxid in the batch. A manifest/library
    problem does not mean the taxid itself is unacquirable, so falling through
    to NCBI (rather than immediately recording a hard `AcquisitionFailure`) is
    the more useful behavior -- NCBI still gets a chance to resolve the taxid.
    """
    if not manifest_path.exists() or not library_root.exists():
        return None

    try:
        with manifest_path.open(newline="") as fh:
            rows = [row for row in csv.DictReader(fh) if row.get("NCBI_TAXONID") == str(taxid)]
        if not rows:
            return None
        rows.sort(key=lambda row: (not row["ASMID"].startswith("GCF_"), row["ASMID"]))
        chosen = rows[0]
        asmid = chosen["ASMID"]
    except (OSError, KeyError, csv.Error) as exc:
        logger.warning(
            "local library manifest %s could not be read/parsed for taxid=%s (%s); "
            "falling through to NCBI",
            manifest_path, taxid, exc,
        )
        return None

    unmasked_path = library_root / f"{asmid}.fa.gz"
    masked_path = library_root / f"{asmid}.masked.fasta.gz"
    if unmasked_path.exists():
        fasta_path = unmasked_path
    elif masked_path.exists():
        fasta_path = masked_path
    else:
        logger.warning(
            "local library manifest lists %s for taxid=%s but no genome file found under %s",
            asmid, taxid, library_root,
        )
        return None

    accession_match = _ASMID_ACCESSION_RE.match(asmid)
    accession = accession_match.group(1) if accession_match else asmid
    species = chosen.get("SPECIES") or chosen.get("SPECIES_IN") or asmid

    return AcquiredGenome(taxid=taxid, accession=accession, fasta_path=fasta_path, species=species)


def _acquire_one(taxid: int, out_dir: Path, runner: DatasetsRunner) -> AcquiredGenome:
    taxid_dir = out_dir / str(taxid)
    taxid_dir.mkdir(parents=True, exist_ok=True)
    zip_path = taxid_dir / "ncbi_dataset.zip"

    result = runner(
        [
            "datasets",
            "download",
            "genome",
            "taxon",
            str(taxid),
            "--reference",
            "--include",
            "genome",
            "--filename",
            str(zip_path),
        ]
    )
    if result.returncode != 0 or not zip_path.exists():
        stderr = (result.stderr or "").strip()
        raise AcquisitionError(
            f"datasets download failed (exit {result.returncode}): {stderr or 'no output'}"
        )

    extract_dir = taxid_dir / "extracted"
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(extract_dir)
    except zipfile.BadZipFile as exc:
        raise AcquisitionError(f"downloaded package for taxid {taxid} is not a valid zip: {exc}") from exc

    report_path = extract_dir / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
    if not report_path.exists():
        raise AcquisitionError(f"no assembly_data_report.jsonl in package for taxid {taxid}")

    first_line = report_path.read_text().splitlines()[:1]
    if not first_line or not first_line[0].strip():
        raise AcquisitionError(f"empty assembly_data_report.jsonl for taxid {taxid}")
    report = json.loads(first_line[0])

    accession = report.get("accession")
    species = (report.get("organism") or {}).get("organismName")
    if not accession or not species:
        raise AcquisitionError(f"assembly report missing accession/organism for taxid {taxid}")

    assembly_dir = extract_dir / "ncbi_dataset" / "data" / accession
    fasta_candidates = sorted(assembly_dir.glob("*_genomic.fna")) if assembly_dir.is_dir() else []
    if not fasta_candidates:
        raise AcquisitionError(f"no genomic FASTA found for {accession} (taxid {taxid})")

    return AcquiredGenome(
        taxid=taxid,
        accession=accession,
        fasta_path=fasta_candidates[0],
        species=species,
    )
