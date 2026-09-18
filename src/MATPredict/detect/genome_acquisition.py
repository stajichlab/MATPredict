"""Genome acquisition for the genome-scale detection rollout's pilot target list.

Acquires one representative/reference genome assembly per taxid via the NCBI
`datasets` command-line tool. `datasets` is not a pixi dependency of this project;
on UCR HPCC it is provided by the environment module `ncbi_datasets/18.30.1`
(`module load ncbi_datasets/18.30.1`), which must be loaded in the shell/job that
calls `acquire_genomes` before this module's `_default_runner` can find `datasets`
on PATH. `datasets` handles assembly discovery, download, and unzip-ready packaging
in one tool, so no `NcbiClient` (esummary/FTP) fallback path is implemented here.
"""
from __future__ import annotations

import json
import logging
import subprocess
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from MATPredict.db.ncbi_client import NcbiClient

logger = logging.getLogger(__name__)


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
) -> list[AcquiredGenome]:
    """Acquire one representative/reference genome assembly per taxid.

    For each taxid, shells out to `datasets download genome taxon <taxid>
    --reference --include genome`, unzips the resulting package under `out_dir`,
    and reads the package's `assembly_data_report.jsonl` to recover the resolved
    accession and organism name, then locates that assembly's `*_genomic.fna`.

    `out_dir` is caller-supplied (e.g. a `$SCRATCH`-based path from Task 3's batch
    orchestrator) and is created if it does not already exist. This function never
    hardcodes a scratch or shared-storage path itself.

    `ncbi` is accepted for interface parity with other MATPredict lookup call
    sites and is reserved for a future `NcbiClient`-based (esummary/FTP) fallback
    path if `datasets` ever becomes unavailable; the `datasets`-CLI path
    implemented here does not use it.

    A taxid that cannot be resolved or downloaded is skipped, not raised: it is
    recorded in `failures` (if a list is passed) and logged as a warning, so a
    caller iterating a pilot list still gets every genome that *could* be
    acquired, plus visibility into which ones did not come through.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genomes: list[AcquiredGenome] = []
    for taxid in taxids:
        try:
            genomes.append(_acquire_one(taxid, out_dir, runner))
        except AcquisitionError as exc:
            logger.warning("genome acquisition failed for taxid=%s: %s", taxid, exc)
            if failures is not None:
                failures.append(AcquisitionFailure(taxid=taxid, reason=str(exc)))
    return genomes


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
