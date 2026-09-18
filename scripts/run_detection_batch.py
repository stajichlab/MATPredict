"""Run one pre-planned detection batch (one SLURM job's worth of genomes).

Reads a JSON file describing the batch (a list of
`{taxid, accession, fasta_path, species}` objects -- the same shape as
`MATPredict.detect.genome_acquisition.AcquiredGenome`, produced by splitting
`plan_batches`'s output, one JSON file per batch), builds the reference FASTA
once, and calls `run_batch` for that batch.

Usage:
    pixi run python scripts/run_detection_batch.py \\
        --genomes-json batch_003.json \\
        --db-root db \\
        --out-dir /bigdata/stajichlab/.../detection_rollout/batch_003

This is the per-job driver `scripts/run_detection_batch.slurm` submits; it is
a thin wrapper because `run_batch` itself (src/MATPredict/detect/batch_runner.py)
holds all of the actual per-genome logic (decompression into $SCRATCH,
run_pipeline call, report writing, per-genome failure isolation).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from MATPredict.detect.batch_runner import GenomeRunFailure, run_batch
from MATPredict.detect.genome_acquisition import AcquiredGenome
from MATPredict.detect.reference_fasta import build_reference_fasta


def _load_genomes(genomes_json: Path) -> list[AcquiredGenome]:
    entries = json.loads(genomes_json.read_text())
    return [
        AcquiredGenome(
            taxid=int(entry["taxid"]),
            accession=entry["accession"],
            fasta_path=Path(entry["fasta_path"]),
            species=entry["species"],
        )
        for entry in entries
    ]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomes-json", required=True, type=Path)
    parser.add_argument("--db-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    genomes = _load_genomes(args.genomes_json)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Built once per job/batch (not once per genome) -- see batch_runner's
    # module docstring for why this matters.
    reference_fasta = build_reference_fasta(args.db_root, out_dir / "_reference.faa")

    failures: list[GenomeRunFailure] = []
    run_batch(genomes, args.db_root, reference_fasta, out_dir, failures=failures)

    print(f"ran {len(genomes)} genomes, {len(failures)} failed -> {out_dir}")
    for failure in failures:
        print(f"FAILED\t{failure.taxid}\t{failure.accession}\t{failure.reason}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
