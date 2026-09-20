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
        --out-dir /bigdata/stajichlab/.../detection_rollout/batch_003 \\
        --phylum Mucoromycota

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
from MATPredict.detect.pipeline import EvidenceFloor
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


def _build_parser() -> argparse.ArgumentParser:
    """The driver's CLI.

    The scope and evidence-floor flags are spelled EXACTLY as
    `matpredict detect`'s (`--phylum`, `--min-hits`, `--min-identity`,
    `--require-core-role` / `--no-require-core-role`). An operator reproducing
    one genome's batch result with a single-genome debug run must not have to
    translate flag names between the two entry points, and a divergent
    spelling would be found only part-way through a rollout.

    The floor defaults are read from `EvidenceFloor()` itself rather than
    restated here, for the same reason the CLI does it: a hard-coded default
    in this file would keep imposing the OLD floor on every flagless batch the
    moment `EvidenceFloor`'s defaults changed. That silent inheritance is how
    the stricter Task 3 defaults reached the batch path unnoticed.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomes-json", required=True, type=Path)
    parser.add_argument("--db-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--phylum",
        required=False,
        help=(
            "Restrict the whole batch to one phylum's curated families, skipping "
            "per-genome taxid routing. This also narrows the tblastn query set to "
            "that phylum's proteins (measured against the live database: 19 for "
            "Mucoromycota versus 181 unrestricted), which is the point -- every "
            "protein outside the routed families is alignment time paid per genome "
            "for hits that can only be spurious. Choices are not validated here; "
            "run_batch rejects a phylum matching no curated family before it "
            "processes any genome."
        ),
    )
    _floor_defaults = EvidenceFloor()
    parser.add_argument(
        "--min-hits", type=int, default=_floor_defaults.min_hits,
        help="Minimum DISTINCT genes a family needs in a cluster to reach Stage 2 "
             f"polishing (default: {_floor_defaults.min_hits})",
    )
    parser.add_argument(
        "--min-identity", type=float, default=_floor_defaults.min_identity,
        help="Minimum best-hit percent identity for a family to reach Stage 2 "
             "polishing (default: no identity cutoff)",
    )
    parser.add_argument(
        "--require-core-role", action=argparse.BooleanOptionalAction,
        default=_floor_defaults.require_core_role,
        help="Require at least one core_MAT hit for a family to reach Stage 2 "
             f"polishing (default: {_floor_defaults.require_core_role})",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    genomes = _load_genomes(args.genomes_json)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Built once per job/batch (not once per genome) -- see batch_runner's
    # module docstring for why this matters. Under `--phylum` it is `run_batch`
    # that does the single build, because only it knows which families the
    # phylum routed to; building an unrestricted one here first would walk
    # every curated proteins.faa to produce a file nothing reads, and risk the
    # query set disagreeing with the routing.
    reference_fasta = out_dir / "_reference.faa"
    if args.phylum is None:
        reference_fasta = build_reference_fasta(args.db_root, reference_fasta)

    failures: list[GenomeRunFailure] = []
    run_batch(
        genomes, args.db_root, reference_fasta, out_dir, failures=failures,
        phylum=args.phylum,
        evidence_floor=EvidenceFloor(
            min_hits=args.min_hits, min_identity=args.min_identity,
            require_core_role=args.require_core_role,
        ),
    )

    print(f"ran {len(genomes)} genomes, {len(failures)} failed -> {out_dir}")
    for failure in failures:
        print(f"FAILED\t{failure.taxid}\t{failure.accession}\t{failure.reason}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
