"""argparse wiring for the `matpredict detect` subcommand."""
from __future__ import annotations

import argparse
from pathlib import Path

from MATPredict.config import MatpredictConfig
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.reference_fasta import build_reference_fasta
from MATPredict.detect.report import write_detection_gff3, write_detection_report


def _cmd_detect(args: argparse.Namespace) -> int:
    config = MatpredictConfig.from_env(repo_root=Path.cwd())
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_fasta = build_reference_fasta(config.db_root, out_dir / "_reference.faa")
    results = run_pipeline(
        genome_fasta=Path(args.genome),
        proteome_fasta=Path(args.proteins) if args.proteins else None,
        taxid=args.taxid,
        db_root=config.db_root,
        reference_fasta=reference_fasta,
    )

    write_detection_gff3(results, out_dir / "detected_loci.gff3")
    write_detection_report(results, out_dir / "detection_report.yaml")
    print(f"detected {len(results)} candidate locus/loci -> {out_dir}")
    return 0


def register_subcommands(subparsers: argparse._SubParsersAction) -> None:
    detect = subparsers.add_parser("detect", help="Detect MAT loci in a genome")
    detect.add_argument("--genome", required=True)
    detect.add_argument("--proteins", required=False)
    detect.add_argument("--taxid", required=False, type=int)
    detect.add_argument("--out-dir", required=True)
    detect.set_defaults(func=_cmd_detect)
