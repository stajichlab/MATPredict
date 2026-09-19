"""argparse wiring for the `matpredict detect` subcommand."""
from __future__ import annotations

import argparse
from pathlib import Path

from MATPredict.config import MatpredictConfig
from MATPredict.detect.benchmark import run_benchmark
from MATPredict.detect.pipeline import EvidenceFloor, run_pipeline
from MATPredict.detect.family_registry import load_all_families
from MATPredict.detect.reference_fasta import build_reference_fasta
from MATPredict.detect.report import write_detection_gff3, write_detection_report
from MATPredict.detect.rollout_aggregate import aggregate_reports, write_rollout_summary
from MATPredict.detect.scope_audit import audit_scope, record_taxids_by_family

_ROLLOUT_REPORT_FILENAME = "detection_report.yaml"


def _cmd_detect(args: argparse.Namespace) -> int:
    if not args.genome or not args.out_dir:
        raise ValueError("--genome and --out-dir are required for `matpredict detect`")
    config = MatpredictConfig.from_env(repo_root=Path.cwd())
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reference_fasta = build_reference_fasta(config.db_root, out_dir / "_reference.faa")
    evidence_floor = EvidenceFloor(
        min_hits=args.min_hits, min_identity=args.min_identity,
        require_core_role=args.require_core_role,
    )
    outcome = run_pipeline(
        genome_fasta=Path(args.genome),
        proteome_fasta=Path(args.proteins) if args.proteins else None,
        taxid=args.taxid,
        db_root=config.db_root,
        reference_fasta=reference_fasta,
        evidence_floor=evidence_floor,
        evidence_diagnostics_path=Path(args.evidence_diagnostics) if args.evidence_diagnostics else None,
    )

    write_detection_gff3(outcome, out_dir / "detected_loci.gff3")
    write_detection_report(outcome, out_dir / "detection_report.yaml")
    print(
        f"detected {len(outcome.results)} candidate locus/loci "
        f"({len(outcome.families_attempted)} families attempted) -> {out_dir}"
    )
    # Sub-floor families are reported, never silently dropped (spec section 3).
    for entry in outcome.not_detected:
        print(
            f"not detected\t{entry.family_key.phylum}:{entry.family_key.locus_name}\t{entry.reason}"
        )
    return 0


def _cmd_detect_benchmark(args: argparse.Namespace) -> int:
    config = MatpredictConfig.from_env(repo_root=Path.cwd())
    results = run_benchmark(config.db_root)
    for r in results:
        sens = f"{r.sensitivity:.2f}" if r.sensitivity is not None else "n/a"
        print(f"{r.family_key.phylum}:{r.family_key.locus_name}\tn={r.n_reference_after_holdout}\tsensitivity={sens}\t{r.note}")
    return 0


def _cmd_detect_rollout_summary(args: argparse.Namespace) -> int:
    if not args.reports_dir or not args.out:
        raise ValueError(
            "--reports-dir and --out are required for `matpredict detect rollout-summary`"
        )
    reports_dir = Path(args.reports_dir)
    if not reports_dir.is_dir():
        raise ValueError(f"--reports-dir {reports_dir} is not a directory")

    # One expected report path per genome directory Task 3's run_batch
    # created (`out_dir/<taxid>_<accession>/`) -- built here rather than
    # inside `aggregate_reports` so a genome whose pipeline failed (leaving
    # an empty directory, per Task 3's known gap) is still counted as an
    # attempted genome rather than silently dropped from `total_genomes`.
    report_paths = sorted(
        genome_dir / _ROLLOUT_REPORT_FILENAME
        for genome_dir in reports_dir.iterdir()
        if genome_dir.is_dir()
    )

    summary = aggregate_reports(report_paths)
    write_rollout_summary(summary, Path(args.out))
    print(
        f"aggregated {summary.total_genomes} genome(s) "
        f"({len(summary.genome_errors)} with no report, "
        f"{len(summary.not_detected)} not-detected entries, "
        f"{len(summary.anomalies)} anomalies) -> {args.out}"
    )
    return 0


def _cmd_audit_scope(args: argparse.Namespace) -> int:
    db_root = Path(args.db_root)
    families = load_all_families(db_root)
    record_taxids = record_taxids_by_family(db_root)
    results = audit_scope(families, record_taxids)
    any_uncovered = False
    for result in sorted(results, key=lambda r: (-len(r.uncovered_taxids), r.family_key.phylum, r.family_key.locus_name)):
        if not result.uncovered_taxids:
            continue
        any_uncovered = True
        pct = 100.0 * len(result.uncovered_taxids) / result.total_records
        print(
            f"{result.family_key.phylum}:{result.family_key.locus_name}: "
            f"{len(result.uncovered_taxids)}/{result.total_records} records ({pct:.0f}%) "
            f"uncovered by taxonomic_scope -- recommend {result.recommended_scope_taxid}"
        )
    if not any_uncovered:
        print("every family's own curated records are covered by its taxonomic_scope")
    return 0


def register_subcommands(subparsers: argparse._SubParsersAction) -> None:
    """Register `detect` and its `benchmark` action onto the top-level parser.

    `matpredict detect --genome ... --out-dir ...` (no action) remains the default
    detection behavior for backward compatibility with Task 11's CLI shape; `matpredict
    detect benchmark` is a nested action alongside it, mirroring curate-db's `action`
    subparser pattern in `src/MATPredict/db/cli.py`. --genome/--out-dir are declared
    optional at the argparse level (since `benchmark` doesn't need them) and enforced
    by `_cmd_detect` itself when no action is given.
    """
    detect = subparsers.add_parser("detect", help="Detect MAT loci in a genome")
    detect.add_argument("--genome", required=False)
    detect.add_argument("--proteins", required=False)
    detect.add_argument("--taxid", required=False, type=int)
    detect.add_argument("--out-dir", required=False)
    detect.add_argument("--evidence-diagnostics", required=False)
    detect.add_argument("--min-hits", type=int, default=1)
    detect.add_argument("--min-identity", type=float, default=None)
    detect.add_argument("--require-core-role", action="store_true")
    detect.set_defaults(func=_cmd_detect)

    action = detect.add_subparsers(dest="detect_action")
    benchmark = action.add_parser("benchmark", help="Run the leave-one-out Sn/Sp benchmark suite")
    benchmark.set_defaults(func=_cmd_detect_benchmark)

    rollout_summary = action.add_parser(
        "rollout-summary",
        help="Consolidate a batch's per-genome detection_report.yaml files into one summary",
    )
    rollout_summary.add_argument("--reports-dir", required=False)
    rollout_summary.add_argument("--out", required=False)
    rollout_summary.set_defaults(func=_cmd_detect_rollout_summary)

    audit_scope_parser = action.add_parser(
        "audit-scope",
        help="Audit every curated record's taxid against its own family's taxonomic_scope",
    )
    audit_scope_parser.add_argument("--db-root", default="db")
    audit_scope_parser.set_defaults(func=_cmd_audit_scope)
