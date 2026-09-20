"""argparse wiring for the `matpredict detect` subcommand."""
from __future__ import annotations

import argparse
from pathlib import Path

from MATPredict import logger
from MATPredict.config import MatpredictConfig
from MATPredict.detect.benchmark import run_benchmark
from MATPredict.detect.pipeline import EvidenceFloor, run_pipeline
from MATPredict.detect.family_registry import available_phyla, load_all_families, route
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

    # Route FIRST, then build the reference FASTA from the routed families
    # only. The order matters: the reference FASTA is the tblastn query set, so
    # building it before routing (as this did) means every run pays to align
    # every phylum's curated proteins no matter how narrowly it routed. The
    # same `RoutingDecision` object is handed to `run_pipeline` so it cannot
    # re-route to a different set than the query set was built for, and so the
    # taxonomy lookup happens once per run.
    routing = route(args.taxid, load_all_families(config.db_root), phylum=args.phylum)
    # An explicitly requested phylum that matches no curated family is a usage
    # error, not a valid empty result. Left to run it would build an empty
    # reference FASTA, search nothing, and exit 0 with an empty
    # `families_attempted` -- which reads as "looked and found nothing" rather
    # than "never looked". `--phylum`'s argparse choices normally prevent this,
    # but they degrade to unconstrained when the database root cannot be read
    # (see `_phylum_choices`), so the check is enforced here too.
    if args.phylum and not routing.families:
        raise ValueError(
            f"--phylum {args.phylum} matches no curated family in {config.db_root}; "
            f"available phyla: {', '.join(available_phyla(config.db_root)) or 'none'}"
        )
    reference_fasta = build_reference_fasta(
        config.db_root, out_dir / "_reference.faa",
        family_keys={f.key for f in routing.families},
    )
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
        routing=routing,
    )

    # `genome_fasta` is passed ONLY when asked for: it is what makes
    # `write_detection_gff3` additionally emit CDS features with
    # `translation=` attributes and a companion FASTA, and the companion
    # FASTA holds every referenced contig's FULL sequence, which is large.
    gff3_kwargs = {"genome_fasta": Path(args.genome)} if args.emit_cds_fasta else {}
    write_detection_gff3(outcome, out_dir / "detected_loci.gff3", **gff3_kwargs)
    write_detection_report(outcome, out_dir / "detection_report.yaml")
    print(
        f"detected {len(outcome.results)} candidate locus/loci "
        f"({len(outcome.families_attempted)} families attempted, "
        f"routing={routing.routing_mode}) -> {out_dir}"
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


def _phylum_choices() -> list[str] | None:
    """The `--phylum` choices, read from the configured database root when the
    parser is built -- never a hardcoded list, so a phylum added to `db/`
    becomes selectable with no code change here.

    Returns None (argparse: accept any string) rather than raising if the
    database root cannot be read at all. `matpredict --help` and `matpredict
    --version` must keep working in a directory with no database, and refusing
    to build the parser here would break every OTHER subcommand too -- this
    runs at parser-build time, before argparse has even seen which subcommand
    was asked for. `available_phyla` already skips an individual unparseable
    `order.yml`; this broader guard covers anything else config resolution or
    the glob can raise, and logs rather than swallowing.

    Returning None does NOT make an unknown `--phylum` harmless: `_cmd_detect`
    rejects a phylum that matches no curated family regardless of whether
    argparse was able to constrain the choices.
    """
    try:
        return available_phyla(MatpredictConfig.from_env(repo_root=Path.cwd()).db_root) or None
    except Exception as err:  # noqa: BLE001 - parser construction must never fail here
        logger.warning(
            "could not read the phylum list from the database root (%s) -- "
            "`matpredict detect --phylum` will not validate its argument against "
            "the curated phyla", err,
        )
        return None


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
    detect.add_argument(
        "--phylum",
        required=False,
        choices=_phylum_choices(),
        help=(
            "Restrict detection to one phylum's curated families outright, skipping "
            "taxid-based routing entirely. Use when the genome's phylum is known but "
            "its taxid routes badly -- without it, a taxid no family's taxonomic_scope "
            "covers falls back to searching every family in every phylum. The choices "
            "are read from the database root at startup, so they always match what is "
            "actually curated."
        ),
    )
    detect.add_argument(
        "--emit-cds-fasta",
        action="store_true",
        help=(
            "Also emit real CDS features with translation= attributes in the GFF3, and "
            "write a companion detected_loci.fasta next to it. The companion FASTA "
            "contains the FULL sequence of every contig the detection results "
            "reference, not just the locus spans, so it can be large (a whole "
            "chromosome-scale contig per referenced contig). Off by default."
        ),
    )
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
