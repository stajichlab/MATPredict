"""argparse wiring for the `matpredict detect` subcommand."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from MATPredict import logger
from MATPredict.config import MatpredictConfig
from MATPredict.detect.benchmark import run_benchmark
from MATPredict.detect.pipeline import (
    DEFAULT_MAX_POLISHED_CLUSTERS_PER_FAMILY, DetectionOutcome, EvidenceFloor,
    not_searched_reason, run_pipeline,
)
from MATPredict.detect.family_registry import available_phyla, load_all_families, route
from MATPredict.detect.reference_fasta import build_reference_fasta
from MATPredict.detect.report import write_detection_gff3, write_detection_report
from MATPredict.detect.rollout_aggregate import aggregate_reports, write_rollout_summary
from MATPredict.detect.scope_audit import audit_scope, record_taxids_by_family
from MATPredict.detect.suppress import default_suppress_paths, filter_rows, load_suppress_list

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
    routing = route(args.taxid, load_all_families(config.db_root), phylum=args.phylum,
                    exhaustive=getattr(args, "exhaustive", False))
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
    # Withheld records leave BOTH the query set and every db-derived quantity;
    # see `detect.holdout`. Passed to the two calls from one parsed set so the
    # reference FASTA and the pipeline can never disagree about what is hidden,
    # which would silently turn a recall measurement back into a
    # self-consistency one.
    if routing.routing_mode == "not_searched":
        # Curator's ruling 2026-09-26: no curated family covers this genome,
        # so nothing is searched and no query set is built. The same outputs
        # are still written, so a panel's per-genome directory is uniform and
        # the report says WHY it holds no call.
        outcome = DetectionOutcome(
            results=[], routing_mode=routing.routing_mode,
            routing_error=routing.routing_error,
            not_searched_reason=not_searched_reason(routing),
        )
        write_detection_gff3(outcome, out_dir / "detected_loci.gff3")
        write_detection_report(outcome, out_dir / "detection_report.yaml")
        print(f"not searched -> {out_dir}: {outcome.not_searched_reason}")
        return 0
    exclude_record_ids = frozenset(
        r.strip() for r in (args.exclude_records or "").split(",") if r.strip()
    )
    reference_fasta = build_reference_fasta(
        config.db_root, out_dir / "_reference.faa",
        family_keys={f.key for f in routing.families},
        exclude_record_ids=exclude_record_ids,
    )
    evidence_floor = EvidenceFloor(
        min_hits=args.min_hits, min_identity=args.min_identity,
        require_core_role=args.require_core_role,
    )
    outcome = run_pipeline(
        genetic_code=_resolve_genetic_code(args),
        genome_fasta=Path(args.genome),
        proteome_fasta=Path(args.proteins) if args.proteins else None,
        taxid=args.taxid,
        db_root=config.db_root,
        reference_fasta=reference_fasta,
        evidence_floor=evidence_floor,
        evidence_diagnostics_path=Path(args.evidence_diagnostics) if args.evidence_diagnostics else None,
        routing=routing,
        exclude_record_ids=exclude_record_ids,
        max_polished_clusters_per_family=_polish_cap_from_args(args),
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


def _polish_cap_from_args(args) -> int | None:
    """The per-family polish cap for this run: the flag's value, 0 meaning no
    cap. An args object without the attribute (a test double) gets the
    pipeline's default, never an implicit "off"."""
    cap = getattr(args, "max_polished_clusters_per_family",
                  DEFAULT_MAX_POLISHED_CLUSTERS_PER_FAMILY)
    return None if cap is None or cap <= 0 else cap


def _cmd_suppress_filter(args: argparse.Namespace) -> int:
    """Filter a panel list (`ASMID<TAB>...` per line) against the suppress
    lists: kept rows to stdout, the skipped count (and each skipped ASMID) to
    stderr. The runners pipe their panel list through this before fanning out."""
    paths = [] if args.no_default_lists else default_suppress_paths(
        MatpredictConfig.from_env(repo_root=Path.cwd()).db_root)
    paths += [Path(p) for p in args.suppress or []]
    suppressed = load_suppress_list(paths)
    source = Path(args.list).read_text() if args.list else sys.stdin.read()
    rows = [r for r in source.splitlines() if r.strip()]
    kept, skipped = filter_rows(rows, suppressed)
    for row in kept:
        print(row)
    for row in skipped:
        print(f"suppressed: {row.split(chr(9))[0]}", file=sys.stderr)
    print(f"skipped {len(skipped)} suppressed genome(s) of {len(rows)} "
          f"(lists: {', '.join(str(p) for p in paths) or 'none'})", file=sys.stderr)
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


def _resolve_genetic_code(args) -> int | None:
    """The translation table for this run, or None to let the pipeline derive it.

    Precedence, matching how `--phylum` overrides taxid routing: an explicit
    `--genetic-code` wins; then `--genetic-code-map` keyed on the genome
    FASTA's basename; then None, which the pipeline resolves from the taxid and
    finally falls back to table 1. Never guessed.

    The map exists because the code is a per-GENOME property, not a per-run
    one: the local BFD sample sheet carries 22,528 genomes at table 1 and
    1,154 at table 12, and a batch spanning both cannot be described by a
    single flag.
    """
    explicit = getattr(args, "genetic_code", None)
    if explicit is not None:
        return int(explicit)
    path = getattr(args, "genetic_code_map", None)
    if not path:
        return None
    wanted = Path(args.genome).name
    for suffix in (".gz", ".fna", ".fa", ".fasta"):
        if wanted.endswith(suffix):
            wanted = wanted[: -len(suffix)]
    for line in Path(path).read_text().splitlines():
        parts = [c.strip() for c in line.replace(",", "\t").split("\t") if c.strip()]
        if len(parts) < 2:
            continue
        key = parts[0]
        for suffix in (".gz", ".fna", ".fa", ".fasta"):
            if key.endswith(suffix):
                key = key[: -len(suffix)]
        if key == wanted:
            try:
                return int(parts[1])
            except ValueError:
                return None  # a header row, or a malformed value: derive instead
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
    detect.add_argument(
        "--exclude-records", default="",
        help=("comma-separated curated record_ids to withhold from this run, for "
              "leave-one-out recall. They leave the query set AND every db-derived "
              "quantity, so the genome a withheld record came from becomes a genuine "
              "held-out test rather than a self-consistency check. See "
              "MATPredict.detect.holdout for choosing a set by taxonomic radius."))
    detect.add_argument("--taxid", required=False, type=int)
    detect.add_argument("--out-dir", required=False)
    detect.add_argument("--evidence-diagnostics", required=False)
    # Defaults are read from `EvidenceFloor` itself, never restated here: a
    # hard-coded CLI default would silently re-impose the old permissive floor
    # on every flagless run the moment the two drifted apart.
    _floor_defaults = EvidenceFloor()
    detect.add_argument("--min-hits", type=int, default=_floor_defaults.min_hits,
                        help="Minimum DISTINCT genes a family needs in a cluster "
                             f"to reach Stage 2 polishing (default: {_floor_defaults.min_hits})")
    detect.add_argument("--min-identity", type=float, default=_floor_defaults.min_identity,
                        help="Minimum best-hit percent identity for a family to reach "
                             "Stage 2 polishing (default: no identity cutoff)")
    # BooleanOptionalAction, not store_true: the default is now True, so an
    # explicit `--no-require-core-role` off switch has to exist. The positive
    # `--require-core-role` spelling keeps working unchanged.
    detect.add_argument("--require-core-role", action=argparse.BooleanOptionalAction,
                        default=_floor_defaults.require_core_role,
                        help="Require at least one core_MAT hit for a family to reach "
                             f"Stage 2 polishing (default: {_floor_defaults.require_core_role})")
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
        "--genetic-code", type=int, default=None,
        help="NCBI translation table for this genome (e.g. 12 for the CUG-Ser1 "
             "clade). Omit to derive it from --taxid, which reads the same "
             "cached taxonomy document routing already fetches; falls back to "
             "1 when that lookup cannot answer. An explicit value always wins.",
    )
    detect.add_argument(
        "--genetic-code-map", default=None,
        help="TSV/CSV file mapping a genome id to its translation table, for "
             "batch runs where each genome differs. Two columns, id and code, "
             "with an optional header; the id is matched against the genome "
             "FASTA's basename with extensions stripped. Used only when "
             "--genetic-code is not given. The BFD sample sheet's "
             "ASMID/TRANSL_TABLE columns are exactly this shape.",
    )
    detect.add_argument(
        "--max-polished-clusters-per-family", type=int,
        default=DEFAULT_MAX_POLISHED_CLUSTERS_PER_FAMILY,
        help="Polish at most N admitted clusters per family, ranked before "
             "polishing by distinct genes, then best identity, then hit count. "
             f"Default: {DEFAULT_MAX_POLISHED_CLUSTERS_PER_FAMILY} (curator's ruling "
             "2026-09-26; measured -45%% compute, 4 of 613 calls lost). 0 = no cap.",
    )
    detect.add_argument(
        "--exhaustive", action="store_true",
        help="Search every family in every phylum when no curated family covers "
             "the genome (no taxid, an unresolvable lineage, or an uncurated "
             "phylum). Default: such a genome is NOT searched and its report says "
             "so (curator's ruling 2026-09-26; the chytrid control ran ~31 min "
             "per genome under exhaustive and called nothing). Never widens a "
             "route that matched.",
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

    suppress_parser = action.add_parser(
        "suppress-filter",
        help="Drop suppressed genomes from a panel list (ASMID<TAB>... per line)",
    )
    suppress_parser.add_argument("--list", default=None,
                                 help="panel list to filter; default: stdin")
    suppress_parser.add_argument(
        "--suppress", action="append", default=[],
        help="an extra suppress list (repeatable), added to the defaults: the "
             "BFD list and <db root>/suppress.txt")
    suppress_parser.add_argument("--no-default-lists", action="store_true",
                                 help="use only the --suppress lists given")
    suppress_parser.set_defaults(func=_cmd_suppress_filter)

    audit_scope_parser = action.add_parser(
        "audit-scope",
        help="Audit every curated record's taxid against its own family's taxonomic_scope",
    )
    audit_scope_parser.add_argument("--db-root", default="db")
    audit_scope_parser.set_defaults(func=_cmd_audit_scope)
