"""argparse wiring for the `matpredict detect` subcommand."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from MATPredict import logger
from MATPredict.config import MatpredictConfig
from MATPredict.detect.benchmark import run_benchmark
from MATPredict.detect.confidence_cap import CONFIDENCE_CAP_POLICIES, DEFAULT_CONFIDENCE_CAP_POLICY
from MATPredict.detect.strong_single_gene import (
    DEFAULT_STRONG_SINGLE_GENE_POLICY, STRONG_SINGLE_GENE_POLICIES,
)
from MATPredict.detect.homothallic import DEFAULT_HOMOTHALLIC_POLICY, HOMOTHALLIC_POLICIES
from MATPredict.detect.partial_core import DEFAULT_PARTIAL_CORE_POLICY, PARTIAL_CORE_POLICIES
from MATPredict.detect.supported_span import DEFAULT_SUPPORTED_MIN_BITSCORE
from MATPredict.detect.pipeline import (
    DEFAULT_MAX_POLISHED_CLUSTERS_PER_FAMILY, DEFAULT_POLISH_CORE_PAIR, DEFAULT_POLISH_STRONG_IDENTITY, DetectionOutcome,
    EvidenceFloor,
    not_searched_reason, run_pipeline,
)
from MATPredict.detect.family_registry import (
    DEFAULT_SCOPED_ROSTER_POLICY, SCOPED_ROSTER_POLICIES, available_phyla, load_all_families, route,
)
from MATPredict.detect.genbank_input import is_genbank, prepare_genbank_input
from MATPredict.detect.reference_fasta import build_reference_fasta
from MATPredict.detect.provenance import RunClock, code_stamp, run_block
from MATPredict.detect.report import _taxonomy_source, write_detection_gff3, write_detection_report
from MATPredict.detect.rollout_aggregate import aggregate_reports, write_rollout_summary
from MATPredict.detect.scope_audit import audit_gene_scopes, audit_scope, record_taxids_by_family
from MATPredict.detect.suppress import default_suppress_paths, filter_rows, load_suppress_list

_ROLLOUT_REPORT_FILENAME = "detection_report.yaml"


_COMPRESSION_MAGIC = ((b"\x1f\x8b", "gzip"), (b"\x28\xb5\x2f\xfd", "zstd"))


def _refuse_compressed_genome(genome: str) -> None:
    """`detect` reads plain FASTA. A compressed file would otherwise surface minutes later as a BLAST or
    index error, so say what is wrong and what to do. Recognised by the first bytes, not the file name."""
    path = Path(genome)
    try:
        with open(path, "rb") as fh:
            head = fh.read(4)
    except OSError:
        return  # missing or unreadable: the later stages report that
    for magic, kind in _COMPRESSION_MAGIC:
        if head.startswith(magic):
            tool = "zcat" if kind == "gzip" else "zstd -dc"
            raise ValueError(
                f"{path.name} is {kind}-compressed; `matpredict detect` reads plain FASTA. "
                f"Decompress it first, for example `{tool} {path} > genome.fna` "
                "(into $SCRATCH on a SLURM node), then pass that file as --genome.")


def _cmd_detect(args: argparse.Namespace) -> int:
    if not args.genome or not args.out_dir:
        raise ValueError("--genome and --out-dir are required for `matpredict detect`")
    _refuse_compressed_genome(args.genome)
    config = MatpredictConfig.from_env(repo_root=Path.cwd())
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    clock = RunClock()
    code_stamp()    # read the checkout's identity now, so the report names the code that was there when the run began
    _use_genbank_genome(args, out_dir)

    # Route FIRST, then build the reference FASTA from the routed families
    # only. The order matters: the reference FASTA is the tblastn query set, so
    # building it before routing (as this did) means every run pays to align
    # every phylum's curated proteins no matter how narrowly it routed. The
    # same `RoutingDecision` object is handed to `run_pipeline` so it cannot
    # re-route to a different set than the query set was built for, and so the
    # taxonomy lookup happens once per run.
    routing = route(args.taxid, load_all_families(config.db_root), phylum=args.phylum,
                    exhaustive=getattr(args, "exhaustive", False),
                    scoped_roster_policy=_scoped_roster_policy_from_args(args))
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
        write_detection_report(outcome, out_dir / "detection_report.yaml",
                               run=_run_doc(args, config, routing, clock))
        _write_html_report(args, out_dir)
        print(f"not searched -> {out_dir}: {outcome.not_searched_reason}")
        return 0
    exclude_record_ids = frozenset(
        r.strip() for r in (args.exclude_records or "").split(",") if r.strip()
    )
    reference_fasta = build_reference_fasta(
        config.db_root, out_dir / "_reference.faa",
        family_keys={f.key for f in routing.families},
        exclude_record_ids=exclude_record_ids,
        families=routing.families,
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
        polish_strong_identity=_polish_strong_identity_from_args(args),
        polish_core_pair=getattr(args, "polish_core_pair", DEFAULT_POLISH_CORE_PAIR),
        partial_core_policy=_partial_core_policy_from_args(args),
        homothallic_policy=_homothallic_policy_from_args(args),
        confidence_cap_policy=_confidence_cap_policy_from_args(args),
        strong_single_gene_policy=_strong_single_gene_policy_from_args(args),
        supported_min_bitscore=getattr(args, "supported_min_bitscore", DEFAULT_SUPPORTED_MIN_BITSCORE),
    )

    # `genome_fasta` is passed ONLY when asked for: it is what makes
    # `write_detection_gff3` additionally emit CDS features with
    # `translation=` attributes and a companion FASTA, and the companion
    # FASTA holds every referenced contig's FULL sequence, which is large.
    gff3_kwargs = {"genome_fasta": Path(args.genome)} if args.emit_cds_fasta else {}
    write_detection_gff3(outcome, out_dir / "detected_loci.gff3", **gff3_kwargs)
    write_detection_report(outcome, out_dir / "detection_report.yaml",
                           run=_run_doc(args, config, routing, clock))
    _write_tracks(args, out_dir)
    _write_html_report(args, out_dir)
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


def _use_genbank_genome(args, out_dir: Path) -> None:
    """When --genome is a GenBank file, run on what it holds: its sequence is the genome, its CDS proteins
    are the gene calls (the proteome route, unless --proteins was given), and the file is the annotation
    track in the report. A GenBank file given as --annotation is left as it is."""
    if not is_genbank(args.genome):
        return
    gb = prepare_genbank_input(args.genome, out_dir / "_genbank")
    print(f"GenBank input {gb.source.name}: {gb.converted.contigs} contig(s), "
          f"{gb.converted.proteins} CDS protein(s) -> {out_dir / '_genbank'}")
    args.annotation = getattr(args, "annotation", None) or str(gb.source)
    args.genome = str(gb.converted.genome_fasta)
    if not getattr(args, "proteins", None):
        args.proteins = str(gb.converted.proteome_fasta)


def _helixerlite_command(args) -> list[str] | None:
    exe = getattr(args, "helixerlite", None) or os.environ.get("MATPREDICT_HELIXERLITE")
    return [exe] if exe else None


def _write_tracks(args, out_dir: Path) -> None:
    """tracks.yaml: helixer-lite on the locus windows and/or the GenBank annotation, compared with the
    homology models. Report only. Never fails the run, for the same reason the HTML report cannot."""
    helixer = _helixerlite_command(args)
    annotation = getattr(args, "annotation", None)
    if helixer is None and not annotation:
        return
    from MATPredict.detect.track_report import TRACKS_YAML, build_tracks_doc, write_tracks
    from MATPredict.report.cli import load_report
    try:
        doc, _ = load_report(out_dir)
        tracks = build_tracks_doc(
            doc, Path(args.genome), helixer_cmd=helixer,
            helixer_lineage=getattr(args, "helixer_lineage", "fungi"),
            helixer_cpus=getattr(args, "helixer_cpus", 4),
            genbank=Path(annotation) if annotation else None)
        if tracks is not None:
            write_tracks(out_dir / TRACKS_YAML, tracks)
            print(f"tracks -> {out_dir / TRACKS_YAML}")
    except Exception as exc:  # noqa: BLE001 - tracks are report-only; the detection result stands
        logger.warning("tracks not written (%s: %s); the detection results are complete",
                       type(exc).__name__, exc)


def html_report_enabled(args) -> bool:
    """`--html`/`--no-html` when given; otherwise on, unless `$MATPREDICT_HTML` is 0/false/no/off.

    The variable exists for batch scripts: they may run `detect` from an older
    frozen worktree that would reject an unknown `--no-html` flag, while an
    unknown environment variable is simply ignored."""
    flag = getattr(args, "html", None)
    if flag is not None:
        return bool(flag)
    return os.environ.get("MATPREDICT_HTML", "1").strip().lower() not in ("0", "false", "no", "off")


def _write_html_report(args, out_dir: Path) -> None:
    """The readable report next to the YAML. Never fails the run: the YAML and GFF3 are the results; a report
    error (or a missing PDF engine) is logged and the run still exits 0."""
    if not (html_report_enabled(args) or getattr(args, "pdf", False)):
        return
    from MATPredict.report.cli import write_run_report
    try:
        for path in write_run_report(out_dir):
            print(f"report -> {path}")
    except Exception as exc:  # noqa: BLE001 - a report must never cost the detection result
        logger.warning("report not written (%s: %s); the detection results are complete",
                       type(exc).__name__, exc)
        return
    if getattr(args, "pdf", False):
        try:
            for path in write_run_report(out_dir, pdf=True)[1:]:
                print(f"pdf -> {path}")
        except Exception as exc:  # noqa: BLE001
            logger.warning("PDF not written (%s: %s); report.html is complete", type(exc).__name__, exc)


def _run_doc(args, config, routing, clock: RunClock) -> dict:
    """The report's `run` block (`detect.provenance`) for this invocation."""
    return run_block(
        clock=clock,
        genome=Path(args.genome) if args.genome else None,
        proteins=Path(args.proteins) if args.proteins else None,
        db_root=config.db_root,
        taxonomy_source=_taxonomy_source(),
        sample=getattr(args, "sample", None),
        organism=getattr(args, "organism", None),
        taxid=args.taxid,
        phylum=getattr(routing, "phylum", None),
        parameters={
            "phylum": args.phylum,
            "exhaustive": getattr(args, "exhaustive", False),
            "genetic_code": _resolve_genetic_code(args),
            "exclude_records": args.exclude_records or None,
            "annotation": Path(args.annotation).name if getattr(args, "annotation", None) else None,
            "min_hits": args.min_hits,
            "min_identity": args.min_identity,
            "require_core_role": args.require_core_role,
            "max_polished_clusters_per_family": _polish_cap_from_args(args),
            "polish_strong_identity": _polish_strong_identity_from_args(args),
            "polish_core_pair": getattr(args, "polish_core_pair", DEFAULT_POLISH_CORE_PAIR),
            "partial_core_policy": _partial_core_policy_from_args(args),
            "homothallic_policy": _homothallic_policy_from_args(args),
            "confidence_cap_policy": _confidence_cap_policy_from_args(args),
            "strong_single_gene_policy": _strong_single_gene_policy_from_args(args),
            "scoped_roster_policy": _scoped_roster_policy_from_args(args),
            "scoped_roster": getattr(routing, "scoped_roster", None),
            "supported_min_bitscore": getattr(args, "supported_min_bitscore", DEFAULT_SUPPORTED_MIN_BITSCORE),
        },
    )


def _polish_cap_from_args(args) -> int | None:
    """The per-family polish cap for this run: the flag's value, 0 meaning no
    cap. An args object without the attribute (a test double) gets the
    pipeline's default, never an implicit "off"."""
    cap = getattr(args, "max_polished_clusters_per_family",
                  DEFAULT_MAX_POLISHED_CLUSTERS_PER_FAMILY)
    return None if cap is None or cap <= 0 else cap


def _partial_core_policy_from_args(args) -> str:
    """The partial-coverage vote policy for this run; an args object without the
    attribute (a test double) gets the default, `off`."""
    return getattr(args, "partial_core_policy", DEFAULT_PARTIAL_CORE_POLICY)


def _homothallic_policy_from_args(args) -> str:
    """The both-genes homothallic policy for this run; an args object without the attribute
    (a test double) gets the default, `off`."""
    return getattr(args, "homothallic_policy", DEFAULT_HOMOTHALLIC_POLICY)


def _confidence_cap_policy_from_args(args) -> str:
    """The confidence-cap policy for this run; an args object without the attribute
    (a test double) gets the default, `off`."""
    return getattr(args, "confidence_cap_policy", DEFAULT_CONFIDENCE_CAP_POLICY)
def _strong_single_gene_policy_from_args(args) -> str:
    """The strong-single-gene policy for this run; an args object without the
    attribute (a test double) gets the default, `off`."""
    return getattr(args, "strong_single_gene_policy", DEFAULT_STRONG_SINGLE_GENE_POLICY)


def _scoped_roster_policy_from_args(args) -> str:
    """The taxon-scoped roster policy for this run; an args object without the
    attribute (a test double) gets the default, `off`."""
    return getattr(args, "scoped_roster_policy", DEFAULT_SCOPED_ROSTER_POLICY)


def _polish_strong_identity_from_args(args) -> float | None:
    """The identity tier above the polish cap for this run: the flag's value as a
    percentage, 0 meaning no tier. An args object without the attribute (a test
    double) gets the pipeline's default."""
    value = getattr(args, "polish_strong_identity", DEFAULT_POLISH_STRONG_IDENTITY)
    return None if value is None or value <= 0 else float(value)


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
    for record_id, gene, taxid, scope in audit_gene_scopes(db_root, families):
        print(f"{record_id}: carries scoped gene {gene} but its taxid {taxid} is outside the gene's "
              f"taxonomic_scope {scope} -- the protein is never searched")
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
    detect.add_argument("--html", action=argparse.BooleanOptionalAction, default=None,
                        help="Write report.html in --out-dir (default: on; $MATPREDICT_HTML=0 turns the default off)")
    detect.add_argument("--pdf", action="store_true",
                        help="Also write report.pdf (needs WeasyPrint or Chrome/Chromium; off by default)")
    detect.add_argument("--sample", required=False,
                        help="Sample name for the report (default: the genome file name up to the first dot)")
    detect.add_argument("--organism", required=False,
                        help="Organism name for the report, e.g. 'Fusarium oxysporum f. sp. lactucae'. "
                             "Shown as given; routing still uses --taxid")
    detect.add_argument("--out-dir", required=False)
    detect.add_argument("--annotation", required=False,
                        help="GenBank annotation of the genome, shown as a gene-model track in the report. "
                             "A GenBank file given as --genome is used the same way, and its CDS proteins "
                             "also serve as the gene calls.")
    detect.add_argument("--helixerlite", required=False,
                        help="helixer-lite executable (default: $MATPREDICT_HELIXERLITE). When set, it runs "
                             "on the locus windows only and its models are a report track; calls are unchanged.")
    detect.add_argument("--helixer-lineage", default="fungi", help="helixer-lite lineage (default: fungi)")
    detect.add_argument("--helixer-cpus", type=int, default=4, help="CPUs for helixer-lite (default: 4)")
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
        "--supported-min-bitscore", type=float, default=DEFAULT_SUPPORTED_MIN_BITSCORE,
        help="Bitscore floor for the report-only `supported_span`: a locus's own hits at or above it, plus its "
             "modelled genes, define the supported extent; weaker hits do not. Does not change calls or the cluster "
             f"span. Default: {DEFAULT_SUPPORTED_MIN_BITSCORE:g}.",
    )
    detect.add_argument(
        "--polish-strong-identity", type=float, default=DEFAULT_POLISH_STRONG_IDENTITY,
        help="Also polish every admitted cluster whose best identity is at or above "
             "this percentage, even past the cap; the rest of the cap is filled from "
             "the usual rank. Stops a true locus that hits few genes at high identity "
             "being outranked by noise that hits more genes at low identity. "
             f"Default: {DEFAULT_POLISH_STRONG_IDENTITY:g}. 0 = the plain cap.",
    )
    detect.add_argument(
        "--polish-core-pair", action=argparse.BooleanOptionalAction, default=DEFAULT_POLISH_CORE_PAIR,
        help="Also polish every admitted cluster with two distinct live core genes of one idiomorph "
             "(e.g. MTLA1 + MTLA2) whatever its rank, on top of the cap (it takes no slot). The property is the cluster's own, so adding "
             "families that admit more clusters cannot push such a cluster out of the cap. "
             "Default: on. --no-polish-core-pair restores the earlier behaviour.",
    )
    detect.add_argument(
        "--partial-core-policy", choices=PARTIAL_CORE_POLICIES, default=DEFAULT_PARTIAL_CORE_POLICY,
        help="off: the best-bitscore idiomorph vote stands (default). coverage_profile: a core hit "
             "that covers < 80%% of its reference does not beat a full-length (>= 90%%) core gene of the "
             "other idiomorph when the alpha-box/HMG-box profiles agree; where they disagree or are "
             "unavailable the call is `undetermined` (see `idiomorph_vote_note`). A partial hit near "
             "a contig end is treated as truncation and the vote is kept. Independent of --idiomorph-policy.",
    )
    detect.add_argument(
        "--homothallic-policy", choices=HOMOTHALLIC_POLICIES, default=DEFAULT_HOMOTHALLIC_POLICY,
        help="off: a locus with an alpha-box and an HMG-box gene is reported as one idiomorph unless the "
             "reference-coverage rule of `classify_locus` passes (default). both_genes: a locus whose modelled "
             "core genes include a COMPLETE alpha-box and a COMPLETE HMG-box protein, each confirmed by the "
             "gene-identity profiles (not an HMG paralog or MAT1-1-3), within the family's separation limit, "
             "is classed `homothallic_candidate` (see `homothallic_note`). A candidate, never an assertion; "
             "the idiomorph is reported `undetermined` (the vote's winner is kept in the note as `vote_idiomorph`) "
             "and the confidence is unchanged.",
    )
    detect.add_argument(
        "--confidence-cap-policy", choices=CONFIDENCE_CAP_POLICIES, default=DEFAULT_CONFIDENCE_CAP_POLICY,
        help="off: no extra cap (default). partial_weak: a `high` call drops to `medium` when its best core "
             "hit is a contig-edge partial (< 80%% of its reference, within 500 bp of a contig end, no "
             "full-length core gene) or a weak hit (< 45%% identity) with no flank gene. weak_only: the second "
             "criterion alone. Changes confidence only; the reason is written as `confidence_cap`.",
    )
    detect.add_argument(
        "--strong-single-gene-policy", choices=STRONG_SINGLE_GENE_POLICIES,
        default=DEFAULT_STRONG_SINGLE_GENE_POLICY,
        help="off: a lone core gene with no flank and no second modelled gene is withheld (default). "
             "on: in a family with no other report, one modelled core gene at >= 85%% identity covering "
             ">= 90%% of its reference is reported as `idiomorph_gene_only`, confidence `low`, flagged "
             "`strong_single_gene`; two such loci of opposite idiomorph are both reported and the genome "
             "gets a `two_idiomorphs` statement. Independent of --idiomorph-policy and --partial-core-policy.",
    )
    detect.add_argument(
        "--scoped-roster-policy", choices=SCOPED_ROSTER_POLICIES, default=DEFAULT_SCOPED_ROSTER_POLICY,
        help="off: a roster gene that declares its own `taxonomic_scope` in order.yml is not searched "
             "(default; the run is the one it was before such genes were added). on: it is searched, and "
             "counts, for a genome whose lineage reaches its scope, and nowhere else (prototype: GAP1/BGL1 "
             "as optional bonus flanks of Ascomycota MAT in Pleosporaceae).",
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
