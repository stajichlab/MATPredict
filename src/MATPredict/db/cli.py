"""argparse wiring for the `matpredict curate-db` subcommand group."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests
import yaml
from Bio import SeqIO

from MATPredict.config import MatpredictConfig
from MATPredict.db import gff_export
from MATPredict.db.build_duckdb import build as build_duckdb
from MATPredict.db.curate import accept_candidate, propose_candidate, reject_candidate
from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.schema import validate_gene_vocabulary
from MATPredict.db.uniprot_client import UniprotClient
from MATPredict.db.validate import _client_for, _independent_translation, validate_record


def _config(args: argparse.Namespace) -> MatpredictConfig:
    return MatpredictConfig.from_env(repo_root=Path.cwd())


def _http_transport(url: str, max_attempts: int = 4, backoff_seconds: float = 2.0) -> str:
    """GET a URL, retrying on rate-limit/server errors so a transient failure never gets
    written to CachedFetcher's on-disk cache (only a value this function *returns* is
    cached, and it only returns after a genuinely successful response)."""
    last_error: Exception | None = None
    for attempt in range(max_attempts):
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        last_error = requests.HTTPError(f"{response.status_code} for {url}: {response.text[:200]}")
        if response.status_code in (429, 500, 502, 503, 504) and attempt < max_attempts - 1:
            time.sleep(backoff_seconds * (attempt + 1))
            continue
        break
    raise last_error


def _make_clients(config: MatpredictConfig) -> tuple[NcbiClient, UniprotClient]:
    """Build the NCBI/UniProt client pair used by both `validate` and `build-gff`."""
    fetcher = CachedFetcher(cache_dir=config.cache_dir, transport=_http_transport)
    ncbi = NcbiClient(email=config.ncbi_email, api_key=config.ncbi_api_key, fetcher=fetcher)
    uniprot = UniprotClient(fetcher=fetcher)
    return ncbi, uniprot


def _cmd_propose(args: argparse.Namespace) -> int:
    config = _config(args)
    record = yaml.safe_load(Path(args.record_file).read_text())
    candidate_dir = propose_candidate(config.db_root, args.phylum, record)
    print(f"proposed candidate at {candidate_dir}")
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    config = _config(args)
    record_path = config.db_root / "candidates" / args.phylum / args.record_id / "metadata.yaml"
    record = yaml.safe_load(record_path.read_text())

    ncbi, uniprot = _make_clients(config)

    result = validate_record(record, ncbi=ncbi, uniprot=uniprot)
    record["validation"].update(result)

    # A curated gene whose name is not declared by its locus in order.yml is silently
    # dropped by detect.search._attribute, so it can never be found in any genome. Report
    # it here rather than letting it pass review unnoticed. The errors are not written
    # into the record (they describe order.yml, not the record's own data).
    order_path = config.db_root / args.phylum / "order.yml"
    vocabulary_errors: list[str] = []
    if not order_path.exists():
        # Nothing to reconcile against. Say so rather than pass silently: the
        # database-wide guard in tests/db/test_schema.py is the permanent net.
        print(f"gene-vocabulary check skipped: no {order_path}")
    elif not (record.get("mating_type") or {}).get("locus_name"):
        print("gene-vocabulary check skipped: record declares no mating_type.locus_name")
    else:
        order_doc = yaml.safe_load(order_path.read_text())
        vocabulary_errors = validate_gene_vocabulary(record, order_doc)
        for error in vocabulary_errors:
            print(f"gene-vocabulary error: {error}")

    # A re-validation that reveals a problem must downgrade the record's status so it gets
    # re-reviewed, unless it was explicitly rejected already (never silently un-reject).
    check_failed = (
        result.get("accession_resolved") is False
        or (result.get("sequence_match") or {}).get("status") == "fail"
        or result.get("taxonomy_current") is False
        or bool(vocabulary_errors)
    )
    if check_failed and record["validation"].get("status") != "rejected":
        record["validation"]["status"] = "needs_review"

    record_path.write_text(yaml.safe_dump(record, sort_keys=False))
    print(json.dumps(result, indent=2, default=str))
    return 0


def _cmd_accept(args: argparse.Namespace) -> int:
    config = _config(args)
    accepted_dir = accept_candidate(config.db_root, args.phylum, args.order_or_family, args.record_id)
    print(f"accepted {args.record_id} -> {accepted_dir}")
    return 0


def _cmd_reject(args: argparse.Namespace) -> int:
    config = _config(args)
    reject_candidate(config.db_root, args.phylum, args.record_id, reason=args.reason)
    print(f"rejected {args.record_id}: {args.reason}")
    return 0


def find_records_missing_proteins_faa(db_root: Path) -> list[tuple[str, str, str]]:
    """(phylum, order_or_family, record_id) for every accepted (non-candidate)
    record whose proteins.faa does not exist on disk yet, OR exists but is
    empty (zero FASTA entries) -- an empty file is what build_gff_for_record
    correctly writes for a record whose genes all lack a protein_accession,
    but it must not be mistaken for "already handled": when that record is
    eventually curated with real accessions, backfill-gff must pick it up
    again rather than silently skipping it because a file happens to exist."""
    missing = []
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        parts = meta_path.relative_to(db_root).parts
        if parts[0] == "candidates":
            continue
        record_dir = meta_path.parent
        proteins_path = record_dir / "proteins.faa"
        if not proteins_path.exists() or ">" not in proteins_path.read_text():
            missing.append((parts[0], parts[1], parts[2]))
    return missing


def find_records_with_stale_locus_gbk(db_root: Path) -> list[tuple[str, str, str]]:
    """(phylum, order_or_family, record_id) for every accepted (non-candidate)
    record whose locus.gbk predates the write_genbank upgrade that added real
    nucleotide sequence and real CDS/translation features.

    Staleness test: the record has at least one `present` gene AND its locus.gbk
    is missing, or parses to zero `CDS` features across all its SeqRecords. The
    old generator emitted only `gene` features, so a zero-CDS file is exactly the
    pre-upgrade output; a regenerated file for a record with a present gene whose
    sequence could be obtained carries at least one `CDS`.

    The `present` gene precondition matters: a record whose genes are all
    `present: false` would still produce a zero-CDS locus.gbk after a perfectly
    successful regeneration, so selecting it would make the sweep report it as
    stale forever with nothing to fix.

    Parsing is done with Bio.SeqIO, not a substring search for "     CDS  ":
    GenBank feature indentation is column-sensitive and a text match on it would
    silently mis-classify records if the writer's spacing ever changed. A file
    Bio.SeqIO cannot parse counts as STALE rather than raising -- an unreadable
    locus.gbk is not evidence of freshness, and regenerating it is the fix.
    """
    stale = []
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        parts = meta_path.relative_to(db_root).parts
        if parts[0] == "candidates":
            continue
        record = yaml.safe_load(meta_path.read_text()) or {}
        if not any(gene.get("present", True) for gene in record.get("genes", [])):
            continue
        gbk_path = meta_path.parent / "locus.gbk"
        if not gbk_path.exists():
            stale.append((parts[0], parts[1], parts[2]))
            continue
        try:
            cds_count = sum(
                1
                for seq_record in SeqIO.parse(gbk_path, "genbank")
                for feature in seq_record.features
                if feature.type == "CDS"
            )
        except Exception:  # noqa: BLE001 -- an unparseable file is stale, not fatal
            cds_count = 0
        if cds_count == 0:
            stale.append((parts[0], parts[1], parts[2]))
    return stale


def placeholder_segments(gbk_path: Path) -> list[str]:
    """The ids of every SeqRecord in this locus.gbk whose sequence is entirely "N",
    i.e. a segment whose real nucleotide sequence could not be fetched and that
    write_genbank therefore left as its all-"N" placeholder.

    Reported per record after a backfill so the sweep's own output names exactly
    which records (and which segments within them) still lack real sequence,
    instead of deferring that to a separate manual audit of the written files.
    Returns [] rather than raising on an unreadable file: this is a reporting aid
    and must never turn a succeeded record into a failure."""
    placeholders = []
    try:
        for seq_record in SeqIO.parse(gbk_path, "genbank"):
            sequence = str(seq_record.seq)
            if sequence and set(sequence.upper()) == {"N"}:
                placeholders.append(seq_record.id)
    except Exception:  # noqa: BLE001 -- a reporting aid must never fail the sweep
        return []
    return placeholders


def build_gff_for_record(
    db_root: Path, phylum: str, order_or_family: str, record_id: str,
    ncbi: NcbiClient, uniprot: UniprotClient,
) -> None:
    """Fetch every present gene's protein sequence and write
    locus.gff3/locus.gbk/proteins.faa for one accepted record. A gene with a
    real `protein_accession` is fetched directly; a gene with none but real
    curated genomic coordinates (an unannotated MAG assembly with no NCBI
    protein record to cite) has its sequence independently derived via
    `_independent_translation`, which fetches and translates from the
    record's own recorded coordinates -- the same function `validate.py`
    already uses to cross-check a claimed accession's sequence, reused here
    as the sequence SOURCE rather than a cross-check target. A gene with
    neither a protein_accession nor derivable coordinates is skipped, same
    as before. The single-record CLI command and the batch backfill command
    both call this so there is exactly one place this logic lives."""
    record_dir = db_root / phylum / order_or_family / record_id
    record = yaml.safe_load((record_dir / "metadata.yaml").read_text())

    sequences: dict[int, str] = {}
    for gene in record.get("genes", []):
        if not gene.get("present", True):
            continue
        if gene.get("protein_accession"):
            client, bare_accession = _client_for(gene["protein_accession"], ncbi, uniprot)
            sequences[gene["gene_index"]] = client.fetch_protein_sequence(bare_accession)
            continue
        derived = _independent_translation(record, gene, ncbi)
        if derived:
            sequences[gene["gene_index"]] = derived

    gff_export.write_gff3(record, out_path=record_dir / "locus.gff3")
    gff_export.write_genbank(record, sequences, out_path=record_dir / "locus.gbk", ncbi=ncbi)
    gff_export.write_proteins_fasta(record, sequences, out_path=record_dir / "proteins.faa")


def _backfill_records(
    db_root: Path, records: list[tuple[str, str, str]], ncbi: NcbiClient, uniprot: UniprotClient,
) -> tuple[list[tuple[str, str, str]], list[tuple[tuple[str, str, str], str]]]:
    """Run build_gff_for_record for each already-selected record, isolating each
    record's failure so one live-fetch error (a suppressed accession, a transient
    NCBI outage) never aborts the rest of the batch -- the same discipline this
    project's genome-acquisition and batch-runner fixes already established.

    Selection is the caller's job. Both backfill selectors (missing proteins.faa
    and stale locus.gbk) feed this one loop, so the failure-isolation discipline
    exists in exactly one place and cannot drift between them."""
    succeeded: list[tuple[str, str, str]] = []
    failed: list[tuple[tuple[str, str, str], str]] = []
    for phylum, order_or_family, record_id in records:
        identifier = (phylum, order_or_family, record_id)
        try:
            build_gff_for_record(db_root, phylum, order_or_family, record_id, ncbi, uniprot)
        except Exception as exc:  # noqa: BLE001 -- recorded, never swallowed silently
            failed.append((identifier, str(exc)))
            continue
        succeeded.append(identifier)
    return succeeded, failed


def backfill_missing_proteins_faa(
    db_root: Path, ncbi: NcbiClient, uniprot: UniprotClient,
) -> tuple[list[tuple[str, str, str]], list[tuple[tuple[str, str, str], str]]]:
    """Run build_gff_for_record for every record find_records_missing_proteins_faa
    reports. Kept as its own named function (rather than folded into its two
    callers) because it is the default `backfill-gff` behavior and is imported
    by name elsewhere; the loop itself lives in _backfill_records."""
    return _backfill_records(db_root, find_records_missing_proteins_faa(db_root), ncbi, uniprot)


def _cmd_build_gff(args: argparse.Namespace) -> int:
    config = _config(args)
    ncbi, uniprot = _make_clients(config)
    build_gff_for_record(config.db_root, args.phylum, args.order_or_family, args.record_id, ncbi, uniprot)
    record_dir = config.db_root / args.phylum / args.order_or_family / args.record_id
    print(f"wrote {record_dir / 'locus.gff3'}, {record_dir / 'locus.gbk'}, {record_dir / 'proteins.faa'}")
    return 0


def _cmd_backfill_gff(args: argparse.Namespace) -> int:
    config = _config(args)
    ncbi, uniprot = _make_clients(config)
    if args.stale_gbk:
        records = find_records_with_stale_locus_gbk(config.db_root)
    else:
        records = find_records_missing_proteins_faa(config.db_root)
    succeeded, failed = _backfill_records(config.db_root, records, ncbi, uniprot)
    for phylum, order_or_family, record_id in succeeded:
        record_dir = config.db_root / phylum / order_or_family / record_id
        proteins_path = record_dir / "proteins.faa"
        for segment_id in placeholder_segments(record_dir / "locus.gbk"):
            print(
                f"WARNING {phylum}/{order_or_family}/{record_id} segment {segment_id} "
                "has no real nucleotide sequence (all-N placeholder)"
            )
        if ">" not in proteins_path.read_text():
            print(
                f"WARNING wrote 0 sequences for {phylum}/{order_or_family}/{record_id} "
                "(no gene has a protein_accession -- proteins.faa is empty)"
            )
        else:
            print(f"backfilled {phylum}/{order_or_family}/{record_id}")
    for (phylum, order_or_family, record_id), message in failed:
        print(f"FAILED {phylum}/{order_or_family}/{record_id}: {message}")
    print(f"{len(succeeded)} succeeded, {len(failed)} failed")
    return 0


def _cmd_draw_locus(args: argparse.Namespace) -> int:
    # Imported here, not at module scope: `draw` imports `pygenomeviz` at import
    # time, which is slow and (unlike this module's other, core dependencies) not
    # declared as a hard dependency in pyproject.toml -- only `draw-locus`/
    # `draw-synteny` need it, so every other `matpredict` command (including
    # `detect`) must not pay for loading it.
    from MATPredict.db import draw

    config = _config(args)
    record_dir = config.db_root / args.phylum / args.order_or_family / args.record_id
    gbk_path = record_dir / "locus.gbk"
    out_path = Path(args.out) if args.out else record_dir / "locus.png"
    draw.draw_locus(gbk_path, out_path)
    print(f"wrote {out_path}")
    return 0


def _cmd_draw_synteny(args: argparse.Namespace) -> int:
    # See _cmd_draw_locus's comment: kept out of module scope for the same reason
    # (synteny.py itself imports draw.py, which imports pygenomeviz).
    from MATPredict.db import synteny

    config = _config(args)
    out_path = Path(args.out)
    synteny.draw_synteny(args.record_ids, config.db_root, out_path)
    print(f"wrote {out_path}")
    return 0


def _cmd_build_duckdb(args: argparse.Namespace) -> int:
    config = _config(args)
    out_path = Path(args.out) if args.out else config.db_root / "matpredict.duckdb"
    build_duckdb(db_root=config.db_root, out_path=out_path)
    print(f"built {out_path}")
    return 0


def register_subcommands(subparsers: argparse._SubParsersAction) -> None:
    """Register `curate-db` and its actions onto the top-level parser."""
    curate_db = subparsers.add_parser("curate-db", help="Curate the MAT locus reference database")
    action = curate_db.add_subparsers(dest="action", required=True)

    propose = action.add_parser("propose")
    propose.add_argument("--phylum", required=True)
    propose.add_argument("--record-file", required=True)
    propose.set_defaults(func=_cmd_propose)

    validate = action.add_parser("validate")
    validate.add_argument("--phylum", required=True)
    validate.add_argument("--record-id", required=True)
    validate.set_defaults(func=_cmd_validate)

    accept = action.add_parser("accept")
    accept.add_argument("--phylum", required=True)
    accept.add_argument("--order-or-family", required=True)
    accept.add_argument("--record-id", required=True)
    accept.set_defaults(func=_cmd_accept)

    reject = action.add_parser("reject")
    reject.add_argument("--phylum", required=True)
    reject.add_argument("--record-id", required=True)
    reject.add_argument("--reason", required=True)
    reject.set_defaults(func=_cmd_reject)

    build_gff = action.add_parser("build-gff")
    build_gff.add_argument("--phylum", required=True)
    build_gff.add_argument("--order-or-family", required=True)
    build_gff.add_argument("--record-id", required=True)
    build_gff.set_defaults(func=_cmd_build_gff)

    backfill_gff = action.add_parser("backfill-gff")
    backfill_gff.add_argument(
        "--stale-gbk", action="store_true",
        help="Select records whose locus.gbk has zero CDS features (pre-upgrade output) "
             "instead of records with a missing/empty proteins.faa",
    )
    backfill_gff.set_defaults(func=_cmd_backfill_gff)

    draw_locus = action.add_parser("draw-locus")
    draw_locus.add_argument("--phylum", required=True)
    draw_locus.add_argument("--order-or-family", required=True)
    draw_locus.add_argument("--record-id", required=True)
    draw_locus.add_argument("--out", required=False, help="Output image path (default: <record_dir>/locus.png)")
    draw_locus.set_defaults(func=_cmd_draw_locus)

    draw_synteny = action.add_parser("draw-synteny")
    draw_synteny.add_argument("--record-ids", required=True, nargs="+", help="2 or more curated record ids to compare")
    draw_synteny.add_argument("--out", required=True, help="Output clinker plot HTML path")
    draw_synteny.set_defaults(func=_cmd_draw_synteny)

    build_db = action.add_parser("build-duckdb")
    build_db.add_argument("--out", required=False)
    build_db.set_defaults(func=_cmd_build_duckdb)

    release = action.add_parser("release")
    release.set_defaults(func=lambda a: (_ for _ in ()).throw(NotImplementedError("wire release cut CLI in Task 14")))
