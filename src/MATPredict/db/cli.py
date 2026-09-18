"""argparse wiring for the `matpredict curate-db` subcommand group."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests
import yaml

from MATPredict.config import MatpredictConfig
from MATPredict.db import gff_export
from MATPredict.db.build_duckdb import build as build_duckdb
from MATPredict.db.curate import accept_candidate, propose_candidate, reject_candidate
from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.schema import validate_gene_vocabulary
from MATPredict.db.uniprot_client import UniprotClient
from MATPredict.db.validate import _client_for, validate_record


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


def _cmd_build_gff(args: argparse.Namespace) -> int:
    config = _config(args)
    record_dir = config.db_root / args.phylum / args.order_or_family / args.record_id
    record_path = record_dir / "metadata.yaml"
    record = yaml.safe_load(record_path.read_text())

    ncbi, uniprot = _make_clients(config)

    sequences: dict[int, str] = {}
    for gene in record.get("genes", []):
        if not gene.get("present", True) or not gene.get("protein_accession"):
            continue
        client, bare_accession = _client_for(gene["protein_accession"], ncbi, uniprot)
        sequences[gene["gene_index"]] = client.fetch_protein_sequence(bare_accession)

    gff3_path = record_dir / "locus.gff3"
    gbk_path = record_dir / "locus.gbk"
    proteins_path = record_dir / "proteins.faa"

    gff_export.write_gff3(record, out_path=gff3_path)
    gff_export.write_genbank(record, sequences, out_path=gbk_path)
    gff_export.write_proteins_fasta(record, sequences, out_path=proteins_path)

    print(f"wrote {gff3_path}, {gbk_path}, {proteins_path}")
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

    build_db = action.add_parser("build-duckdb")
    build_db.add_argument("--out", required=False)
    build_db.set_defaults(func=_cmd_build_duckdb)

    release = action.add_parser("release")
    release.set_defaults(func=lambda a: (_ for _ in ()).throw(NotImplementedError("wire release cut CLI in Task 14")))
