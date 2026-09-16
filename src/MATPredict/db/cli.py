"""argparse wiring for the `matpredict curate-db` subcommand group."""
from __future__ import annotations

import argparse


def _cmd_placeholder(args: argparse.Namespace) -> int:
    """Filled in by later tasks (propose/validate/accept/reject/build-gff/build-duckdb/release)."""
    print(f"curate-db {args.action}: not yet implemented")
    return 0


def register_subcommands(subparsers: argparse._SubParsersAction) -> None:
    """Register `curate-db` and its actions onto the top-level parser."""
    curate_db = subparsers.add_parser("curate-db", help="Curate the MAT locus reference database")
    curate_db_sub = curate_db.add_subparsers(dest="action", required=True)

    for action in ("propose", "validate", "accept", "reject", "build-gff", "build-duckdb", "release"):
        p = curate_db_sub.add_parser(action)
        p.set_defaults(func=_cmd_placeholder, action=action)
