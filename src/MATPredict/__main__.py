"""MATPredict CLI entry point."""
from __future__ import annotations

import argparse
import sys

from MATPredict import __author__, __version__, logger
from MATPredict.db.cli import register_subcommands


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level `matpredict` argument parser."""
    parser = argparse.ArgumentParser(prog="matpredict", description="MATPredict: fungal MAT locus tooling")
    parser.add_argument("-V", "--version", action="version", version=__version__)
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose/debug logging")
    subparsers = parser.add_subparsers(dest="command", required=True)
    register_subcommands(subparsers)
    return parser


def main(args: list[str] | None = None) -> int:
    """Tool for building and querying the MAT locus reference database."""
    parser = build_parser()
    argv = args if args is not None else sys.argv[1:]
    if not argv or "help" in argv or "-h" in argv or "--help" in argv:
        parser.print_help()
        return 0
    parsed = parser.parse_args(argv)
    if parsed.verbose:
        logger.setLevel("DEBUG")
    try:
        return parsed.func(parsed)
    except Exception as err:  # noqa: BLE001 - top-level CLI error boundary
        logger.error(err)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
