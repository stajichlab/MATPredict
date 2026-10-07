"""`matpredict reads-type`: idiomorph typing from raw reads."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from MATPredict import logger
from MATPredict.reads.panel import Panel
from MATPredict.reads.typing import type_reads


def register_subcommands(subparsers) -> None:
    p = subparsers.add_parser("reads-type", help="Call the MAT idiomorph(s) of strains from raw reads")
    p.add_argument("--idiomorph", action="append", required=True, metavar="NAME=FASTA",
                   help="Idiomorph reference, repeat once per idiomorph (at least two), e.g. MAT1-1=AB011379.2.fa")
    p.add_argument("--k", type=int, default=31, help="k-mer length (default 31)")
    p.add_argument("--reads", nargs="+", help="FASTQ file(s) of one strain (.gz/.zst accepted)")
    p.add_argument("--sample", default=None, help="Sample name for --reads (default: first file's stem)")
    p.add_argument("--samples", help="TSV with columns sample, reads (comma-separated files); one row per strain")
    p.add_argument("--max-reads", type=int, default=None, help="Use only the first N reads per strain")
    p.add_argument("--out", required=True, help="Output TSV")
    p.set_defaults(func=_cmd_reads_type)


COLUMNS = ["sample", "call", "reads_used", "shared_depth", "flags"]


def _parse_idiomorphs(items: list[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in items:
        name, sep, path = item.partition("=")
        if not sep or not name or not path:
            raise ValueError(f"--idiomorph needs NAME=FASTA, got {item!r}")
        out[name] = path
    return out


def _cmd_reads_type(args: argparse.Namespace) -> int:
    fastas = _parse_idiomorphs(args.idiomorph)
    panel = Panel.from_fastas(fastas, k=args.k)
    jobs: list[tuple[str, list[str]]] = []
    if args.samples:
        with open(args.samples, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                jobs.append((row["sample"], row["reads"].split(",")))
    if args.reads:
        jobs.append((args.sample or Path(args.reads[0]).name.split(".")[0], list(args.reads)))
    if not jobs:
        raise ValueError("give --reads or --samples")
    names = list(fastas)
    header = COLUMNS[:2] + [f"{n}_breadth" for n in names] + [f"{n}_depth" for n in names] + COLUMNS[2:]
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(header)
        for sample, files in jobs:
            res = type_reads(panel, files, max_reads=args.max_reads)
            logger.info("%s: %s", sample, res.call)
            w.writerow([sample, res.call] + [f"{res.breadth[n]:.4f}" for n in names]
                       + [f"{res.depth[n]:.2f}" for n in names]
                       + [res.reads_used, f"{res.shared_depth:.2f}", ",".join(res.flags) or "-"])
            out.flush()
    return 0
