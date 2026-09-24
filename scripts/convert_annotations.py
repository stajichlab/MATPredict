"""Convert annotation GenBank files into the (genome, proteome) pairs `detect` needs.

Thin driver: all the logic lives in `MATPredict.detect.annotation_export`, which
is where the defline-format contract is documented and tested.

    pixi run python scripts/convert_annotations.py \
        --out-dir $SCRATCH/zygolife/fasta \
        /bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/annotate/*/annotate_results/*.gbk

Writes to $SCRATCH, never to shared storage: 899 ZygoLife genomes is tens of GB
of FASTA. Already-converted genomes are skipped, so a killed run resumes.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from MATPredict.detect.annotation_export import convert_genbank


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gbk", nargs="+", help="annotation .gbk files")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--force", action="store_true",
        help="re-convert genomes that already have both output files",
    )
    args = parser.parse_args(argv)

    out_dir = Path(args.out_dir)
    done = skipped = failed = 0
    for path in args.gbk:
        gbk = Path(path)
        pair = (out_dir / f"{gbk.stem}.fna", out_dir / f"{gbk.stem}.faa")
        if not args.force and all(p.exists() and p.stat().st_size for p in pair):
            skipped += 1
            continue
        try:
            result = convert_genbank(gbk, out_dir)
        except Exception as exc:  # one bad annotation must not stop the batch
            print(f"FAIL {gbk.stem}: {type(exc).__name__}: {exc}", file=sys.stderr)
            failed += 1
            continue
        done += 1
        print(f"{gbk.stem}\t{result.contigs} contigs\t{result.proteins} proteins", flush=True)
    print(f"converted {done}, skipped {skipped}, failed {failed}", file=sys.stderr)
    return 1 if failed and not done else 0


if __name__ == "__main__":
    raise SystemExit(main())
