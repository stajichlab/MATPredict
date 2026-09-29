#!/usr/bin/env python3
"""Write a classifier's HMG-paralog negative set (`paralog_negatives.faa`).

The MAT-gene gate threshold is the 95th percentile of these proteins' best
classifier scores against each build's HMMs (curator ruling 2026-09-29; see
MATPredict.detect.classifier_build.gate_threshold). The set is the one the
2026-09-28 fragment validation used as negatives
(results/2026-09-28_validation_f3_f4/f3_fragment_loo.py): the non-locus HMG
copies that fell OUTSIDE the sexM and sexP clades ("other_HMG") in the
2026-09-27 FastTree HMG-box tree, full-length proteins. They are never used to
train the HMMs.

Usage: make_paralog_negatives.py TIP_NAMES.tsv ALL_PROTEINS.faa OUT.faa
"""
import csv
import sys
from pathlib import Path


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main(tips_path, proteins_path, out_path):
    tips = [r for r in csv.DictReader(open(tips_path), delimiter="\t")
            if r["status"] == "nonlocus" and r["tree_clade"] == "other_HMG"]
    proteins = read_fasta(proteins_path)
    keep = {r["source"]: proteins[r["source"]] for r in tips if r["source"] in proteins}
    with open(out_path, "w") as fh:
        for name in sorted(keep):
            seq = keep[name].replace("*", "")
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    print(f"{len(tips)} other_HMG non-locus tips; {len(keep)} with a protein -> {out_path}",
          file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:4])
