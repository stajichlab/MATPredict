"""F4: false-positive rate of the strict-CAAX rule (T within 10 kb).

Sets, each counted as loci (or windows) flagged / total:
  gt_mating        curated mating receptors in the 6 ground-truth genomes
  gt_nonmating     other STE3-like copies in those genomes
  agaricales_rcpt  STE3-like loci in 15 Agaricales CAAX-panel genomes (mix)
  agaricales_rand  random windows (same size, away from STE3) in those genomes
  basidio_rand     random windows in the ground-truth + 52 uncurated genomes
  asco_rcpt        STE3-like loci in 30 Pezizomycotina genomes
  asco_rand        random windows in those genomes
  rust_rcpt        STE3-like loci in rust genomes (earlier 4 + new)
  rust_rand        random windows in rust genomes
Column used: T_10kb (count of strict-CAAX ORF candidates within +-10 kb).
"""
import csv
import glob
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(HERE, "..", "2026-09-27_pheromone_positional")


def asms(path):
    return {l.strip() for l in open(path) if l.strip()}


def rows(pattern):
    out = []
    for f in glob.glob(pattern):
        out += list(csv.DictReader(open(f), delimiter="\t"))
    return out


def flagged(rs, col="T_10kb"):
    return sum(int(r[col]) > 0 for r in rs), len(rs)


sets = {}
gt = list(csv.DictReader(open(os.path.join(PP, "gt_loci_labelled.tsv")), delimiter="\t"))
sets["gt_mating"] = flagged([r for r in gt if r["mating"] in ("1", "True", "mating", "yes")])
sets["gt_nonmating"] = flagged([r for r in gt if r["mating"] not in ("1", "True", "mating", "yes")])
sets["basidio_rand(gt+uncurated52)"] = flagged(rows(os.path.join(PP, "out_gt", "*.random.tsv"))
                                               + rows(os.path.join(PP, "out_uncur", "*.random.tsv")))
sets["uncurated52_rcpt"] = flagged(rows(os.path.join(PP, "out_uncur", "*.loci.tsv")))

neg = os.path.join(HERE, "out_neg")
groups = {"asco": asms(os.path.join(HERE, "asco.txt")),
          "rust_new": asms(os.path.join(HERE, "rust.txt")),
          "agaricales": asms(os.path.join(HERE, "agaricales.txt"))}
for g, names in groups.items():
    L = [r for a in names for r in rows(os.path.join(neg, f"{a}.loci.tsv"))]
    Rw = [r for a in names for r in rows(os.path.join(neg, f"{a}.random.tsv"))]
    done = sum(os.path.exists(os.path.join(neg, f"{a}.loci.tsv")) for a in names)
    sets[f"{g}_rcpt (genomes {done}/{len(names)})"] = flagged(L)
    sets[f"{g}_rand (genomes {done}/{len(names)})"] = flagged(Rw)
old_rust_L = rows(os.path.join(PP, "out_rust", "*.loci.tsv"))
old_rust_R = rows(os.path.join(PP, "out_rust", "*.random.tsv"))
sets["rust_old4_rcpt"] = flagged(old_rust_L)
sets["rust_old4_rand"] = flagged(old_rust_R)

lines = ["set\tflagged\ttotal\trate"]
for k, (a, n) in sets.items():
    lines.append(f"{k}\t{a}\t{n}\t{(a / n if n else float('nan')):.3f}")
txt = "\n".join(lines)
open(os.path.join(HERE, "f4_summary.tsv"), "w").write(txt + "\n")
print(txt)
