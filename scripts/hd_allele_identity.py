#!/usr/bin/env python3
"""Pairwise protein identity of HD1-like and HD2 genes: within species, within genus, between genera.

Inputs: stage-1 output (control class HD1 / HD2 per gene) for the Glomeromycota genomes, and the
curated Basidiomycota HD1 / HD2 proteins (db/Basidiomycota, gene names HD1, HD2) for comparison.
Method: one mafft alignment per gene class and group; identity = identical columns / columns where
both sequences have a residue (aligned identity), and also identical / shorter length.
Output (in --out): hd_pairs_<set>_<class>.tsv, hd_allele_summary.tsv.
Species = first two name tokens of the genome directory (Genus_species); "sp." names are their own species.
"""
import argparse
import collections
import itertools
import os
import re
import statistics
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homeodomain_tree import Prot, table  # noqa: E402


def fasta_full(path):
    h, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if h is not None:
                yield h, "".join(buf)
            h, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if h is not None:
        yield h, "".join(buf)


def mafft(seqs, tmp, tag, threads):
    fa = os.path.join(tmp, tag + ".faa")
    with open(fa, "w") as fo:
        for k, s in seqs.items():
            fo.write(f">{k}\n{s}\n")
    r = subprocess.run(["mafft", "--auto", "--thread", str(threads), "--quiet", fa], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stderr[-500:])
    aln, k = {}, None
    for line in r.stdout.splitlines():
        if line.startswith(">"):
            k = line[1:].strip(); aln[k] = []
        else:
            aln[k].append(line.strip())
    return {k: "".join(v).upper() for k, v in aln.items()}


def ident(a, b):
    both = same = 0
    for x, y in zip(a, b):
        if x != "-" and y != "-":
            both += 1
            same += x == y
    la, lb = sum(c != "-" for c in a), sum(c != "-" for c in b)
    return (same / both if both else 0.0), (same / min(la, lb) if min(la, lb) else 0.0), both


def species_of(name):
    t = name.split("_")
    return t[0] + "_" + t[1] if len(t) > 1 and t[1] != "sp." else name


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--db", default="/bigdata/stajichlab/jstajich/projects/MATPredict/db/Basidiomycota")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    D = f"{a.base}/stage1_all/per_genome"
    P = Prot()

    sets = {"AMF": {"HD1": {}, "HD2": {}}, "DIK": {"HD1": {}, "HD2": {}}}
    meta = {}  # label -> (genome, species, genus)
    for g in [x.strip() for x in open(f"{a.base}/list_glomeromycota.txt") if x.strip()]:
        for c in table(f"{D}/{g}.candidates.tsv"):
            if c[5] in ("HD1", "HD2") and float(c[7]) <= 1e-20:
                lab = f"{g}|{c[0]}"
                s = P.get(g).get(c[0])
                if s:
                    sets["AMF"][c[5]][lab] = s
                    meta[lab] = (g, species_of(g), g.split("_")[0])
    for root, _, files in os.walk(a.db):
        if "proteins.faa" in files and "HD" in os.path.basename(root):
            rec = os.path.basename(root)
            sp = None
            for line in open(os.path.join(root, "metadata.yaml")):
                m = re.match(r"\s+species:\s*(.+)", line)
                if m:
                    sp = m.group(1).strip().replace(" ", "_"); break
            for h, s in fasta_full(os.path.join(root, "proteins.faa")):
                m = re.search(r"name=(HD[12])\b", h)
                if m:
                    lab = f"{rec}|{m.group(1)}"
                    sets["DIK"][m.group(1)][lab] = s
                    meta[lab] = (rec, sp, sp.split("_")[0])

    rows = []
    for setname, classes in sets.items():
        for cl, seqs in classes.items():
            if len(seqs) < 2:
                continue
            aln = mafft(seqs, a.tmp, f"{setname}_{cl}", a.threads)
            path = os.path.join(a.out, f"hd_pairs_{setname}_{cl}.tsv")
            with open(path, "w") as fo:
                fo.write("a\tb\tspecies_a\tspecies_b\trelation\taligned_identity\tidentity_over_shorter\taligned_cols\n")
                for x, y in itertools.combinations(aln, 2):
                    ma, mb = meta[x], meta[y]
                    rel = "same_species" if ma[1] == mb[1] else "same_genus" if ma[2] == mb[2] else "diff_genus"
                    i1, i2, n = ident(aln[x], aln[y])
                    fo.write(f"{x}\t{y}\t{ma[1]}\t{mb[1]}\t{rel}\t{i1:.4f}\t{i2:.4f}\t{n}\n")
                    rows.append((setname, cl, rel, i1))
    with open(os.path.join(a.out, "hd_allele_summary.tsv"), "w") as fo:
        fo.write("set\tclass\trelation\tpairs\tmedian_identity\tmin\tmax\tq25\tq75\n")
        by = collections.defaultdict(list)
        for s, c, r, i in rows:
            by[(s, c, r)].append(i)
        for k, v in sorted(by.items()):
            v = sorted(v)
            q = statistics.quantiles(v, n=4) if len(v) >= 2 else [v[0], v[0], v[0]]
            fo.write("\t".join(map(str, [*k, len(v), f"{statistics.median(v):.3f}", f"{v[0]:.3f}", f"{v[-1]:.3f}",
                                          f"{q[0]:.3f}", f"{q[2]:.3f}"])) + "\n")
    print(open(os.path.join(a.out, "hd_allele_summary.tsv")).read())
    print({s: {c: len(v) for c, v in cl.items()} for s, cl in sets.items()})


if __name__ == "__main__":
    main()
