#!/usr/bin/env python3
"""Flank-first idiomorph probe for a small set of related genomes (for example the three Allomyces genomes).

Read-only discovery tool. Needs stage-1 output of flank_first_scan.py (genes_all.tsv, candidates.tsv per genome) and
the funannotate proteomes.

For every HMG_box / Homeodomain candidate gene X in genome A and every other genome B:
  1. reciprocal best hits (diamond blastp, E <= --evalue) between all genes of A and B;
  2. take the genes within --flank genes of X (not X); keep those with a reciprocal best hit in B;
  3. synteny: >= --min-flank of them lie within --span genes of each other on one contig of B (the largest such group);
  4. class of X for B:
       conserved     X has a reciprocal best hit inside the syntenic region
       other_domain  the region holds an HMG/Homeodomain candidate that is not the reciprocal best hit of X
                     (a divergent or replaced gene: idiomorph-like)
       absent        syntenic region found, no HMG/Homeodomain gene in it (X missing in B: idiomorph-like or gene loss)
       no_synteny    fewer than --min-flank flank genes map as one block (not testable)
The class is a lead, not a call. Gene loss, assembly gaps and annotation errors give the same result.

Usage: flank_idiomorph_probe.py --stage1 DIR --out OUT --genomes G1 G2 [G3 ...]
Output: OUT/baseline.tsv (random non-candidate genes through the same test: conserved / not_conserved / no_synteny), OUT/probe.tsv (genome_A, gene, domain, genome_B, class, flank_mapped, region_contig, region_span), OUT/summary.tsv
"""
import argparse
import collections
import gzip
import os
import random
import subprocess
import sys

ANN = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation"


def table(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    return rows[1:]


def read_genes(stage1, g):
    """-> (genes_in_order: list of (gene, contig, idx), idx_of: gene->(contig, idx))"""
    rows = table(f"{stage1}/{g}.genes_all.tsv")
    genes = [(r[0], r[1], int(r[5])) for r in rows]
    return genes, {r[0]: (r[1], int(r[5])) for r in rows}


def read_candidates(stage1, g):
    return {r[0]: r[3] for r in table(f"{stage1}/{g}.candidates.tsv") if r[3] in ("HMG_box", "HMG_box_2", "Homeodomain")}


def rbh(a, b, tmp, evalue, threads):
    """Reciprocal best hits between proteomes a and b (fasta paths, gz ok) -> dict gene_a -> gene_b."""
    def prep(path, tag):
        out = os.path.join(tmp, tag + ".faa")
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fi, open(out, "w") as fo:
            for line in fi:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    name = name.rsplit("-T", 1)[0] if "-T" in name else name
                    fo.write(f">{name}\n")
                else:
                    fo.write(line)
        return out
    fa, fb = prep(a, "A"), prep(b, "B")
    best = {}
    for q, s, tag in ((fa, fb, "ab"), (fb, fa, "ba")):
        subprocess.run(["diamond", "makedb", "--in", s, "-d", os.path.join(tmp, "db_" + tag), "--quiet"], check=True)
        out = os.path.join(tmp, tag + ".tsv")
        subprocess.run(["diamond", "blastp", "-q", q, "-d", os.path.join(tmp, "db_" + tag), "-o", out, "--quiet",
                        "--max-target-seqs", "1", "-e", str(evalue), "--threads", str(threads),
                        "--outfmt", "6", "qseqid", "sseqid", "bitscore"], check=True)
        best[tag] = {}
        for line in open(out):
            f = line.split("\t")
            best[tag].setdefault(f[0], f[1])
    return {x: y for x, y in best["ab"].items() if best["ba"].get(y) == x}


def classify_gene(x, flank_genes, orth, idx_b, cand_b, min_flank, span):
    """Return (class, n_flank_mapped, contig, region_span). Pure function, used by the tests.
    x: gene id; flank_genes: genes near x in A; orth: A->B reciprocal best hits;
    idx_b: gene->(contig, idx) in B; cand_b: set of HMG/Homeodomain genes in B."""
    pos = [idx_b[orth[f]] for f in flank_genes if f in orth and orth[f] in idx_b]
    by = collections.defaultdict(list)
    for c, i in pos:
        by[c].append(i)
    best = (0, None, None, None)
    for c, v in by.items():
        v.sort()
        lo = 0
        for hi in range(len(v)):
            while v[hi] - v[lo] > span:
                lo += 1
            n = hi - lo + 1
            if n > best[0]:
                best = (n, c, v[lo], v[hi])
    n, c, lo, hi = best
    if n < min_flank:
        return "no_synteny", len(pos), "", ""
    inside = [g for g in cand_b if idx_b[g][0] == c and lo - 2 <= idx_b[g][1] <= hi + 2]
    if x in orth and orth[x] in idx_b and idx_b[orth[x]][0] == c and lo - 2 <= idx_b[orth[x]][1] <= hi + 2:
        return "conserved", len(pos), c, f"{lo}-{hi}"
    if inside:
        return "other_domain", len(pos), c, f"{lo}-{hi}"
    return "absent", len(pos), c, f"{lo}-{hi}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage1", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--genomes", nargs="+", required=True)
    ap.add_argument("--flank", type=int, default=5)
    ap.add_argument("--min-flank", type=int, default=3)
    ap.add_argument("--span", type=int, default=25)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--ann", default=ANN, help="folder holding <genome>/predict_results/<genome>.proteins.fa.gz")
    ap.add_argument("--baseline", type=int, default=500, help="random non-candidate genes per genome pair for the baseline")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    G = {g: read_genes(a.stage1, g) for g in a.genomes}
    C = {g: read_candidates(a.stage1, g) for g in a.genomes}
    rows, base = [], []
    rng = random.Random(a.seed)
    for ga in a.genomes:
        for gb in a.genomes:
            if ga == gb:
                continue
            tmp = os.path.join(a.tmp, f"fip_{os.getpid()}_{a.genomes.index(ga)}_{a.genomes.index(gb)}")
            os.makedirs(tmp, exist_ok=True)
            orth = rbh(f"{a.ann}/{ga}/predict_results/{ga}.proteins.fa.gz", f"{a.ann}/{gb}/predict_results/{gb}.proteins.fa.gz",
                       tmp, a.evalue, a.threads)
            genes_a, idx_a = G[ga]
            by_contig = collections.defaultdict(dict)
            for g, c, i in genes_a:
                by_contig[c][i] = g
            idx_b = G[gb][1]
            for x, dom in C[ga].items():
                c, i = idx_a[x]
                fl = [by_contig[c][j] for j in range(i - a.flank, i + a.flank + 1) if j != i and j in by_contig[c]]
                cls, n, rc, sp = classify_gene(x, fl, orth, idx_b, set(C[gb]), a.min_flank, a.span)
                rows.append((ga, x, dom, gb, cls, n, rc, sp))
            others = [g for g, _c, _i in genes_a if g not in C[ga]]
            for x in rng.sample(others, min(a.baseline, len(others))):
                c, i = idx_a[x]
                fl = [by_contig[c][j] for j in range(i - a.flank, i + a.flank + 1) if j != i and j in by_contig[c]]
                cls, n, rc, sp = classify_gene(x, fl, orth, idx_b, set(), a.min_flank, a.span)
                base.append((ga, gb, "conserved" if cls == "conserved" else "not_conserved" if cls == "absent" else cls))
    with open(os.path.join(a.out, "probe.tsv"), "w") as fo:
        fo.write("genome_A\tgene\tdomain\tgenome_B\tclass\tflank_mapped\tregion_contig\tregion_span\n")
        for r in rows:
            fo.write("\t".join(map(str, r)) + "\n")
    cnt = collections.Counter((r[0], r[3], r[4]) for r in rows)
    with open(os.path.join(a.out, "summary.tsv"), "w") as fo:
        fo.write("genome_A\tgenome_B\tclass\tgenes\n")
        for k, v in sorted(cnt.items()):
            fo.write("\t".join(map(str, k)) + f"\t{v}\n")
    bc = collections.Counter(base)
    with open(os.path.join(a.out, "baseline.tsv"), "w") as fo:
        fo.write("genome_A\tgenome_B\tclass\tgenes\n")
        for k, v in sorted(bc.items()):
            fo.write("\t".join(map(str, k)) + f"\t{v}\n")
    print(open(os.path.join(a.out, "summary.tsv")).read())
    print(open(os.path.join(a.out, "baseline.tsv")).read())


if __name__ == "__main__":
    sys.exit(main())
