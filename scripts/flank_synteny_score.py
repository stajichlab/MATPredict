#!/usr/bin/env python3
"""Flank-first scan, stage 2: is a candidate locus anchored by collinear flanking genes?

Read-only discovery tool. Needs stage-1 output (flank_first_scan.py).

For each group of genomes (genus by default) one reference genome is chosen (highest
N50 from the scaffolds .fai). Every other genome is compared with it by reciprocal best
diamond hits (RBH). For each stage-1 locus:
  - flanks   = the genes in flanks.tsv with in_locus = 0 (up to 2 x window)
  - mapped   = flanks with an RBH ortholog in the reference
  - block    = size of the largest set of mapped flanks whose reference orthologs lie on
               one reference contig with gaps of <= --gap reference genes
  - cand_rbh = how many locus genes have an RBH ortholog in the reference
  - null     = the same block size for the flanks of EVERY gene in the genome (window
               2 x W, locus gene removed). The null shows how collinear an ordinary gene
               neighbourhood is. A locus is "anchored" when its block is >= --min-block.
  - percentile = share of null windows with a block size smaller than the locus block.
Strains of one species share gene order, so a high block alone is expected. The null and
the cand_rbh column are what separate a locus whose own genes diverge (no RBH) from
conserved ones.

Output (in --out): <group>.loci_score.tsv, <group>.null_hist.tsv, <group>.refmap.tsv.
"""
import argparse
import collections
import gzip
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flank_first_scan import parse_gff, read_fasta, run  # noqa: E402

ANN = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation"


def n50(fai):
    lens = sorted((int(l.split("\t")[1]) for l in open(fai)), reverse=True)
    tot, acc = sum(lens), 0
    for x in lens:
        acc += x
        if acc * 2 >= tot:
            return x
    return 0


def proteome(g, tmp):
    """Write one longest protein per gene; return {gene: (contig, index)} and the FASTA path."""
    p = f"{ANN}/{g}/predict_results"
    order, m2g = parse_gff(f"{p}/{g}.gff3.gz")
    pos = {x[0]: (x[1], i) for i, x in enumerate(order)}
    best = {}
    for pid, seq in read_fasta(f"{p}/{g}.proteins.fa.gz"):
        gene = m2g.get(pid, re.sub(r"-T\d+$", "", pid))
        if gene in pos and (gene not in best or len(seq) > len(best[gene])):
            best[gene] = seq.rstrip("*")
    fa = os.path.join(tmp, f"{g}.faa")
    with open(fa, "w") as fo:
        for gene, s in best.items():
            fo.write(f">{gene}\n{s}\n")
    return pos, order, fa


def rbh(qfa, rfa, tmp, threads, evalue):
    def hits(q, d, tag):
        out = os.path.join(tmp, tag + ".tsv")
        run(["diamond", "blastp", "-q", q, "-d", d, "--evalue", str(evalue), "--max-target-seqs", "1",
             "--threads", str(threads), "--query-cover", "50", "--subject-cover", "50",
             "--outfmt", "6", "qseqid", "sseqid", "bitscore", "-o", out, "--tmpdir", tmp, "--quiet"])
        return {l.split("\t")[0]: l.split("\t")[1] for l in open(out)}
    fw = hits(qfa, rfa + ".dmnd", "fw")
    run(["diamond", "makedb", "--in", qfa, "-d", qfa + ".dmnd", "--quiet"])
    bw = hits(rfa, qfa + ".dmnd", "bw")
    return {q: r for q, r in fw.items() if bw.get(r) == q}


def block(refpos, gap):
    """Largest set of (contig, index) points on one contig with consecutive gaps <= gap."""
    if not refpos:
        return 0, None
    pts = sorted(refpos)
    best, cur, bestspan = 0, [pts[0]], None
    for p in pts[1:] + [None]:
        if p is not None and p[0] == cur[-1][0] and p[1] - cur[-1][1] <= gap:
            cur.append(p)
            continue
        if len(cur) > best:
            best, bestspan = len(cur), (cur[0][0], cur[0][1], cur[-1][1])
        cur = [p] if p is not None else []
    return best, bestspan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage1", required=True, help="stage-1 per_genome directory")
    ap.add_argument("--list", required=True, help="genome directory names, one per line")
    ap.add_argument("--out", required=True)
    ap.add_argument("--group-by", choices=["genus", "all"], default="genus")
    ap.add_argument("--min-group", type=int, default=3)
    ap.add_argument("--ref", help="reference genome (only with --group-by all)")
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--gap", type=int, default=5, help="max reference-gene gap inside a block")
    ap.add_argument("--min-block", type=int, default=6)
    ap.add_argument("--evalue", type=float, default=1e-10)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    genomes = [l.strip() for l in open(a.list) if l.strip()]
    groups = collections.defaultdict(list)
    for g in genomes:
        groups["all" if a.group_by == "all" else g.split("_")[0]].append(g)

    for grp, members in sorted(groups.items()):
        if len(members) < a.min_group:
            print(f"skip group {grp}: {len(members)} genome(s)")
            continue
        if a.ref and a.group_by == "all":
            ref = a.ref
        else:
            ref = max(members, key=lambda g: n50(f"{ANN}/{g}/predict_results/{g}.scaffolds.fa.gz.fai"))
        tmp = os.path.join(a.tmp, f"ffs2_{grp}")
        os.makedirs(tmp, exist_ok=True)
        rpos, rorder, rfa = proteome(ref, tmp)
        run(["diamond", "makedb", "--in", rfa, "-d", rfa + ".dmnd", "--quiet"])
        rinv = {v: k for k, v in rpos.items()}  # (contig,index) -> gene
        W = a.window
        hist = collections.defaultdict(collections.Counter)  # genome -> block size -> n windows
        with open(f"{a.out}/{grp}.loci_score.tsv", "w") as fo, open(f"{a.out}/{grp}.refmap.tsv", "w") as fm:
            fo.write("group\tref\tgenome\tlocus\tflag\tclasses\tpfam\tn_flank\tn_mapped\tblock\tblock_ref_contig\t"
                     "block_ref_start\tblock_ref_end\tn_locus_genes\tlocus_genes_with_rbh\tnull_percentile\n")
            fm.write("genome\tgene\tref_gene\n")
            for g in members:
                if not os.path.exists(f"{a.stage1}/{g}.loci.tsv"):
                    print(f"no stage-1 output for {g}")
                    continue
                pos, order, fa = proteome(g, tmp)
                m = rbh(fa, rfa, tmp, a.threads, a.evalue) if g != ref else {x: x for x in pos}
                for q, r in m.items():
                    fm.write(f"{g}\t{q}\t{r}\n")
                # null: block size for the flanks of every gene
                null = collections.Counter()
                for i, gene in enumerate(o[0] for o in order):
                    pts = []
                    for j in range(max(0, i - W), min(len(order), i + W + 1)):
                        if j == i or order[j][1] != order[i][1]:
                            continue
                        r = m.get(order[j][0])
                        if r in rpos:
                            pts.append(rpos[r])
                    null[block(pts, a.gap)[0]] += 1
                hist[g] = null
                tot = sum(null.values())
                # loci
                loci = {}
                for l in list(open(f"{a.stage1}/{g}.loci.tsv"))[1:]:
                    f = l.rstrip("\n").split("\t")
                    loci[f[0]] = f
                members_of = collections.defaultdict(list)
                flanks = collections.defaultdict(list)
                for l in list(open(f"{a.stage1}/{g}.flanks.tsv"))[1:]:
                    lid, gene, off, inl = l.rstrip("\n").split("\t")
                    (members_of if inl == "1" else flanks)[lid].append(gene)
                for lid, f in loci.items():
                    pts = [rpos[m[x]] for x in flanks[lid] if x in m and m[x] in rpos]
                    b, span = block(pts, a.gap)
                    below = sum(v for k, v in null.items() if k < b)
                    rb = sum(1 for x in members_of[lid] if x in m)
                    fo.write("\t".join(map(str, [grp, ref, g, lid, f[7], f[5], f[6], len(flanks[lid]), len(pts), b,
                                                 span[0] if span else "", span[1] if span else "",
                                                 span[2] if span else "", len(members_of[lid]), rb,
                                                 f"{below / tot:.4f}"])) + "\n")
                print(f"{grp}: {g} done ({len(loci)} loci)", flush=True)
        with open(f"{a.out}/{grp}.null_hist.tsv", "w") as fh:
            fh.write("genome\tblock\tn_windows\n")
            for g, c in hist.items():
                for k, v in sorted(c.items()):
                    fh.write(f"{g}\t{k}\t{v}\n")
        print(f"group {grp}: ref={ref}, members={len(members)}")


if __name__ == "__main__":
    main()
