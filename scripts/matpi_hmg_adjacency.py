#!/usr/bin/env python3
"""How often does an HMG-box gene lie near a matPi-like (homeodomain) region on the same contig?

Read-only discovery tool. Combines
  - matpi_hits.tsv from matpi_tblastn.py (genome-level homeodomain-like regions, annotated or not), and
  - scan_genes.tsv from pm_lead_followup.py scan (all genes with an HMG_box / HMG_box_2 Pfam hit),
with gene coordinates from the BFD GFF3. For every HMG gene it finds the nearest region on the same contig.
Output (--out): adjacency_pairs.tsv (genome, hmg_gene, contig, region_rank, region_identity, distance_bp) for pairs within
--max-dist, and adjacency_summary.tsv (genome, hmg_genes, hmg_genes_near_region, regions_tested).
The summary gives the background rate: in a genome with several HMG genes and several homeodomain regions, how often a pair
falls within --max-dist by chance on one contig. Assembly fragmentation hides pairs, so a missing pair is not evidence of
absence.
With --perms N the same number of HMG genes is placed at random gene positions N times in every genome; adjacency_null.tsv gives the
observed number of genomes with an HMG gene within --null-dist of a region, the expected number and p (share of draws >= observed).
Usage: matpi_hmg_adjacency.py --hits matpi_hits.tsv --genes SCAN_GENES.tsv [SCAN_GENES.tsv ...] --genomes GENOMES.tsv --out OUT
"""
import argparse
import csv
import gzip
import os
import random
import sys


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_gene_coords(gff, wanted):
    out = {}
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene":
                gid = f[8].split("ID=")[1].split(";")[0]
                if gid in wanted:
                    out[gid] = (f[0], int(f[3]), int(f[4]))
    return out


def read_all_genes(gff):
    out = []
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene":
                out.append((f[0], int(f[3]), int(f[4])))
    return out


def gap(a0, a1, b0, b1):
    return max(0, max(a0, b0) - min(a1, b1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True)
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--genomes", required=True, help="TSV name, fasta, gff3, matMc gene")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-dist", type=int, default=50000)
    ap.add_argument("--perms", type=int, default=2000, help="random placements of the HMG genes for the null (0 = skip)")
    ap.add_argument("--null-dist", type=int, default=20000, help="distance used for the permutation test")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    gff = {l.split("\t")[0]: l.split("\t")[2] for l in open(a.genomes)}
    hits = {}
    for r in csv.DictReader(open(a.hits), delimiter="\t"):
        hits.setdefault(r["genome"], []).append(r)
    hmg = {}
    for p in a.genes:
        for r in csv.DictReader(open(p), delimiter="\t"):
            if "HMG" in r["domain"]:
                hmg.setdefault(r["genome"], set()).add(r["gene"])
    fp = open(os.path.join(a.out, "adjacency_pairs.tsv"), "w")
    fs = open(os.path.join(a.out, "adjacency_summary.tsv"), "w")
    fp.write("genome\thmg_gene\tcontig\tregion_rank\tregion_identity\tdistance_bp\n")
    fs.write("genome\thmg_genes\thmg_genes_near_region\tregions_tested\n")
    for g, regs in hits.items():
        coords = read_gene_coords(gff[g], hmg.get(g, set()))
        near = set()
        for gene, (c, s, e) in coords.items():
            for r in regs:
                if r["contig"] != c:
                    continue
                d = gap(s, e, int(r["start"]), int(r["end"]))
                if d <= a.max_dist:
                    near.add(gene)
                    fp.write(f"{g}\t{gene}\t{c}\t{r['rank']}\t{r['best_identity']}\t{d}\n")
        fs.write(f"{g}\t{len(coords)}\t{len(near)}\t{len(regs)}\n")
    fp.close()
    fs.close()
    if a.perms:
        rng = random.Random(a.seed)
        obs, exp = 0, [0] * a.perms
        used = 0
        for g, regs in hits.items():
            allg = read_all_genes(gff[g])
            k = len(hmg.get(g, ()))
            if not k or not allg:
                continue
            used += 1
            by = {}
            for r in regs:
                by.setdefault(r["contig"], []).append((int(r["start"]), int(r["end"])))
            def near(genes):
                return any(c in by and any(gap(s0, e0, x0, x1) <= a.null_dist for x0, x1 in by[c]) for c, s0, e0 in genes)
            coords = read_gene_coords(gff[g], hmg.get(g, set()))
            obs += near(list(coords.values()))
            for i in range(a.perms):
                exp[i] += near(rng.sample(allg, k))
        p_val = (sum(e >= obs for e in exp) + 1) / (a.perms + 1)
        with open(os.path.join(a.out, "adjacency_null.tsv"), "w") as fn:
            fn.write("max_dist_bp\tgenomes\tobserved_genomes\texpected_mean\texpected_max\tp\n")
            fn.write(f"{a.null_dist}\t{used}\t{obs}\t{sum(exp) / a.perms:.2f}\t{max(exp)}\t{p_val:.4f}\n")
        print(open(os.path.join(a.out, "adjacency_null.tsv")).read())
    print(open(os.path.join(a.out, "adjacency_summary.tsv")).read())
    print(open(os.path.join(a.out, "adjacency_pairs.tsv")).read())


if __name__ == "__main__":
    sys.exit(main())
