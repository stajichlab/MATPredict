#!/usr/bin/env python3
"""Apply the positional rule to uncurated Agaricomycete genomes.

Rule (from the ground-truth test): an STE3-like locus is flagged when a
pheromone-precursor candidate lies within +-20 kb:
  T2  = stop-anchored ORF, in-frame Met 20-130 codons upstream of the stop,
        ending C[VITE][IVT][AVMG]; or
  R   = >=2 short ORFs with an identical 8-aa C-terminal tail ending C-x-x-x.
Also reported: strict T (C[VI][IV][AVMG]) at +-10 kb.

Maps each locus to its BFD funannotate protein (GFF3 overlap) to place it in
the step-1 receptor tree, and reports tree clustering of flagged copies.
"""
import csv
import glob
import os
import sys
from collections import defaultdict, Counter

from Bio import Phylo

D = "out_uncur"
S1 = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_receptor_explore"
samp = {r[0]: r for r in csv.reader(open("uncurated_sample.tsv"), delimiter="\t")}


def flagged(r, strict=False):
    if strict:
        return int(r["T_10kb"]) > 0
    return int(r["T2_20kb"]) > 0 or int(r["R_20kb"]) >= 2


def gff_mrnas(path):
    out = defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        pid = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x).get("ID", "")
        out[f[0]].append((int(f[3]), int(f[4]), pid))
    return out


rows = []
per_genome = []
rand_rate = []
for f in sorted(glob.glob(os.path.join(D, "*.loci.tsv"))):
    asm = os.path.basename(f).replace(".loci.tsv", "")
    order, sp, st, prot = samp[asm][1], samp[asm][2], samp[asm][3], samp[asm][4]
    gff = prot.replace(".proteins.fa", ".gff3")
    mr = gff_mrnas(gff) if os.path.exists(gff) else {}
    loci = list(csv.DictReader(open(f), delimiter="\t"))
    rnd = list(csv.DictReader(open(f.replace(".loci.tsv", ".random.tsv")), delimiter="\t"))
    rr = sum(flagged(x) for x in rnd) / max(1, len(rnd))
    rrs = sum(flagged(x, True) for x in rnd) / max(1, len(rnd))
    rand_rate.append((rr, rrs))
    nf = nfs = 0
    for r in loci:
        s, e = int(r["start"]), int(r["end"])
        pids = [p for a, b, p in mr.get(r["contig"], []) if a <= e and b >= s]
        r.update(asm=asm, order=order, species=sp, strain=st, protein=";".join(pids),
                 flag=int(flagged(r)), flag_strict=int(flagged(r, True)))
        nf += r["flag"]
        nfs += r["flag_strict"]
        rows.append(r)
    # clusters: flagged loci within 50 kb of each other on one contig = one P/R locus
    fl = sorted([(r["contig"], int(r["start"])) for r in loci if flagged(r)])
    clusters = 0
    last = (None, -10**9)
    for c, s in fl:
        if c != last[0] or s - last[1] > 50000:
            clusters += 1
        last = (c, s)
    per_genome.append((asm, order, sp, st, len(loci), nf, nfs, clusters, rr, len(loci) * rr))

with open("uncurated_loci.tsv", "w") as fo:
    keys = ["asm", "order", "species", "strain", "contig", "strand", "start", "end", "protein", "flag", "flag_strict",
            "best_ident", "best_ref", "best_ref_ident", "T_10kb", "T2_20kb", "R_20kb", "Hx_20kb"]
    fo.write("\t".join(keys) + "\n")
    for r in rows:
        fo.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
with open("uncurated_per_genome.tsv", "w") as fo:
    fo.write("asm\torder\tspecies\tstrain\tste3_loci\tflagged\tflagged_strict\tflagged_clusters\trandom_window_rate\texpected_by_chance\n")
    for p in per_genome:
        fo.write("\t".join(str(x if not isinstance(x, float) else round(x, 3)) for x in p) + "\n")

print(f"genomes {len(per_genome)}; STE3 loci {len(rows)}")
by = defaultdict(list)
for p in per_genome:
    by[p[1]].append(p)
print("order\tgenomes\tloci/genome(med)\tflagged/genome(med)\tgenomes>=1 flag\tclusters/genome(med)\texpected flags by chance/genome(med)")
import statistics as stt
for o, ps in sorted(by.items()):
    print(o, len(ps), stt.median(p[4] for p in ps), stt.median(p[5] for p in ps), sum(p[5] > 0 for p in ps),
          stt.median(p[7] for p in ps), round(stt.median(p[9] for p in ps), 2), sep="\t")
tot_f = sum(p[5] for p in per_genome)
tot_e = sum(p[9] for p in per_genome)
print(f"total flagged {tot_f}; expected by chance {tot_e:.1f}; strict flagged {sum(p[6] for p in per_genome)}")
print(f"median random-window rate: T2|R@20kb {stt.median(x[0] for x in rand_rate):.3f}; T@10kb {stt.median(x[1] for x in rand_rate):.3f}")

# tree: are flagged copies closer to curated mating receptors than unflagged ones?
tree = Phylo.read(os.path.join(S1, "tree_ft.treefile"), "newick")
tips = {t.name: t for t in tree.get_terminals()}
agar_refs = [t for n, t in tips.items() if n.startswith("REF|5346_") or n.startswith("REF|5334_")]
res = {"flag": [], "noflag": []}
for r in rows:
    for pid in r["protein"].split(";"):
        name = f"{r['asm']}|{pid}"
        if name in tips:
            d = min(tree.distance(tips[name], a) for a in agar_refs)
            res["flag" if r["flag"] else "noflag"].append(d)
for k, v in res.items():
    if v:
        print(f"tree: {k} copies placed n={len(v)}; median distance to nearest Coprinopsis/Schizophyllum mating receptor {stt.median(v):.2f}")
