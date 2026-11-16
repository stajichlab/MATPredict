#!/usr/bin/env python3
"""Homeodomain tree: do AMF HD1/HD2 and adjacent homeodomain pairs in other lineages group together?

Read-only discovery tool (needs stage-1 output of flank_first_scan.py).

Steps
 1. Collect proteins (labels carry the group):
      DIK_HD1 / DIK_HD2   curated Basidiomycota HD1/HD2 (db/Basidiomycota, gene names HD1, HD2)
      AMF_HD1 / AMF_HD2   Glomeromycota genes with a Ropars HD1-like / HD2 control hit
      AMFbg               all other Pfam Homeodomain genes of one Rhizophagus reference genome
      ZOO_all             every Pfam Homeodomain gene of Zoopagomycota
      <LIN>_pair          adjacent homeodomain pair genes (two Homeodomain genes with no more than
                          one gene between them, no third within 6 genes) from other lineages
      MORbg               Pfam Homeodomain genes of a few Mortierellomycota genomes (background)
 2. hmmsearch (Pfam Homeodomain, GA cutoff), keep the best domain per protein.
 3. hmmalign the domains to the Pfam model, drop columns with > 50% gaps.
 4. Output FASTA/alignment files for tree building (done by the SLURM wrapper).
"""
import argparse
import collections
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flank_first_scan import parse_gff, read_fasta, run  # noqa: E402

ANN = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation"
LINS = {"zoopagomycota": "ZOO", "basidiobolomycota": "BAS", "chytridiomycota": "CHY",
        "entomophthoromycota": "ENT", "kickxellomycota": "KIC"}


class Prot:
    def __init__(self):
        self.cache = {}

    def get(self, g):
        if g not in self.cache:
            p = f"{ANN}/{g}/predict_results"
            _, m2g = parse_gff(f"{p}/{g}.gff3.gz")
            best = {}
            for pid, seq in read_fasta(f"{p}/{g}.proteins.fa.gz"):
                gene = m2g.get(pid, re.sub(r"-T\d+$", "", pid))
                if gene not in best or len(seq) > len(best[gene]):
                    best[gene] = seq.rstrip("*")
            self.cache[g] = best
        return self.cache[g]


def table(path, skip=1):
    return [l.rstrip("\n").split("\t") for l in list(open(path))[skip:]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="results/2026-11-13_flank_first_scan")
    ap.add_argument("--db", default="/bigdata/stajichlab/jstajich/projects/MATPredict/db/Basidiomycota")
    ap.add_argument("--hmm", required=True, help="Pfam Homeodomain HMM (single model)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref-amf", default="Rhizophagus_irregularis_DAOM-197198")
    ap.add_argument("--mor-n", type=int, default=6)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    D = f"{a.base}/stage1_all/per_genome"
    P = Prot()
    seqs = {}

    def add(label, s):
        if s and len(s) >= 40:
            seqs[label] = s

    # Dikarya curated
    n_dik = 0
    for root, _, files in os.walk(a.db):
        if "proteins.faa" in files and "HD" in os.path.basename(root):
            rec = os.path.basename(root)
            for h, s in read_fasta(os.path.join(root, "proteins.faa")):
                m = re.search(r"name=(HD[12])", h + " ")
                if m is None:
                    # read_fasta keeps only the first token; recover full header
                    pass
            for line_h, s in _fasta_full(os.path.join(root, "proteins.faa")):
                m = re.search(r"name=(HD[12])\b", line_h)
                if m:
                    add(f"DIK_{m.group(1)}|{rec}", s); n_dik += 1

    # AMF HD hits, Glomeromycota
    glo = [x.strip() for x in open(f"{a.base}/list_glomeromycota.txt") if x.strip()]
    for g in glo:
        for c in table(f"{D}/{g}.candidates.tsv"):
            if c[5] in ("HD1", "HD2"):
                add(f"AMF_{c[5]}|{g}|{c[0]}", P.get(g).get(c[0]))
    # AMF background homeodomain genes in the reference genome
    hd_ctl = {c[0] for c in table(f"{D}/{a.ref_amf}.candidates.tsv") if c[5] in ("HD1", "HD2")}
    for c in table(f"{D}/{a.ref_amf}.candidates.tsv"):
        if c[3] == "Homeodomain" and c[0] not in hd_ctl:
            add(f"AMFbg|{a.ref_amf}|{c[0]}", P.get(a.ref_amf).get(c[0]))

    # other lineages
    for lin, tag in LINS.items():
        gs = [x.strip() for x in open(f"{a.base}/list_{lin}.txt") if x.strip()]
        for g in gs:
            gene_tab = table(f"{D}/{g}.genes_all.tsv")
            st = {r[0]: (r[1], r[4], int(r[5])) for r in gene_tab}
            hom = sorted([c[0] for c in table(f"{D}/{g}.candidates.tsv") if c[3] == "Homeodomain"],
                         key=lambda x: st[x][2])
            if tag == "ZOO":
                for h in hom:
                    add(f"ZOO_all|{g}|{h}", P.get(g).get(h))
            for x, y in zip(hom, hom[1:]):
                if st[x][0] != st[y][0] or st[y][2] - st[x][2] > 2:
                    continue
                near = [h for h in hom if st[h][0] == st[x][0] and abs(st[h][2] - st[x][2]) <= 6]
                if len(near) != 2:
                    continue
                s1, s2 = st[x][1], st[y][1]
                rel = "div" if (s1, s2) == ("-", "+") else "conv" if (s1, s2) == ("+", "-") else "same"
                for h in (x, y):
                    add(f"{tag}_pair{rel}|{g}|{h}", P.get(g).get(h))
    # Mortierellomycota background
    mor = [x.strip() for x in open(f"{a.base}/list_mortierellomycota.txt") if x.strip()]
    step = max(1, len(mor) // a.mor_n)
    for g in mor[::step][: a.mor_n]:
        for c in table(f"{D}/{g}.candidates.tsv"):
            if c[3] == "Homeodomain":
                add(f"MORbg|{g}|{c[0]}", P.get(g).get(c[0]))

    # write FASTA, run hmmsearch
    allfa = os.path.join(a.out, "homeodomain_proteins.faa")
    with open(allfa, "w") as fo:
        for k, s in seqs.items():
            fo.write(f">{k}\n{s}\n")
    counts = collections.Counter(k.split("|")[0] for k in seqs)
    print("proteins collected:", dict(counts))
    dom = os.path.join(a.out, "homeodomain.domtbl")
    run(["hmmsearch", "--cut_ga", "--domtblout", dom, "-o", os.devnull, a.hmm, allfa])
    best = {}
    for line in open(dom):
        if line[0] == "#":
            continue
        f = line.split()
        e = float(f[12])  # i-Evalue
        if f[0] not in best or e < best[f[0]][0]:
            best[f[0]] = (e, int(f[19]), int(f[20]))  # env from/to
    with open(os.path.join(a.out, "homeodomain_domains.faa"), "w") as fo:
        kept = collections.Counter()
        for k, (e, s0, s1) in best.items():
            fo.write(f">{k}\n{seqs[k][s0 - 1:s1]}\n"); kept[k.split("|")[0]] += 1
    print("with Pfam Homeodomain domain:", dict(kept))
    missing = collections.Counter(k.split("|")[0] for k in seqs if k not in best)
    print("no Pfam Homeodomain hit (kept out of the domain tree):", dict(missing))


def _fasta_full(path):
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


if __name__ == "__main__":
    main()
