#!/usr/bin/env python3
"""Build an HD1-like profile HMM from AMF homologs and search proteomes with it.

Read-only discovery tool. Needs stage-1 output (flank_first_scan.py) and the Ropars 2016 control FASTA.

build  : training set = unique sequences of (a) Ropars 2016 proteins titled HD1-like and (b) HD1-class genes
         (control hit E <= --train-evalue) from the Glomeromycota genomes in the stage-1 output.
         mafft --localpair alignment; columns with < --occupancy residues removed; hmmbuild.
search : hmmsearch of the HMM against every proteome in a genome list. For each hit the table gives
         score, E-value, contig, gene index, strand, and the nearest Pfam Homeodomain gene
         (genes between, strand relation) from the stage-1 candidates.

Usage: hd1_hmm.py build  --base BASE --out OUT
       hd1_hmm.py search --base BASE --out OUT --lists LIST [LIST ...]
"""
import argparse
import collections
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homeodomain_tree import Prot, table  # noqa: E402
from flank_first_scan import run  # noqa: E402

ANN = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation"


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


def build(a):
    os.makedirs(a.out, exist_ok=True)
    D = f"{a.base}/stage1_all/per_genome"
    P = Prot()
    train = {}
    for h, s in fasta_full(f"{a.base}/ropars/ropars_proteins.faa"):
        if re.search(r"hd1", h, re.I):
            train[f"ropars|{h.split()[0]}"] = s
    for g in [x.strip() for x in open(f"{a.base}/list_glomeromycota.txt") if x.strip()]:
        for c in table(f"{D}/{g}.candidates.tsv"):
            if c[5] == "HD1" and float(c[7]) <= a.train_evalue:
                s = P.get(g).get(c[0])
                if s:
                    train[f"{g}|{c[0]}"] = s
    uniq, seen = {}, {}
    for k, s in train.items():
        if s in seen:
            continue
        seen[s] = k; uniq[k] = s
    print(f"training proteins: {len(train)} total, {len(uniq)} unique sequences")
    fa = os.path.join(a.out, "hd1_train.faa")
    with open(fa, "w") as fo:
        for k, s in uniq.items():
            fo.write(f">{k}\n{s}\n")
    r = subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", str(a.threads), "--quiet", fa],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(r.stderr[-400:])
    aln, k = {}, None
    for line in r.stdout.splitlines():
        if line.startswith(">"):
            k = line[1:].strip(); aln[k] = []
        else:
            aln[k].append(line.strip())
    aln = {k: "".join(v).upper() for k, v in aln.items()}
    L = len(next(iter(aln.values())))
    keep = [i for i in range(L) if sum(s[i] != "-" for s in aln.values()) >= a.occupancy * len(aln)]
    with open(os.path.join(a.out, "hd1_train_trim.afa"), "w") as fo:
        for k, s in aln.items():
            fo.write(f">{k}\n{''.join(s[i] for i in keep)}\n")
    print(f"alignment columns {L}, kept {len(keep)} (occupancy >= {a.occupancy})")
    run(["hmmbuild", "--amino", "-n", "AMF_HD1like", os.path.join(a.out, "AMF_HD1like.hmm"),
         os.path.join(a.out, "hd1_train_trim.afa")])
    print("wrote", os.path.join(a.out, "AMF_HD1like.hmm"))


def search(a):
    hmm = os.path.join(a.out, "AMF_HD1like.hmm")
    tmp = os.environ.get("SCRATCH", "/tmp")
    D = f"{a.base}/stage1_all/per_genome"
    out = open(os.path.join(a.out, "hd1_hmm_hits.tsv"), "w")
    out.write("genome\tgene\tscore\tevalue\tcontig\tindex\tstrand\tnearest_homeodomain_gene\tgenes_between\tstrand_relation\n")
    n_done = 0
    for lst in a.lists:
        for g in [x.strip() for x in open(lst) if x.strip()]:
            p = f"{ANN}/{g}/predict_results/{g}.proteins.fa.gz"
            if not os.path.exists(p) or not os.path.exists(f"{D}/{g}.genes_all.tsv"):
                continue
            tbl = os.path.join(tmp, "hd1.tbl")
            r = subprocess.run(f"zcat {p} | hmmsearch --cpu 2 -E {a.evalue} --tblout {tbl} -o /dev/null {hmm} -",
                               shell=True, capture_output=True, text=True)
            if r.returncode != 0:
                print("FAIL", g, r.stderr[-200:]); continue
            hits = []
            for line in open(tbl):
                if line[0] == "#":
                    continue
                f = line.split()
                hits.append((re.sub(r"-T\d+$", "", f[0]), float(f[5]), float(f[4])))
            if not hits:
                n_done += 1; continue
            st = {r_[0]: (r_[1], r_[4], int(r_[5])) for r_ in table(f"{D}/{g}.genes_all.tsv")}
            hom = [(c[0], st[c[0]]) for c in table(f"{D}/{g}.candidates.tsv") if c[3] == "Homeodomain" and c[0] in st]
            for gene, score, e in hits:
                if gene not in st:
                    continue
                ctg, strand, idx = st[gene]
                best = None
                for h, (c2, s2, i2) in hom:
                    if c2 == ctg and h != gene and (best is None or abs(i2 - idx) < abs(best[1][2] - idx)):
                        best = (h, (c2, s2, i2))
                if best and abs(best[1][2] - idx) <= a.near:
                    d = abs(best[1][2] - idx) - 1
                    left = (strand, best[1][1]) if idx < best[1][2] else (best[1][1], strand)
                    rel = "divergent" if left == ("-", "+") else "convergent" if left == ("+", "-") else "same"
                    out.write(f"{g}\t{gene}\t{score}\t{e}\t{ctg}\t{idx}\t{strand}\t{best[0]}\t{d}\t{rel}\n")
                else:
                    out.write(f"{g}\t{gene}\t{score}\t{e}\t{ctg}\t{idx}\t{strand}\t\t\t\n")
            n_done += 1
    out.close()
    print("genomes searched:", n_done)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["build", "search"])
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--lists", nargs="*", default=[])
    ap.add_argument("--train-evalue", type=float, default=1e-20)
    ap.add_argument("--occupancy", type=float, default=0.5)
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--near", type=int, default=5)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    build(a) if a.cmd == "build" else search(a)


if __name__ == "__main__":
    main()
