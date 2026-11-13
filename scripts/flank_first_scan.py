#!/usr/bin/env python3
"""Flank-first scan, stage 1: candidate genes and gene neighbourhoods for one genome.

Read-only discovery tool. It does not use or change the MATPredict detector.

Inputs (per genome): funannotate proteins FASTA and GFF3 (plain or .gz).
Steps:
 1. hmmsearch (Pfam GA cutoff) for HMG_box / HMG_box_2 / Homeodomain.
 2. diamond blastp of the proteome against a control protein set. The class of a
    control comes from its FASTA header (HD1, HD2, choline, PGAM).
 3. Gene order from the GFF3. A candidate is any gene with a Pfam hit or a control hit.
 4. Candidates within WINDOW genes of each other on one contig form a locus.
 5. For each locus, write the flanking genes (WINDOW on each side).

Outputs in --out, prefix <name>:
  <name>.genes_all.tsv   gene, contig, start, end, strand, index (all genes, in order)
  <name>.candidates.tsv  gene, contig, index, pfam, pfam_evalue, control_class, control_id, control_evalue
  <name>.loci.tsv        locus id, contig, first/last index, genes, classes, flags
  <name>.flanks.tsv      locus id, gene, offset (negative = upstream index), in_locus
  <name>.flank_proteins.faa  proteins of the flanking genes, id = <name>|<gene>
"""
import argparse
import collections
import gzip
import os
import re
import subprocess
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    name, seq = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def parse_gff(path):
    """Return ordered gene list [(gene, contig, start, end, strand)] and mRNA->gene map."""
    genes, m2g = {}, {}
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "gene":
                genes[attr["ID"]] = (attr["ID"], f[0], int(f[3]), int(f[4]), f[6])
            elif f[2] == "mRNA":
                m2g[attr["ID"]] = attr.get("Parent", attr["ID"])
    order = sorted(genes.values(), key=lambda g: (g[1], g[2]))
    return order, m2g


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"command failed ({r.returncode}): {' '.join(cmd)}\n{r.stderr[-800:]}")
    return r


def control_class(header):
    h = header.lower()
    if "|" in header:  # curated Mucoromycota MAT set: >id|gene|class
        return header.strip().split("|")[-1].split()[0]
    if "hd1" in h:
        return "HD1"
    if "hd2" in h:
        return "HD2"
    if "choline" in h:
        return "choline"
    if "phosphoglycerate" in h:
        return "PGAM"
    return "other"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--pfam-hmm", required=True, help="HMM file with the Pfam models to search")
    ap.add_argument("--controls-dmnd", required=True, help="diamond db of control proteins")
    ap.add_argument("--controls-faa", required=True, help="control FASTA (for class names)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=5, help="genes on each side (default 5)")
    ap.add_argument("--control-evalue", type=float, default=1e-5)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    tmp = os.path.join(a.tmp, f"ffs_{a.name}")
    os.makedirs(tmp, exist_ok=True)
    pre = os.path.join(a.out, a.name)

    order, m2g = parse_gff(a.gff3)
    idx = {g[0]: i for i, g in enumerate(order)}
    contig_of = {g[0]: g[1] for g in order}

    prot = {}
    for pid, seq in read_fasta(a.proteins):
        gene = m2g.get(pid, re.sub(r"-T\d+$", "", pid))
        if gene not in idx:
            continue
        if gene not in prot or len(seq) > len(prot[gene]):
            prot[gene] = seq.rstrip("*")
    faa = os.path.join(tmp, "proteome.faa")
    with open(faa, "w") as fo:
        for g, s in prot.items():
            fo.write(f">{g}\n{s}\n")

    with open(pre + ".genes_all.tsv", "w") as fo:
        fo.write("gene\tcontig\tstart\tend\tstrand\tindex\n")
        for i, g in enumerate(order):
            fo.write("\t".join(map(str, (*g, i))) + "\n")

    # 1. Pfam
    dom = os.path.join(tmp, "dom.tbl")
    run(["hmmsearch", "--cut_ga", "--cpu", "1", "--domtblout", dom, "-o", os.devnull, a.pfam_hmm, faa])
    pfam = {}
    with open(dom) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split()
            g, d, e = f[0], f[3], float(f[6])
            if g not in pfam or e < pfam[g][1]:
                pfam[g] = (d, e)

    # 2. control homologs
    ctl_name = {n: control_class(open_h) for n, open_h in _headers(a.controls_faa).items()}
    tab = os.path.join(tmp, "ctl.tsv")
    run(["diamond", "blastp", "-q", faa, "-d", a.controls_dmnd, "--evalue", str(a.control_evalue),
         "--max-target-seqs", "5", "--threads", str(a.threads), "--outfmt", "6", "qseqid", "sseqid",
         "evalue", "bitscore", "-o", tab, "--tmpdir", tmp, "--quiet"])
    ctl = {}
    with open(tab) as fh:
        for line in fh:
            q, s, e, b = line.rstrip("\n").split("\t")
            b = float(b)
            if q not in ctl or b > ctl[q][3]:
                ctl[q] = (ctl_name.get(s, "other"), s, float(e), b)

    cands = sorted(set(pfam) | set(ctl), key=lambda g: idx[g])
    with open(pre + ".candidates.tsv", "w") as fo:
        fo.write("gene\tcontig\tindex\tpfam\tpfam_evalue\tcontrol_class\tcontrol_id\tcontrol_evalue\n")
        for g in cands:
            p = pfam.get(g, ("", ""))
            c = ctl.get(g, ("", "", ""))
            fo.write(f"{g}\t{contig_of[g]}\t{idx[g]}\t{p[0]}\t{p[1]}\t{c[0]}\t{c[1]}\t{c[2]}\n")

    # 3. loci: candidates within window genes on one contig
    loci, cur = [], []
    for g in cands:
        if cur and contig_of[g] == contig_of[cur[-1]] and idx[g] - idx[cur[-1]] <= a.window:
            cur.append(g)
        else:
            if cur:
                loci.append(cur)
            cur = [g]
    if cur:
        loci.append(cur)

    with open(pre + ".loci.tsv", "w") as fl, open(pre + ".flanks.tsv", "w") as ff, \
            open(pre + ".flank_proteins.faa", "w") as fp:
        fl.write("locus\tcontig\tfirst_index\tlast_index\tn_genes\tclasses\tpfam\tflag\n")
        ff.write("locus\tgene\toffset\tin_locus\n")
        written = set()
        for n, L in enumerate(loci, 1):
            lid = f"{a.name}_L{n:04d}"
            cls = sorted({ctl[g][0] for g in L if g in ctl})
            pf = sorted({pfam[g][0] for g in L if g in pfam})
            has_mat = bool({"HD1", "HD2"} & set(cls)) or bool(pf)
            has_flank = bool({"choline", "PGAM"} & set(cls))
            flag = "mat_and_flank" if (has_mat and has_flank) else ("mat_only" if has_mat else "flank_only")
            first, last = idx[L[0]], idx[L[-1]]
            fl.write(f"{lid}\t{contig_of[L[0]]}\t{first}\t{last}\t{len(L)}\t{','.join(cls)}\t{','.join(pf)}\t{flag}\n")
            lo, hi = max(0, first - a.window), min(len(order) - 1, last + a.window)
            for i in range(lo, hi + 1):
                g = order[i][0]
                if order[i][1] != contig_of[L[0]]:
                    continue
                off = i - first if i < first else (i - last if i > last else 0)
                ff.write(f"{lid}\t{g}\t{off}\t{int(g in L)}\n")
                if g in prot and (a.name, g) not in written:
                    written.add((a.name, g))
                    fp.write(f">{a.name}|{g}\n{prot[g]}\n")
    print(f"{a.name}: genes={len(order)} candidates={len(cands)} loci={len(loci)}")


def _headers(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                d[line[1:].split()[0]] = line[1:]
    return d


if __name__ == "__main__":
    main()
