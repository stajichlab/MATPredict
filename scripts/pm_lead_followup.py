#!/usr/bin/env python3
"""Structural follow-up of PM-family detector leads, and a rescan of genomes with no call.

Read-only discovery tool. It does not change the detector.

lead   : for every `PM` locus in detection_report.yaml files, take the gene_evidence entries for matMc / matPi
         (contig, start, end, strand), find the predicted gene (BFD funannotate GFF3) that overlaps each, test that
         protein with Pfam HMG_box / HMG_box_2 / Homeodomain (hmmsearch, GA cutoff), and write the gene order of
         the locus +/- --flank genes with strands. Output: lead_genes.tsv, lead_summary.tsv.
scan   : for genomes with no PM call, run hmmsearch (HMG_box, HMG_box_2, Homeodomain) over the proteome and diamond blastp
         against the curated PM proteins (matMc, matPi, matMi) and the curated S. pombe mat1 proteins. Output: scan_genes.tsv,
         scan_pairs.tsv (candidate genes within --near genes on one contig).
Usage: pm_lead_followup.py lead --runs RUNS_DIR --genomes LIST --out OUT --hmm HMM
       pm_lead_followup.py scan --genomes LIST --out OUT --hmm HMM --controls FAA
Pfam and diamond are needed on PATH (module hmmer/3.4, diamond).
"""
import argparse
import collections
import gzip
import os
import re
import subprocess
import sys

ANN = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation"
DOMAINS = ("HMG_box", "HMG_box_2", "Homeodomain")


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_gff_genes(gff):
    """-> {contig: [(start, end, strand, gene_id)]} sorted by start."""
    out = collections.defaultdict(list)
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"ID=([^;]+)", f[8])
            out[f[0]].append((int(f[3]), int(f[4]), f[6], m.group(1) if m else ""))
    for v in out.values():
        v.sort()
    return out


def overlapping_gene(genes, contig, s, e):
    """Gene with the largest overlap with [s, e] on contig -> (index, gene tuple) or (None, None)."""
    best, bi = 0, None
    for i, g in enumerate(genes.get(contig, [])):
        ov = min(e, g[1]) - max(s, g[0]) + 1
        if ov > best:
            best, bi = ov, i
    return (bi, genes[contig][bi]) if bi is not None else (None, None)


def hmm_domains(proteome, hmm, tmp, tag):
    """{protein_id: [(domain, evalue)]} from hmmsearch --cut_ga."""
    tbl = os.path.join(tmp, tag + ".tbl")
    faa = os.path.join(tmp, tag + ".faa")
    subprocess.run(f"zcat -f {proteome} > {faa}", shell=True, check=True)
    # hmmsearch cannot rewind a pipe when the HMM file holds several models, so search a plain file
    subprocess.run(["hmmsearch", "--cut_ga", "--cpu", "2", "--domtblout", tbl, "-o", "/dev/null", hmm, faa], check=True)
    out = collections.defaultdict(list)
    for line in open(tbl):
        if line[0] == "#":
            continue
        f = line.split()
        out[re.sub(r"-T\d+$", "", f[0])].append((f[3], float(f[6])))
    return out


def layout_row(rows_idx, genes_c, hits, flank):
    """Gene order text for the genes spanning the evidence indices +/- flank."""
    lo, hi = max(0, min(rows_idx) - flank), min(len(genes_c) - 1, max(rows_idx) + flank)
    out = []
    for i in range(lo, hi + 1):
        g = genes_c[i]
        dom = ",".join(sorted({d for d, _ in hits.get(g[3], []) if d in DOMAINS}))
        out.append((i, g[3], g[2], g[0], g[1], dom))
    return out


def do_lead(a):
    import yaml
    os.makedirs(a.out, exist_ok=True)
    tmp = os.environ.get("SCRATCH", "/tmp")
    fo = open(os.path.join(a.out, "lead_genes.tsv"), "w")
    fs = open(os.path.join(a.out, "lead_summary.tsv"), "w")
    fo.write("genome\tlocus\tgene_index\tgene\tstrand\tstart\tend\tdomain\tevidence\n")
    fs.write("genome\tcontig\tmatMc_gene\tmatMc_domain\tmatPi_gene\tmatPi_domain\tgenes_between\tstrand_relation\n")
    for line in open(a.genomes):
        asm, _, name = line.rstrip("\n").partition("\t")
        rep = os.path.join(a.runs, asm, "detection_report.yaml")
        if not os.path.exists(rep):
            continue
        d = yaml.safe_load(open(rep))
        pm = [x for x in (d.get("detected") or []) if x["family"].endswith(":PM")]
        if not pm:
            continue
        P = f"{ANN}/{name}/predict_results"
        genes = read_gff_genes(f"{P}/{name}.gff3.gz")
        hits = hmm_domains(f"{P}/{name}.proteins.fa.gz", a.hmm, tmp, "lead")
        for li, locus in enumerate(pm):
            contig = locus["contig"]
            ev, idx = {}, []
            for e in locus["gene_evidence"]:
                if e["gene"] in ("matMc", "matPi", "matMi"):
                    i, g = overlapping_gene(genes, contig, e["start"], e["end"])
                    if g:
                        ev[e["gene"]] = (i, g, f"{e['identity']}% id, {e['coverage']}% cov")
                        idx.append(i)
            if not idx:
                continue
            for r in layout_row(idx, genes[contig], hits, a.flank):
                tag = ";".join(f"{k}:{v[2]}" for k, v in ev.items() if v[0] == r[0])
                fo.write(f"{name}\t{contig}#{li}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{tag}\n")

            def dom(k):
                return ",".join(sorted({x for x, _ in hits.get(ev[k][1][3], []) if x in DOMAINS})) if k in ev else "NA"
            mc, pi = ev.get("matMc"), ev.get("matPi")
            btw = rel = "NA"
            if mc and pi and mc[0] != pi[0]:
                btw = abs(mc[0] - pi[0]) - 1
                left, right = (mc, pi) if mc[0] < pi[0] else (pi, mc)
                rel = {"+-": "convergent", "-+": "divergent"}.get(left[1][2] + right[1][2], "same")
            elif mc and pi:
                btw, rel = "same gene", "fused"
            fs.write(f"{name}\t{contig}\t{mc[1][3] if mc else 'NA'}\t{dom('matMc')}\t{pi[1][3] if pi else 'NA'}\t{dom('matPi')}\t{btw}\t{rel}\n")
    fo.close(); fs.close()
    print(open(os.path.join(a.out, "lead_summary.tsv")).read())


def do_scan(a):
    os.makedirs(a.out, exist_ok=True)
    tmp = os.environ.get("SCRATCH", "/tmp")
    subprocess.run(["diamond", "makedb", "--in", a.controls, "-d", os.path.join(tmp, "ctl"), "--quiet"], check=True)
    fg = open(os.path.join(a.out, "scan_genes.tsv"), "w")
    fp = open(os.path.join(a.out, "scan_pairs.tsv"), "w")
    fg.write("genome\tgene\tcontig\tindex\tstrand\tdomain\tbest_control\tcontrol_evalue\tcontrol_identity\n")
    fp.write("genome\tgene_a\tgene_b\tcontig\tgenes_between\tstrand_relation\tdomains\n")
    for line in open(a.genomes):
        name = line.strip().split("\t")[-1]
        if not name:
            continue
        P = f"{ANN}/{name}/predict_results"
        genes = read_gff_genes(f"{P}/{name}.gff3.gz")
        order = {}
        for c, v in genes.items():
            for i, g in enumerate(v):
                order[g[3]] = (c, i, g[2])
        hits = hmm_domains(f"{P}/{name}.proteins.fa.gz", a.hmm, tmp, "scan")
        out = os.path.join(tmp, "dm.tsv")
        subprocess.run(f"zcat -f {P}/{name}.proteins.fa.gz | diamond blastp -d {os.path.join(tmp, 'ctl')} -o {out} --quiet "
                       f"--max-target-seqs 1 -e 1e-3 --outfmt 6 qseqid sseqid evalue pident --threads 2 -q /dev/stdin",
                       shell=True, check=True)
        best = {}
        for l in open(out):
            f = l.split("\t")
            best.setdefault(re.sub(r"-T\d+$", "", f[0]), (f[1].split("|")[-1] if "name=" not in f[1] else re.search(r"name=([^|]+)", f[1]).group(1), f[2], f[3].strip()))
        cand = []
        for g, hs in hits.items():
            dm = ",".join(sorted({d for d, _ in hs if d in DOMAINS}))
            if dm and g in order:
                b = best.get(g, ("", "", ""))
                c, i, s = order[g]
                cand.append((c, i, s, g, dm))
                fg.write(f"{name}\t{g}\t{c}\t{i}\t{s}\t{dm}\t{b[0]}\t{b[1]}\t{b[2]}\n")
        for g in best:
            if g in order and g not in {x[3] for x in cand}:
                c, i, s = order[g]
                fg.write(f"{name}\t{g}\t{c}\t{i}\t{s}\t\t{best[g][0]}\t{best[g][1]}\t{best[g][2]}\n")
        cand.sort()
        for x, y in zip(cand, cand[1:]):
            if x[0] == y[0] and y[1] - x[1] <= a.near:
                rel = {"+-": "convergent", "-+": "divergent"}.get(x[2] + y[2], "same")
                fp.write(f"{name}\t{x[3]}\t{y[3]}\t{x[0]}\t{y[1] - x[1] - 1}\t{rel}\t{x[4]}|{y[4]}\n")
    fg.close(); fp.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["lead", "scan"])
    ap.add_argument("--runs")
    ap.add_argument("--genomes", required=True, help="TSV: ASMID<TAB>annotation folder name")
    ap.add_argument("--out", required=True)
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--controls")
    ap.add_argument("--flank", type=int, default=5)
    ap.add_argument("--near", type=int, default=5)
    a = ap.parse_args()
    do_lead(a) if a.cmd == "lead" else do_scan(a)


if __name__ == "__main__":
    sys.exit(main())
