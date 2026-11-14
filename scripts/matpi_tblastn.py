#!/usr/bin/env python3
"""Genome-level search for matPi (homeodomain) regions and their distance to the matMc-like gene.

Read-only discovery tool. The curated T. deformans matPi is a 181 aa partial protein whose homeodomain region lies in
unannotated sequence, so an annotation-based search (pm_lead_followup.py scan) cannot find it. This script searches the
assembly itself.

For each genome (TSV: name, genome fasta, BFD gff3, matMc gene id or '-'):
  1. tblastn of the matPi query set (FASTA) against the genome (E <= --evalue, no SEG filter);
  2. per contig, hits closer than --merge bp are merged; the best regions by summed bitscore are kept;
  3. the matMc gene's coordinates come from the GFF3; for each matPi region on the matMc contig the distance is given.
Output (--out): matpi_hits.tsv (genome, rank, contig, start, end, best_identity, best_evalue, bitscore, on_matMc_contig,
distance_to_matMc_bp, overlapping_gene).
"""
import argparse
import collections
import gzip
import os
import subprocess
import sys


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def merge_hits(hits, merge):
    """hits: list of (contig, start, end, pident, evalue, bits) -> merged regions sorted by total bitscore.
    Returns list of dict(contig, start, end, ident, evalue, bits)."""
    by = collections.defaultdict(list)
    for c, s, e, p, ev, b in hits:
        by[c].append((min(s, e), max(s, e), p, ev, b))
    regions = []
    for c, v in by.items():
        v.sort()
        cur = None
        for s, e, p, ev, b in v:
            if cur and s - cur["end"] <= merge:
                cur["end"] = max(cur["end"], e)
                cur["ident"] = max(cur["ident"], p)
                cur["evalue"] = min(cur["evalue"], ev)
                cur["bits"] += b
            else:
                if cur:
                    regions.append(cur)
                cur = {"contig": c, "start": s, "end": e, "ident": p, "evalue": ev, "bits": b}
        if cur:
            regions.append(cur)
    return sorted(regions, key=lambda r: -r["bits"])


def distance(a_start, a_end, b_start, b_end):
    """Gap in bp between two intervals; 0 when they overlap."""
    return max(0, max(a_start, b_start) - min(a_end, b_end))


def gene_coords(gff, gene_id):
    if gene_id in ("", "-"):
        return None
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene" and f"ID={gene_id};" in f[8] + ";":
                return f[0], int(f[3]), int(f[4])
    return None


def genes_on(gff, contig):
    out = []
    with opener(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene" and f[0] == contig:
                out.append((int(f[3]), int(f[4]), f[8].split("ID=")[1].split(";")[0]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--merge", type=int, default=2000)
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fo = open(os.path.join(a.out, "matpi_hits.tsv"), "w")
    fo.write("genome\trank\tcontig\tstart\tend\tbest_identity\tbest_evalue\tbitscore\ton_matMc_contig\tdistance_to_matMc_bp\toverlapping_gene\n")
    for line in open(a.genomes):
        name, fa, gff, mc = line.rstrip("\n").split("\t")
        db = os.path.join(a.tmp, "tb_" + name)
        genome = os.path.join(a.tmp, name + ".fa")
        subprocess.run(f"zcat -f {fa} > {genome}", shell=True, check=True)
        subprocess.run(["makeblastdb", "-in", genome, "-dbtype", "nucl", "-out", db, "-logfile", "/dev/null"], check=True)
        out = os.path.join(a.tmp, name + ".tb")
        subprocess.run(["tblastn", "-query", a.queries, "-db", db, "-seg", "no", "-evalue", str(a.evalue), "-out", out,
                        "-outfmt", "6 sseqid sstart send pident evalue bitscore"], check=True)
        hits = []
        for l in open(out):
            c, s, e, p, ev, b = l.split("\t")
            hits.append((c, int(s), int(e), float(p), float(ev), float(b)))
        regions = merge_hits(hits, a.merge)[:a.top]
        mcc = gene_coords(gff, mc)
        for i, r in enumerate(regions, 1):
            on = "NA"
            dist = "NA"
            if mcc:
                on = "yes" if r["contig"] == mcc[0] else "no"
                if on == "yes":
                    dist = distance(r["start"], r["end"], mcc[1], mcc[2])
            gl = [g[2] for g in genes_on(gff, r["contig"]) if g[1] >= r["start"] and g[0] <= r["end"]]
            fo.write(f"{name}\t{i}\t{r['contig']}\t{r['start']}\t{r['end']}\t{r['ident']}\t{r['evalue']}\t{r['bits']}\t{on}\t{dist}\t{','.join(gl) or '-'}\n")
        for f in (genome, out):
            os.remove(f)
    fo.close()
    print(open(os.path.join(a.out, "matpi_hits.tsv")).read())


if __name__ == "__main__":
    sys.exit(main())
