#!/usr/bin/env python3
"""Transcripts that overlap a locus, with exon structure, strand and Pfam hits on every ORF of the spliced sequence.

Read-only analysis of a StringTie assembly (GTF) and the spliced transcript FASTA from `gffread -w`.
For every transcript on --contig that overlaps [--lo, --hi]: number of exons, span, strand, spliced length, the longest ORF,
and the Pfam domains (HMG_box, HMG_box_2, Homeodomain; GA cutoff) found on any ORF of >= --min-aa residues in the six frames of
the spliced sequence, with the frame (negative = reverse of the transcript as written).
Limit: with unstranded RNA-seq StringTie takes the strand from splice-site motifs, so a homeodomain ORF in a negative frame of a
'+' transcript (or the reverse) is on the strand opposite to the assembled transcript.
Usage: locus_transcripts.py --gtf X.gtf --tx X.tx.fa --contig C --lo N --hi N --hmm HMM
"""
import argparse
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from region_pfam_check import orfs  # noqa: E402


def read_fasta(path):
    d, h = {}, None
    for line in open(path):
        if line[0] == ">":
            h = line[1:].split()[0]
            d[h] = []
        else:
            d[h].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}


def read_gtf_transcripts(path, contig):
    """-> {transcript_id: {'strand': s, 'exons': [(start, end)]}} for one contig."""
    recs = {}
    for line in open(path):
        f = line.split("\t")
        if len(f) < 9 or f[0] != contig:
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        r = recs.setdefault(m.group(1), {"strand": f[6], "exons": []})
        if f[2] == "exon":
            r["exons"].append((int(f[3]), int(f[4])))
    return recs


def overlaps(exons, lo, hi):
    ex = sorted(exons)
    return bool(ex) and ex[-1][1] >= lo and ex[0][0] <= hi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--tx", required=True)
    ap.add_argument("--contig", required=True)
    ap.add_argument("--lo", type=int, required=True)
    ap.add_argument("--hi", type=int, required=True)
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--min-aa", type=int, default=30)
    ap.add_argument("--tmp", default=os.environ.get("SCRATCH", "/tmp"))
    a = ap.parse_args()
    tx = read_fasta(a.tx)
    print("transcript\tstrand\texons\tspan\ttx_len\tlongest_orf_aa\tpfam_hits(aa,frame,domain,evalue)")
    for tid, r in read_gtf_transcripts(a.gtf, a.contig).items():
        if not overlaps(r["exons"], a.lo, a.hi) or tid not in tx:
            continue
        ex = sorted(r["exons"])
        os_ = orfs(tx[tid], a.min_aa)
        faa = os.path.join(a.tmp, "lt.faa")
        with open(faa, "w") as fh:
            for i, o in enumerate(os_):
                fh.write(f">o{i}\n{o[3]}\n")
        tbl = os.path.join(a.tmp, "lt.tbl")
        subprocess.run(["hmmsearch", "--cut_ga", "--domtblout", tbl, "-o", "/dev/null", a.hmm, faa], check=True)
        hits = []
        for line in open(tbl):
            if line[0] != "#":
                f = line.split()
                o = os_[int(f[0][1:])]
                hits.append((len(o[3]), o[0], f[3], f[6]))
        longest = max((len(o[3]) for o in os_), default=0)
        print(f"{tid}\t{r['strand']}\t{len(ex)}\t{ex[0][0]}-{ex[-1][1]}\t{len(tx[tid])}\t{longest}\t{hits or 'none'}")


if __name__ == "__main__":
    main()
