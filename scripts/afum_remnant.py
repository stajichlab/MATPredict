#!/usr/bin/env python3
"""Is the degenerate MAT1-2-1 remnant of the A1163 MAT1-1 locus present, and at the same place, in another assembly?

usage: afum_remnant.py A1163.fasta STRAIN.fasta WORKDIR
For one strain this
  1. aligns the strain to the A1163 assembly with minimap2 (-x asm5 -c) and lifts the A1163 remnant and MAT1-1-1 intervals
     over to the strain (CIGAR walk), keeping only alignments on the A1163 locus scaffold;
  2. runs BLASTN of the A1163 remnant sequence and the A1163 MAT1-1-1 gene against the strain (identity >= 90%);
  3. reports whether the blast hit of the remnant lies in the lifted-over interval ("same location"), how much of the
     1,035-bp remnant is covered, and how far the remnant lies from MAT1-1-1 (A1163: divergent, 117-bp gap);
  4. writes the fraction of 100-bp bins of the A1163 13.3-kb locus covered by alignment blocks (presence strip).
A1163 coordinates (scaffold_3): see tblastn of the database proteins EDP53119.1 and EDP53120.1 (analysis note).
Prints two TSV lines: summary, then strip.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

SCAF = "scaffold_3"
LOCUS = (1653631, 1666890)
REMNANT = (1661652, 1662686)
MAT111 = (1660381, 1661535)
BIN = 100
SLICE = (1540000, 1780000)   # A1163 scaffold_3 slice used as the alignment reference (the whole assembly takes ~7 min per strain)
OFF = SLICE[0] - 1           # PAF target coordinates are slice-relative: add OFF to get scaffold_3 coordinates


def read_scaffold(fa, name):
    seq, on = [], False
    for line in open(fa):
        if line.startswith(">"):
            if on:
                break
            on = line[1:].split()[0] == name
        elif on:
            seq.append(line.strip())
    return "".join(seq)


def lift(paf_lines, interval):
    """Project an A1163 interval onto the strain. Returns (n_mapped_bases, contig, qmin, qmax, strand)."""
    s, e = interval[0] - OFF, interval[1] - OFF   # slice coordinates
    best = None
    for rec in paf_lines:
        f = rec.split("\t")
        if f[5] != SCAF:
            continue
        ts, te = int(f[7]), int(f[8])  # 0-based, end exclusive
        if te < s or ts > e:
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
        if cg is None:
            continue
        strand, qname, qs, qe = f[4], f[0], int(f[2]), int(f[3])
        t, q, mapped, qpos = ts, 0, 0, []
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg):
            n = int(n)
            if op in "M=X":
                lo, hi = max(t, s - 1), min(t + n, e)
                if hi > lo:
                    mapped += hi - lo
                    a, b = q + (lo - t), q + (hi - t)
                    qpos += [a, b]
                t += n; q += n
            elif op in "DN":
                t += n
            elif op == "I":
                q += n
        if mapped and (best is None or mapped > best[0]):
            if strand == "+":
                lo, hi = qs + min(qpos), qs + max(qpos)
            else:
                lo, hi = qe - max(qpos), qe - min(qpos)
            best = (mapped, qname, lo, hi, strand)
    return best


def bins_from_paf(paf_lines):
    n = (LOCUS[1] - LOCUS[0]) // BIN + 1
    cov = [0] * n
    for rec in paf_lines:
        f = rec.split("\t")
        if f[5] != SCAF:
            continue
        ts, te = int(f[7]), int(f[8])
        a, b = max(ts, LOCUS[0] - 1 - OFF), min(te, LOCUS[1] - OFF)
        for i in range(n):
            lo, hi = LOCUS[0] - 1 - OFF + i * BIN, LOCUS[0] - 1 - OFF + (i + 1) * BIN
            ov = min(hi, b) - max(lo, a)
            if ov > 0:
                cov[i] = min(BIN, cov[i] + ov)
    return [c / BIN for c in cov]


def blast_cov(query_fa, db, length):
    out = subprocess.run(["blastn", "-query", query_fa, "-db", db, "-evalue", "1e-10", "-perc_identity", "90",
                          "-outfmt", "6 sseqid pident qstart qend sstart send"], capture_output=True, text=True, check=True).stdout
    iv, hits, best = [], [], 0.0
    for line in out.splitlines():
        sid, pid, qs, qe, ss, se = line.split("\t")
        iv.append((int(qs), int(qe))); hits.append((sid, min(int(ss), int(se)), max(int(ss), int(se)), float(pid)))
    iv.sort(); tot, end = 0, 0
    for a, b in iv:
        a = max(a, end + 1)
        if b >= a:
            tot += b - a + 1; end = max(end, b)
    top = sorted(hits, key=lambda h: -(h[2] - h[1]))[:1]
    return tot / length, top[0] if top else None, hits


def main(a1163, strain, work):
    work = Path(work); work.mkdir(parents=True, exist_ok=True)
    name = Path(strain).name.replace(".sorted.fasta", "")
    scaf = read_scaffold(a1163, SCAF)
    for tag, (s, e) in (("remnant", REMNANT), ("mat111", MAT111)):
        (work / f"{tag}.fa").write_text(f">{tag}\n{scaf[s - 1:e]}\n")
    (work / "ref_slice.fa").write_text(f">{SCAF}\n{scaf[SLICE[0] - 1:SLICE[1]]}\n")
    proc = subprocess.run(["minimap2", "-x", "asm5", "-c", "--secondary=no", "-t", "2", str(work / "ref_slice.fa"), strain], capture_output=True, text=True)
    paf = [l for l in proc.stdout.splitlines() if l.split("\t")[5] == SCAF]
    lr, lm = lift(paf, REMNANT), lift(paf, MAT111)
    db = work / (name + "_db")
    subprocess.run(["makeblastdb", "-in", strain, "-dbtype", "nucl", "-out", str(db)], check=True, capture_output=True)
    rc, rtop, rhits = blast_cov(str(work / "remnant.fa"), str(db), REMNANT[1] - REMNANT[0] + 1)
    mc, mtop, mhits = blast_cov(str(work / "mat111.fa"), str(db), MAT111[1] - MAT111[0] + 1)
    for f in work.glob(db.name + ".n*"):
        f.unlink()
    rlen = REMNANT[1] - REMNANT[0] + 1
    lift_cov = lr[0] / rlen if lr else 0.0
    same = ""
    if lr and rtop:
        same = "yes" if (rtop[0] == lr[1] and rtop[2] >= lr[2] - 1000 and rtop[1] <= lr[3] + 1000) else "no"
    gap = ""
    if rtop and mtop and rtop[0] == mtop[0]:
        gap = str(max(rtop[1], mtop[1]) - min(rtop[2], mtop[2]))
    row = [name, f"{rc:.3f}", f"{rtop[3]:.1f}" if rtop else "", f"{rtop[0]}:{rtop[1]}-{rtop[2]}" if rtop else "",
           f"{mc:.3f}", f"{mtop[0]}:{mtop[1]}-{mtop[2]}" if mtop else "", gap,
           f"{lift_cov:.3f}", f"{lr[1]}:{lr[2]}-{lr[3]}({lr[4]})" if lr else "", same,
           f"{(lm[0] / (MAT111[1] - MAT111[0] + 1)):.3f}" if lm else "0.000"]
    print("\t".join(row))
    print("STRIP\t" + name + "\t" + ",".join(f"{v:.2f}" for v in bins_from_paf(paf)))


if __name__ == "__main__":
    main(*sys.argv[1:4])
