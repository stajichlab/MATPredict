#!/usr/bin/env python3
"""Count long reads by best-matching allele reference (primary alignment) over the HD locus.

Usage: amf_allele_read_count.py BAM MIN_SPAN_FRAC > counts.tsv
A read counts for the reference of its primary alignment when the alignment covers at least MIN_SPAN_FRAC
of that reference. Output: reference, reads, median alignment identity (NM-based).
"""
import collections, statistics, subprocess, sys
bam, frac = sys.argv[1], float(sys.argv[2])
hdr = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True, check=True).stdout
ln = {}
for l in hdr.splitlines():
    if l.startswith("@SQ"):
        f = dict(x.split(":", 1) for x in l.split("\t")[1:])
        ln[f["SN"]] = int(f["LN"])
cnt, ids = collections.Counter(), collections.defaultdict(list)
p = subprocess.Popen(["samtools", "view", "-F", "0x904", bam], stdout=subprocess.PIPE, text=True)
for line in p.stdout:
    f = line.split("\t")
    ref, cigar = f[2], f[5]
    n, rl, m, nm = "", 0, 0, 0
    for ch in cigar:
        if ch.isdigit(): n += ch; continue
        k = int(n); n = ""
        if ch in "MDN=X": rl += k
        if ch in "M=X": m += k
    for t in f[11:]:
        if t.startswith("NM:i:"): nm = int(t[5:])
    if rl >= frac * ln[ref]:
        cnt[ref] += 1; ids[ref].append(1 - nm / max(rl, 1))
print("reference\treads\tmedian_identity")
for r in sorted(ln):
    print(f"{r}\t{cnt[r]}\t{statistics.median(ids[r]):.3f}" if ids[r] else f"{r}\t0\t")
