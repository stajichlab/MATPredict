#!/usr/bin/env python3
"""Cut a MAT1-1 / MAT1-2 record pair the way the repo's panels were cut: core = the gap between the two long colinear
blastn blocks that flank the idiomorph, plus up to FLANK bp of flank each side, kept only while the two source loci align
there at >= MINID identity.  Coordinates 1-based inclusive.
usage: build_record.py REF.fa(idiomorph A locus window, id used as-is) REFNAME CONTIG.fa CONTIGNAME OUT_A.fasta OUT_B.fasta A_label B_label
Both outputs are written on the forward strand of REF; CONTIG is reverse-complemented when it aligns on the minus strand."""
import subprocess, sys, tempfile, os, re
FLANK, MINID, MINBLOCK = 600, 97.0, 500
OFFSET = int(os.environ.get("REF_OFFSET", "0"))   # 1-based start of the REF window minus 1, to report coordinates on the source contig
def read(p):
    d = {}; n = None
    for l in open(p):
        if l[0] == ">": n = l[1:].split()[0]; d[n] = []
        else: d[n].append(l.strip().upper())
    return {k: "".join(v) for k, v in d.items()}
def rc(s): return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
ref_fa, refname, ctg_fa, ctgname, outA, outB, labA, labB = sys.argv[1:9]
ref = read(ref_fa); ctg = read(ctg_fa)
rseq = next(iter(ref.values())); cseq = ctg[ctgname]
out = subprocess.run(["blastn", "-query", ref_fa, "-subject", ctg_fa, "-evalue", "1e-20", "-outfmt", "6 qstart qend sstart send pident length qseq sseq"],
                     capture_output=True, text=True, check=True).stdout.strip().splitlines()
H = []
for l in out:
    f = l.split("\t"); qs, qe, ss, se = map(int, f[:4]); H.append(dict(qs=qs, qe=qe, ss=ss, se=se, pid=float(f[4]), len=int(f[5]), qseq=f[6], sseq=f[7]))
H = [h for h in H if h["len"] >= MINBLOCK]
if not H: sys.exit("no long blocks")
minus = sum(h["len"] for h in H if h["ss"] > h["se"]) > sum(h["len"] for h in H if h["ss"] < h["se"])
L = len(cseq)
if minus:
    cseq = rc(cseq)
    for h in H: h["ss"], h["se"] = L - h["ss"] + 1, L - h["se"] + 1
H = [h for h in H if h["ss"] < h["se"]]; H.sort(key=lambda h: h["qs"])
# adjacent blocks with the largest gap on the REF side = the idiomorph
best = max(range(len(H) - 1), key=lambda i: H[i+1]["qs"] - H[i]["qe"]) if len(H) > 1 else None
A, B = H[best], H[best+1]
refA = (A["qe"] + 1, B["qs"] - 1); ctgB = (A["se"] + 1, B["ss"] - 1)
def flank_ok(h, side):
    """longest flank <= FLANK (counted in ref bases from the core edge) whose alignment identity >= MINID"""
    q, s = h["qseq"], h["sseq"]
    if side == "L": q, s = q[::-1], s[::-1]   # walk away from the core: left flank is the right end of block A
    pos = []   # (ref bases consumed, matches so far) after each alignment column that consumes a ref base
    n = ident = 0
    for x, y in zip(q, s):
        if x != "-": n += 1
        if x == y: ident += 1
        if x != "-": pos.append((n, ident))
        if n >= FLANK: break
    best = 0
    for n in range(min(FLANK, len(pos)), 0, -25):     # longest prefix (steps of 25 bp) with identity >= MINID
        if pos[n-1][1] / n * 100 >= MINID: best = n; break
    return best
# left flank is the end of block A (nearest the core), right flank is the start of block B
lf = flank_ok(A, "L"); rf = flank_ok(B, "R")
rf_len, lf_len = min(rf, FLANK), min(lf, FLANK)
def cut(seq, a, b, nm, coords_offset=0):
    return seq[a-1-lf_len: b+rf_len], (a - lf_len, b + rf_len)
recA, ca = cut(rseq, *refA, refname); recB, cb = cut(cseq, *ctgB, ctgname)
print(f"blocks A q{A['qs']}-{A['qe']} s{A['ss']}-{A['se']} id{A['pid']}; B q{B['qs']}-{B['qe']} s{B['ss']}-{B['se']} id{B['pid']}; minus={minus}", file=sys.stderr)
print(f"{labA} core ref {refA} ({refA[1]-refA[0]+1} bp); {labB} core contig {ctgB} ({ctgB[1]-ctgB[0]+1} bp); flanks L{lf_len} R{rf_len}", file=sys.stderr)
rid = refname if ":" in refname else refname
strand = "(-)" if minus else ""
open(outA, "w").write(f">{refname}:{ca[0]+OFFSET}-{ca[1]+OFFSET}\n{recA}\n")
open(outB, "w").write(f">{ctgname}:{cb[0]}-{cb[1]}{strand}\n{recB}\n")
print(len(recA), len(recB), file=sys.stderr)
