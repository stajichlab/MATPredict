#!/usr/bin/env python3
"""Cut a MAT1-1 / MAT1-2 record pair the way the repo's panels were cut: core = the gap between the two long colinear
blastn blocks that flank the idiomorph, plus up to FLANK bp of flank each side, kept only while the two source loci align
there at >= MINID identity.  Coordinates 1-based inclusive.
usage: build_record.py REF.fa(idiomorph A locus window, id used as-is) REFNAME CONTIG.fa CONTIGNAME OUT_A.fasta OUT_B.fasta A_label B_label
Both outputs are written on the forward strand of REF; CONTIG is reverse-complemented when it aligns on the minus strand."""
import subprocess, sys, tempfile, os, re
FLANK, MINID, MINBLOCK = 600, 97.0, 500
PFX = os.environ.get("CTG_PREFIX", "")
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
def conserved_edge(h, side):
    """Walk from the core edge into block h (side 'L': block A, from its end backwards; 'R': block B, from its start).
    The columns right next to the core are often diverged (the alignment extends a little into the idiomorph), so first find
    the first column c0 where a 100-column window is >= MINID identical; everything before c0 stays in the core.
    Returns (ref_bp_trimmed, ctg_bp_trimmed, flank_len): bases between the block edge and c0 on each sequence, and the
    longest flank (<= FLANK ref bases, 25-bp steps) from c0 with identity >= MINID."""
    q, s_ = h["qseq"], h["sseq"]
    if side == "L": q, s_ = q[::-1], s_[::-1]
    n = len(q)
    match = [1 if (x == y and x != "-") else 0 for x, y in zip(q, s_)]
    c0 = None
    for c in range(0, n - 100):
        if sum(match[c:c+100]) >= MINID:
            c0 = c; break
    if c0 is None: return 0, 0, 0
    trim_q = sum(1 for x in q[:c0] if x != "-"); trim_s = sum(1 for y in s_[:c0] if y != "-")
    nq = ident = 0; best = 0; steps = {}
    for x, y, m in zip(q[c0:], s_[c0:], match[c0:]):
        if x != "-":
            nq += 1
            ident += m
            if nq <= FLANK and nq % 25 == 0 and ident / nq * 100 >= MINID: best = nq
        else:
            ident += 0
        if nq >= FLANK: break
    return trim_q, trim_s, best
tqA, tsA, lf = conserved_edge(A, "L"); tqB, tsB, rf = conserved_edge(B, "R")
lf_len, rf_len = min(lf, FLANK), min(rf, FLANK)
refA = (A["qe"] - tqA + 1, B["qs"] + tqB - 1); ctgB = (A["se"] - tsA + 1, B["ss"] + tsB - 1)
def cut(seq, a, b):
    return seq[a-1-lf_len: b+rf_len], (a - lf_len, b + rf_len)
recA, ca = cut(rseq, *refA); recB, cb = cut(cseq, *ctgB)
print(f"blocks A q{A['qs']}-{A['qe']} s{A['ss']}-{A['se']} id{A['pid']}; B q{B['qs']}-{B['qe']} s{B['ss']}-{B['se']} id{B['pid']}; minus={minus}", file=sys.stderr)
print(f"{labA} core ref {refA} ({refA[1]-refA[0]+1} bp); {labB} core contig {ctgB} ({ctgB[1]-ctgB[0]+1} bp); flanks L{lf_len} R{rf_len}", file=sys.stderr)
rid = refname if ":" in refname else refname
strand = "(-)" if minus else ""
open(outA, "w").write(f">{refname}:{ca[0]+OFFSET}-{ca[1]+OFFSET}\n{recA}\n")
open(outB, "w").write(f">{PFX}{ctgname}:{cb[0]}-{cb[1]}{strand}\n{recB}\n")
print(len(recA), len(recB), file=sys.stderr)
