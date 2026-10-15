"""Compare native tblastn HSPs with miniprot alignments on the same query and target.

Recall: a tblastn HSP is covered when miniprot has an alignment for the same query
on the same contig that overlaps it. Bit scores (tblastn) are fit against miniprot AS.
"""
import math
import statistics
from collections import defaultdict

TH_BITS = 33.0  # supported_min_bitscore default in the pipeline

tb = []  # tblastn HSPs
with open("tblastn.tsv") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        q, c = f[0], f[1]
        pident, length, qlen, qcov = float(f[2]), int(f[3]), int(f[4]), float(f[5])
        bits = float(f[7])
        s, e = sorted((int(f[8]), int(f[9])))
        tb.append(dict(q=q, c=c, pident=pident, length=length, qlen=qlen, qcov=qcov, bits=bits, s=s, e=e))

mp = []  # miniprot alignments (PAF, 0-based half-open target coords)
with open("miniprot.paf") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        tags = {t.split(":")[0]: t.split(":", 2)[2] for t in f[12:] if t.count(":") >= 2}
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
        tstart, tend = int(f[7]) + 1, int(f[8])  # to 1-based inclusive
        alnlen, nmatch = int(f[10]), int(f[9])
        mp.append(dict(q=f[0], c=f[5], s=min(tstart, tend), e=max(tstart, tend), qlen=qlen,
                       qcov=100.0 * (qe - qs) / qlen, ident=100.0 * nmatch / alnlen,
                       AS=float(tags.get("AS", "nan"))))

by_qc = defaultdict(list)
for m in mp:
    by_qc[(m["q"], m["c"])].append(m)


def overlap(a_s, a_e, b_s, b_e):
    return max(0, min(a_e, b_e) - max(a_s, b_s) + 1)


def covered(h):
    best = 0
    for m in by_qc.get((h["q"], h["c"]), []):
        best = max(best, overlap(h["s"], h["e"], m["s"], m["e"]) / (h["e"] - h["s"] + 1))
    return best


strong = [h for h in tb if h["bits"] >= TH_BITS]
print(f"tblastn HSPs: {len(tb)}; with bits>={TH_BITS}: {len(strong)}")
print(f"miniprot alignments: {len(mp)} across {len({m['q'] for m in mp})} queries")

# HSP-level recall at >=50% overlap of the HSP
cov_frac = [covered(h) for h in strong]
n_cov = sum(1 for x in cov_frac if x >= 0.5)
print(f"HSP recall (>=50% of HSP overlapped by a miniprot alignment, same query+contig): "
      f"{n_cov}/{len(strong)} = {100.0 * n_cov / max(1, len(strong)):.1f}%")
n_any = sum(1 for x in cov_frac if x > 0)
print(f"HSP recall (any overlap): {n_any}/{len(strong)} = {100.0 * n_any / max(1, len(strong)):.1f}%")

# Query-level recall: queries with any strong tblastn HSP vs queries with any miniprot alignment
q_strong = {h["q"] for h in strong}
q_mp = {m["q"] for m in mp}
print(f"query recall: {len(q_strong & q_mp)}/{len(q_strong)} strong-tblastn queries have a miniprot alignment")
print(f"miniprot-only queries (no strong tblastn HSP): {len(q_mp - q_strong)}")

# Pairs with overlap: bits vs AS, identity, coverage
pairs = []
for h in tb:
    best = None
    for m in by_qc.get((h["q"], h["c"]), []):
        o = overlap(h["s"], h["e"], m["s"], m["e"])
        if o / (h["e"] - h["s"] + 1) >= 0.5:
            if best is None or abs(m["AS"] - 0) > -1:
                best = m
    if best is not None:
        pairs.append((h, best))
print(f"overlapping HSP/alignment pairs: {len(pairs)}")
if len(pairs) >= 3:
    xs = [m["AS"] for _, m in pairs]
    ys = [h["bits"] for h, _ in pairs]
    mx, my = statistics.mean(xs), statistics.mean(ys)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    slope = sxy / sxx
    icpt = my - slope * mx
    r = sxy / math.sqrt(sxx * syy)
    print(f"bits ~ AS: slope={slope:.3f} intercept={icpt:.2f} r={r:.3f}")
    # Residual spread near the threshold
    near = [(h["bits"], m["AS"]) for h, m in pairs if 20 <= h["bits"] <= 50]
    print(f"pairs with bits in [20,50]: {len(near)}")
    dp = [abs(h["pident"] - m["ident"]) for h, m in pairs]
    print(f"identity |tblastn pident - miniprot identity|: median {statistics.median(dp):.2f}, max {max(dp):.2f}")
    dc = [abs(h["qcov"] - m["qcov"]) for h, m in pairs]
    print(f"coverage |tblastn qcovhsp - miniprot query coverage|: median {statistics.median(dc):.2f}, max {max(dc):.2f}")

# Uncovered strong HSPs: what are they?
miss = [h for h, x in zip(strong, cov_frac) if x < 0.5]
print(f"uncovered strong HSPs: {len(miss)}")
for h in sorted(miss, key=lambda h: -h["bits"])[:8]:
    print(f"  {h['q'][:40]:40s} {h['c']:24s} {h['s']:>7d}-{h['e']:<7d} bits={h['bits']:.1f} pident={h['pident']:.1f} qcov={h['qcov']:.1f}")
