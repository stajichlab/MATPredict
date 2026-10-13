"""Score a candidate localizer against native tblastn (bits >= 33 HSPs).

usage: evaluate.py <name> <candidate.tsv> <format>
  format: mmseqs  -> query,target,pident,alnlen,qlen,qstart,qend,tstart,tend,evalue,bits
          blast   -> tblastn outfmt 6 used by the pipeline
Coordinates on the target are 1-based; start/end are sorted.
"""
import sys
from collections import defaultdict

TH = 33.0
name, path, fmt = sys.argv[1], sys.argv[2], sys.argv[3]


def rng(a, b):
    return (min(a, b), max(a, b))


tb = []
with open("tblastn.tsv") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        s, e = rng(int(f[8]), int(f[9]))
        tb.append(dict(q=f[0], c=f[1], s=s, e=e, bits=float(f[7])))

cand = []
with open(path) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if fmt == "mmseqs":
            s, e = rng(int(f[7]), int(f[8]))
            cand.append(dict(q=f[0], c=f[1], s=s, e=e, bits=float(f[10]), pident=float(f[2])))
        else:
            s, e = rng(int(f[8]), int(f[9]))
            cand.append(dict(q=f[0], c=f[1], s=s, e=e, bits=float(f[7]), pident=float(f[2])))

by_qc = defaultdict(list)
for h in cand:
    by_qc[(h["q"], h["c"])].append(h)


def frac_covered(h):
    best = 0
    for m in by_qc.get((h["q"], h["c"]), []):
        o = max(0, min(h["e"], m["e"]) - max(h["s"], m["s"]) + 1)
        best = max(best, o / (h["e"] - h["s"] + 1))
    return best


strong = [h for h in tb if h["bits"] >= TH]
q_strong = {h["q"] for h in strong}
hsp_hit = sum(1 for h in strong if frac_covered(h) >= 0.5)
q_hit = len({h["q"] for h in cand} & q_strong)
print(f"[{name}] candidates={len(cand)} queries_with_candidates={len({h['q'] for h in cand})}")
print(f"[{name}] HSP recall (>=50% overlap): {hsp_hit}/{len(strong)} = {100.0*hsp_hit/len(strong):.1f}%")
print(f"[{name}] query recall: {q_hit}/{len(q_strong)} = {100.0*q_hit/len(q_strong):.1f}%")
