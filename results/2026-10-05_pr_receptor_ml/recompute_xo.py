#!/usr/bin/env python3
"""Cross-order versions of the precursor-homology and MAT-gene (HD, flank) features.

build_genome.py leaves out curated queries of the genome's own species. A new order has no
curated records at all, so here every curated query from the genome's own ORDER is left out
(record order from taxid_order.tsv). Writes OUT/<asm>.xo.tsv with d_Hx_xo, nHx_10kb_xo,
d_HD_xo, d_FLANK_xo, n_HD_20kb_xo, n_FLANK_50kb_xo per locus.
Usage: recompute_xo.py ASM OUTDIR
"""
import bisect
import importlib.util
import os
import sys
import tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "build_genome.py")).read().replace("\nmain()\n", "\n")
ns = {"__file__": os.path.join(HERE, "build_genome.py"), "__name__": "bg"}
exec(compile(src, "build_genome.py", "exec"), ns)
sg, paf_hits, nearest, OWN_TAX, THREADS, MAT, BIG = (ns[k] for k in ("sg", "paf_hits", "nearest", "OWN_TAX", "THREADS", "MAT", "BIG"))
ORDER_OF = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "taxid_order.tsv")))
asm, out = sys.argv[1], sys.argv[2]
my_order = ORDER_OF[OWN_TAX[asm[:15]]]
bad = {t + "_" for t, o in ORDER_OF.items() if o == my_order}


def excluded(q):
    return any(q.startswith(b) for b in bad)


gz = os.path.join(ns["LIB"], asm + ".fa.gz")
scratch = os.environ.get("SCRATCH", "/bigdata/stajichlab/jstajich/prml_work/tmp")
tmp = tempfile.mkdtemp(prefix="xo_" + asm + "_", dir=scratch)
fa = os.path.join(tmp, asm + ".fa")
seqs = sg.read_fasta(gz)
with open(fa, "w") as fo:
    for k, v in seqs.items():
        fo.write(f">{k}\n{v}\n")
homs = [h for h in sg.homology_hits(fa, tmp) if not excluded(h[2])]
H = defaultdict(list)
for ctg, pos, q, ev in homs:
    if ev <= 1:
        H[ctg].append(pos)
for k in H:
    H[k].sort()
mh = defaultdict(lambda: defaultdict(list))
for q, qcov, t, ts, te, sc in paf_hits(gz, MAT, THREADS):
    if qcov >= 0.5 and sc >= 150 and te - ts < 60000 and not excluded(q):
        mh[q.split("|")[1]][t].append((ts + te) // 2)
for c in mh:
    for t in mh[c]:
        mh[c][t].sort()
rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(out, asm + ".loci.tsv"))]
hdr = rows[0]
ix = {h: i for i, h in enumerate(hdr)}
with open(os.path.join(out, asm + ".xo.tsv"), "w") as fo:
    fo.write("asm\tcontig\tstart\tend\td_Hx_xo\tnHx_10kb_xo\td_HD_xo\td_FLANK_xo\tn_HD_20kb_xo\tn_FLANK_50kb_xo\n")
    for r in rows[1:]:
        c, a, b = r[ix["contig"]], int(r[ix["start"]]), int(r[ix["end"]])

        def cnt(arr, W):
            return bisect.bisect_right(arr, b + W) - bisect.bisect_left(arr, a - W)
        fo.write("\t".join(map(str, [asm, c, a, b, nearest(H, c, a, b), cnt(H.get(c, []), 10000),
                                      nearest(mh["HD"], c, a, b), nearest(mh["FLANK"], c, a, b),
                                      cnt(mh["HD"].get(c, []), 20000), cnt(mh["FLANK"].get(c, []), 50000)])) + "\n")
os.remove(fa)
print(asm, "xo done")
