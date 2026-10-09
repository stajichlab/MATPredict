#!/usr/bin/env python3
"""B-locus clustering assessment. Reuses the array-study outputs (inventory table, scan_out cand.tsv
= strict-CAAX ORFs 'T' and tblastn precursor homology 'Hx' E<=1 from scan_genome.py / array_scan.py).
Usage: b_cluster.py SCAN_OUT_DIR HX_HELDOUT_DIR(optional)"""
import sys, os, glob, re, gzip
import numpy as np, pandas as pd
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); ST = HERE + "/../2026-10-06_agaricomycetes_pr_arrays"
SCAN = sys.argv[1]; REPO = HERE + "/../.."
rng = np.random.default_rng(1)
t = pd.read_csv(ST + "/inventory/ste3_loci_table.tsv.gz", sep="\t", low_memory=False)
q = t[t.qpass].copy()
q["array_n"] = q.array_id.str.split("|").str[1].astype(int)
TOL = 300

def load_cand(asm):
    T, H = defaultdict(list), defaultdict(list)
    for l in open(f"{SCAN}/{asm}.cand.tsv"):
        k, c, p = l.rstrip("\n").split("\t")
        if k == "T": T[c].append(int(p))
        elif k == "Hx": H[c].append(int(p))
    T = {c: np.array(sorted(set(v))) for c, v in T.items()}
    Hc = {}
    for c, v in H.items():
        v = sorted(v); cl = [v[0]]
        for x in v[1:]:
            if x - cl[-1] > TOL: cl.append(x)
        Hc[c] = np.array(cl)
    return T, Hc

def cnt(arr, a, b):
    if arr is None or len(arr) == 0: return 0
    return int(np.searchsorted(arr, b, "right") - np.searchsorted(arr, a, "left"))

def near(arr, p):  # T positions with an H cluster within TOL
    return 0

def feats(T, H, c, a, b, W):
    Tc, Hc = T.get(c, np.array([])), H.get(c, np.array([]))
    ta = Tc[(Tc >= a - W) & (Tc <= b + W)] if len(Tc) else Tc
    ha = Hc[(Hc >= a - W) & (Hc <= b + W)] if len(Hc) else Hc
    nT, nH = len(ta), len(ha)
    nTH = int(sum(np.any(np.abs(ha - x) <= TOL) for x in ta)) if nT and nH else 0
    nHonly = int(sum(not np.any(np.abs(ta - x) <= TOL) for x in ha)) if nH else 0
    return nT, nH, nTH, nT + nHonly

rows = []; loc_rows = []
for asm, g in q.groupby("genome"):
    T, H = load_cand(asm)
    # background T: not within 10 kb of any locus
    bgT = {}
    for c, arr in T.items():
        gl = g[g.contig == c]
        keep = np.ones(len(arr), bool)
        for s, e in zip(gl.start, gl.end): keep &= ~((arr >= s - 10000) & (arr <= e + 10000))
        bgT[c] = arr[keep]
    clen = g.drop_duplicates("contig").set_index("contig").contig_len.to_dict()
    for _, r in g.iterrows():
        nT, nH, nTH, nU = feats(T, H, r.contig, r.start, r.end, 5000)
        L = r.end - r.start
        # contig-uniform null: same locus length, random position on its contig, T from background only
        n_rep = 40; cl = clen[r.contig]
        ok2 = 0
        if cl > L + 10000:
            pos = rng.integers(5000, max(5001, cl - L - 5000), n_rep)
            Tb = bgT.get(r.contig, np.array([]))
            if len(Tb):
                lo = np.searchsorted(Tb, pos - 5000, "left"); hi = np.searchsorted(Tb, pos + L + 5000, "right")
                ok2 = float(np.mean((hi - lo) >= 2))
        loc_rows.append(dict(locus_id=r.locus_id, genome=asm, nT5=nT, nH5=nH, nTH5=nTH, nU5=nU,
                             casA=nT >= 2, casB=(nU >= 2 and nTH >= 1), casC=nTH >= 2, null_casA=ok2))
loc = pd.DataFrame(loc_rows)
q = q.merge(loc, on=["locus_id", "genome"])
q.to_csv(HERE + "/b_loci_cassette.tsv.gz", sep="\t", index=False)

# arrays
ar = []
for aid, g in q.groupby("array_id"):
    asm = g.genome.iloc[0]
    ar.append(dict(array_id=aid, genome=asm, order=g.order.iloc[0], family=g.family.iloc[0], species=g.species.iloc[0], contig=g.contig.iloc[0],
                   start=g.start.min(), end=g.end.max(), n_loci=len(g), span=g.end.max() - g.start.min(),
                   caax_flag=g.caax_flag.any(), in_call=g.pipeline_call.any(), withheld=g.withheld_cluster.any(),
                   max_cas_nT5=g.nT5.max(), any_casA=g.casA.any(), any_casB=g.casB.any(), any_casC=g.casC.any(),
                   n_casA=int(g.casA.sum()), exp_casA=g.null_casA.sum()))
A = pd.DataFrame(ar)
# precursors within array extent +-10kb (per array, distinct)
cache = {}
def arr_prec(r):
    if r.genome not in cache: cache.clear(); cache[r.genome] = load_cand(r.genome)
    T, H = cache[r.genome]
    return feats(T, H, r.contig, r.start, r.end, 10000)
A[["aT10", "aH10", "aTH10", "aU10"]] = A.apply(lambda r: pd.Series(arr_prec(r)), axis=1)
ar0 = pd.read_csv(ST + "/arrays.tsv", sep="\t")
A["arr_n"] = A.array_id.str.split("|").str[1].astype(int)
A = A.merge(ar0[["genome", "array", "d_STE20", "has_STE20"]].rename(columns={"array": "arr_n"}), on=["genome", "arr_n"], how="left")
A["ste20_250"] = A.d_STE20 <= 250000
A.to_csv(HERE + "/b_arrays_features.tsv.gz", sep="\t", index=False)
print("loci", len(q), "arrays", len(A), "genomes", A.genome.nunique())
