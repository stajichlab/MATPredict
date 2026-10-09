#!/usr/bin/env python3
"""Independent re-derivation and controls for the Agaricomycete receptor-array study (full scan, 1,853 genomes).

Inputs (this directory): agari_genomes.tsv, genome_table.tsv, scan_out/ (array_scan.py output, unpacked from
scan_out_full.tar.xz), clen/ (contig lengths, contig_lengths.sh), pr_calls.tsv.gz, hd_loci.tsv.gz.
Nothing here uses a label, a PR-only call or CAAX status as ground truth.  Outputs: v_*.tsv and v_summary.txt.

Sections
 A  sample (per-order genome / species counts at each filter)
 B  array definition re-derived independently; sensitivity to the clustering gap; clustering vs a within-contig null
 C  CAAX flag chance model at array level (within-contig analytic null; genome null of the pipeline for comparison);
    sensitivity to the flag definition (window 5/10/20/50 kb; >=2 CAAX ORFs; + precursor homology)
 D  STE20 (and MIPBF, HD) proximity of flagged vs unflagged arrays with position-matched nulls, within-genome
    comparison, permutation test, leave-one-order-out
 E  scenarios S0-S5 re-derived, with array-level chance for S2
 F  withheld loci inside arrays
Intervals: percentile bootstrap resampling species (cluster unit); one-genome-per-species set as a check.
"""
import glob, os, sys
from collections import defaultdict
import numpy as np
import pandas as pd

RNG = np.random.default_rng(7)
NB = 2000
OUT = open("v_summary.txt", "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); OUT.write(s + "\n"); OUT.flush()
def T(df, name, **kw):
    df.to_csv(name, sep="\t", index=False)
    P(f"\n== {name}"); P(df.to_string(index=False, **kw))

ORDERS = ["Agaricales", "Boletales", "Polyporales", "Cantharellales", "Russulales", "Hymenochaetales", "Auriculariales"]
def grp(o): return o if o in ORDERS else "other orders"

# ------------------------------------------------------------------ load
gen = pd.read_csv("agari_genomes.tsv", sep="\t")
gen["sp"] = gen.species.fillna("").str.split().str[:2].str.join(" ")
gen["ord"] = gen.order.map(grp)
scanned = {os.path.basename(f)[:-11] for f in glob.glob("scan_out/*.chance.tsv")}
gen["scanned"] = gen.genome.isin(scanned)
gen["ok"] = gen.scanned & gen.qpass & (gen.contigs <= 5000)
# one genome per species: best BUSCO then N50
q = gen[gen.ok].sort_values(["busco_complete_pct", "n50_bp"], ascending=False)
one = set(q.drop_duplicates("sp").genome)
gen["one_per_sp"] = gen.genome.isin(one)
meta = gen.set_index("genome")

def rd(suffix):
    out = []
    for f in glob.glob(f"scan_out/*.{suffix}.tsv"):
        d = pd.read_csv(f, sep="\t"); d["genome"] = os.path.basename(f)[:-len(suffix) - 5]
        if len(d): out.append(d)
    return pd.concat(out, ignore_index=True)
loci = rd("loci").drop(columns=["asm"], errors="ignore")
reg = rd("region"); cand = rd("cand"); chance = rd("chance")
clen = {}
for f in glob.glob("clen/*.tsv"):
    g = os.path.basename(f)[:-4]
    d = pd.read_csv(f, sep="\t", header=None, names=["contig", "len"], dtype={"contig": str})
    clen[g] = dict(zip(d.contig, d.len))
P(f"genomes with contig lengths: {len(clen)}")
loci = loci[loci.genome.isin(gen[gen.ok].genome)].copy()
loci["order"] = loci.genome.map(meta["ord"]); loci["sp"] = loci.genome.map(meta.sp)

# ================================================================== A  sample
rows = []
for o in ORDERS + ["other orders", "ALL"]:
    s = gen if o == "ALL" else gen[gen.ord == o]
    r = s[s.ok]
    rows.append(dict(order=o, genomes_all=len(s), species_all=s.sp.nunique(), scanned=int(s.scanned.sum()), quality_pass=int((s.scanned & s.qpass).sum()),
                     plus_contigs_le5000=len(r), species_ok=r.sp.nunique(), one_per_species=int(r.one_per_sp.sum()),
                     busco_median=r.busco_complete_pct.median(), n50_median_kb=round(r.n50_bp.median() / 1000), contigs_median=r.contigs.median()))
SA = pd.DataFrame(rows); T(SA, "v_sample.tsv")
P("note: qpass = BUSCO complete >= 70 and N50 >= 20 kb; contigs <= 5000 removes",
  int((gen.scanned & gen.qpass & (gen.contigs > 5000)).sum()), "genomes")

# ================================================================== B  arrays
def make_arrays(df, gap):
    df = df.sort_values(["genome", "contig", "start"]).copy()
    key = df.genome + "|" + df.contig
    cm = df.groupby(key, sort=False).end.cummax()
    prev = cm.groupby(key, sort=False).shift()
    first = key != key.shift()
    new = first | ((df.start - prev) > gap)
    df["array"] = new.cumsum().values
    return df
Lg = {}
for gap in (10_000, 25_000, 50_000, 100_000, 200_000):
    Lg[gap] = make_arrays(loci, gap)
old = pd.read_csv("loci_all.tsv.gz", sep="\t")
cmp = Lg[50_000].merge(old[["genome", "contig", "start", "end", "array"]], on=["genome", "contig", "start", "end"], suffixes=("", "_old"))
same = cmp.groupby("array").array_old.nunique().max() == 1 and cmp.groupby("array_old").array.nunique().max() == 1
P(f"\nB. my arrays (gap 50 kb) reproduce analyze.py partition on {len(cmp)} loci: {same}")
L50 = Lg[50_000]
# loci whose intervals overlap an adjacent locus on the opposite strand
ls = L50.sort_values(["genome", "contig", "start"]); prevend = ls.groupby(["genome", "contig"]).end.shift()
P(f"loci overlapping the previous locus (opposite strand): {(ls.start < prevend).sum()} of {len(ls)} ({100*(ls.start < prevend).mean():.1f}%)")

def species_boot(df, f, key="sp", nb=NB):
    """f(subdf)->float; resample species with replacement"""
    groups = {k: v for k, v in df.groupby(key)}
    keys = list(groups); pt = f(df)
    if len(keys) < 2: return pt, np.nan, np.nan
    res = []
    for _ in range(nb):
        pick = RNG.integers(0, len(keys), len(keys))
        res.append(f(pd.concat([groups[keys[i]] for i in pick])))
    return pt, *np.nanpercentile(res, [2.5, 97.5])
def fmt(t, d=3):
    return f"{t[0]:.{d}f} ({t[1]:.{d}f}-{t[2]:.{d}f})"

# precompute genome-level aggregates for fast bootstrap: dict sp -> sums
def boot_sums(gdf, cols, nb=NB):
    """gdf one row per genome, with sp; returns function of ratio of column sums with species bootstrap"""
    sp = gdf.groupby("sp")[cols].sum()
    v = sp.values; n = len(v)
    idx = RNG.integers(0, n, size=(nb, n))
    return v, idx
def ratio_ci(gdf, num, den):
    v, idx = boot_sums(gdf, [num, den])
    pt = v[:, 0].sum() / v[:, 1].sum()
    r = v[idx, 0].sum(1) / np.maximum(v[idx, 1].sum(1), 1e-9)
    return pt, *np.percentile(r, [2.5, 97.5])

rows = []
for gap, A in Lg.items():
    arr = A.groupby("array").agg(genome=("genome", "first"), n=("start", "size"), contig=("contig", "first"), start=("start", "min"), end=("end", "max"),
                                 nflag=("nT_10kb", lambda s: int((s > 0).sum())))
    arr["order"] = arr.genome.map(meta["ord"])
    ng = arr.genome.nunique()
    rows.append(dict(gap_kb=gap // 1000, genomes=ng, loci=len(A), arrays=len(arr), arrays_per_genome=round(len(arr) / ng, 2),
                     singleton_arrays_pct=round(100 * (arr.n == 1).mean(), 1), loci_in_arrays_ge2_pct=round(100 * A.groupby("array").size().pipe(lambda s: s[s >= 2].sum()) / len(A), 1),
                     largest_median=arr.groupby("genome").n.max().median(), flagged_arrays=int((arr.nflag > 0).sum()),
                     flagged_pct_of_arrays=round(100 * (arr.nflag > 0).mean(), 1)))
T(pd.DataFrame(rows), "v_gap_sensitivity.tsv")

# clustering vs within-contig null: place the same number of loci uniformly on each contig (locus length kept), same gap
def cluster_null(A, gap, nrep=20):
    A = A.copy()
    A["len"] = A.end - A.start
    obs = (A.groupby("array").size().pipe(lambda s: s[s >= 2].sum())) / len(A)
    res = []
    g = A.groupby(["genome", "contig"])
    ctg = [(k, v.len.values, v.contig_len.iloc[0]) for k, v in g]
    for _ in range(nrep):
        inarr = 0
        for k, ln, L in ctg:
            if len(ln) == 1:
                continue
            st = np.sort(RNG.integers(0, max(1, L - ln.max()), len(ln)))
            en = st + ln[np.argsort(RNG.random(len(ln)))]
            # arrays by gap on this contig
            run = 1; sizes = []
            cur = en[0]
            for i in range(1, len(st)):
                if st[i] - cur > gap:
                    sizes.append(run); run = 1; cur = en[i]
                else:
                    run += 1; cur = max(cur, en[i])
            sizes.append(run)
            inarr += sum(s for s in sizes if s >= 2)
        res.append(inarr / len(A))
    return obs, np.mean(res), np.min(res), np.max(res)
rows = []
for gap in (10_000, 50_000, 200_000):
    for o in ["ALL"] + ORDERS:
        A = Lg[gap] if o == "ALL" else Lg[gap][Lg[gap].order == o]
        ob, mu, lo, hi = cluster_null(A, gap, nrep=10)
        rows.append(dict(gap_kb=gap // 1000, order=o, loci=len(A), obs_frac_loci_in_arrays_ge2=round(ob, 3), null_mean=round(mu, 3), null_min=round(lo, 3), null_max=round(hi, 3)))
T(pd.DataFrame(rows), "v_clustering_vs_null.tsv")

# ================================================================== C  CAAX flag chance (array level)
# T ORF positions per genome/contig, dropping those within 10 kb of any STE3 locus (the signal, not the background)
Tpos = defaultdict(dict)
for g, d in cand[cand.kind == "T"].groupby("genome"):
    if g not in set(loci.genome.unique()):
        continue
    for c, e in d.groupby("contig"):
        Tpos[g][c] = np.sort(e.pos.values)
locsB = {k: v for k, v in L50.groupby(["genome", "contig"])}
def bg_T(g, c, W=10_000):
    t = Tpos[g].get(c)
    if t is None: return np.array([])
    l = locsB.get((g, c))
    if l is None: return t
    keep = np.ones(len(t), bool)
    for a, b in zip(l.start.values, l.end.values):
        keep &= ~((t >= a - W) & (t <= b + W))
    return t[keep]
def union_len(iv, lo, hi):
    """total length of union of intervals iv (n,2) clipped to [lo,hi]"""
    if len(iv) == 0: return 0.0
    iv = iv[np.argsort(iv[:, 0])]
    tot = 0.0; cs, ce = None, None
    for a, b in iv:
        a = max(a, lo); b = min(b, hi)
        if b <= a: continue
        if cs is None: cs, ce = a, b
        elif a <= ce: ce = max(ce, b)
        else:
            tot += ce - cs; cs, ce = a, b
    if cs is not None: tot += ce - cs
    return tot
def p_window(t, L, s, W, k=1):
    """P(a window of span s placed uniformly on a contig of length L has >=k points of t within W of it)"""
    hi = L - s
    if hi <= 0: hi = 1.0
    if len(t) < k: return 0.0
    if k == 1: iv = np.stack([t - s - W, t + W], 1)
    else: iv = np.stack([t[k - 1:] - s - W, t[:len(t) - k + 1] + W], 1)
    iv = iv[iv[:, 1] > iv[:, 0]]
    return union_len(iv, 0, hi) / hi

def build_arrays(A, W=10_000):
    arr = A.groupby("array").agg(genome=("genome", "first"), contig=("contig", "first"), start=("start", "min"), end=("end", "max"), n=("start", "size"),
                                 contig_len=("contig_len", "first"), dT=("d_T", "min"), nT=("nT_10kb", "max"), nHx=("nHx_10kb", "max"), dHx=("d_Hx", "min")).reset_index()
    arr["span"] = arr.end - arr.start
    arr["order"] = arr.genome.map(meta["ord"]); arr["sp"] = arr.genome.map(meta.sp)
    arr["one"] = arr.genome.map(meta.one_per_sp)
    return arr
arr = build_arrays(L50)
for W in (5_000, 10_000, 20_000, 50_000):
    arr[f"flag{W//1000}"] = arr.dT <= W
arr["flag"] = arr.flag10
arr["hx"] = arr.dHx <= 10_000
# chance per array under the within-contig null, several definitions
def chance_cols(arr):
    out = {k: [] for k in ("p5", "p10", "p20", "p50", "p10_k2")}
    for r in arr.itertuples():
        L = r.contig_len
        for W in (5_000, 10_000, 20_000, 50_000):
            t = bg_T(r.genome, r.contig, W)
            out[f"p{W//1000}"].append(p_window(t, L, r.span, W))
            if W == 10_000:
                out["p10_k2"].append(p_window(t, L, r.span, W, k=2))
    for k, v in out.items(): arr[k] = v
chance_cols(arr)
# observed >=2 CAAX ORFs within 10 kb of the array: from loci nT_10kb (count within 10 kb of a locus, so use max)
arr["flag10_k2"] = arr.nT >= 2
arr.to_csv("v_arrays.tsv.gz", sep="\t", index=False)

def order_table(df, defs):
    rows = []
    for o in ORDERS + ["other orders", "ALL"]:
        s = df if o == "ALL" else df[df.order == o]
        for name, (obs, p) in defs.items():
            sub = s.assign(O=s[obs].astype(int), E=s[p])
            gs = sub.groupby(["sp"])[["O", "E"]].sum().reset_index()
            pt, lo, hi = ratio_ci(sub.assign(sp=sub.sp), "E", "O") if sub.O.sum() > 0 else (np.nan,) * 3
            rows.append(dict(order=o, flag_def=name, arrays=len(sub), observed_flagged=int(sub.O.sum()), expected_by_chance=round(sub.E.sum(), 1),
                             chance_share=round(pt, 3), lo=round(lo, 3), hi=round(hi, 3)))
    return pd.DataFrame(rows)
defs = {"T within 10kb (pipeline)": ("flag10", "p10")}
C1 = order_table(arr, defs); T(C1, "v_array_chance_by_order.tsv")
defs = {"T within 5kb": ("flag5", "p5"), "T within 10kb": ("flag10", "p10"), "T within 20kb": ("flag20", "p20"), "T within 50kb": ("flag50", "p50"),
        ">=2 T within 10kb": ("flag10_k2", "p10_k2")}
C2 = order_table(arr, defs)
T(C2[C2.order.isin(["ALL", "Agaricales", "Polyporales", "Boletales", "Cantharellales", "Russulales", "Hymenochaetales"])], "v_flag_definition_sensitivity.tsv")
# one genome per species check
C3 = order_table(arr[arr.one], {"T within 10kb, one genome per species": ("flag10", "p10")})
T(C3, "v_array_chance_one_per_species.tsv")
# comparison of the three chance numbers at locus level for ALL
Gp = chance.set_index("genome").p_random
loci["p_pipe"] = loci.genome.map(Gp); loci["flag"] = loci.nT_10kb > 0
P(f"\nlocus level, ALL: flagged {loci.flag.sum()}, expected by pipeline p_random {loci.p_pipe.sum():.1f} (share {loci.p_pipe.sum()/loci.flag.sum():.3f})")
P(f"array level (within-contig null) flagged {arr.flag.sum()}, expected {arr.p10.sum():.1f} (share {arr.p10.sum()/arr.flag.sum():.3f})")
# flagged by array size
rows = []
for lab, s in (("1", arr[arr.n == 1]), ("2", arr[arr.n == 2]), ("3", arr[arr.n == 3]), ("4+", arr[arr.n >= 4])):
    rows.append(dict(array_size=lab, arrays=len(s), flagged=int(s.flag.sum()), flagged_pct=round(100 * s.flag.mean(), 1), expected_chance=round(s.p10.sum(), 1),
                     chance_pct_per_array=round(100 * s.p10.mean(), 1)))
T(pd.DataFrame(rows), "v_flag_by_array_size.tsv")

# ================================================================== D  region-gene proximity of flagged vs unflagged arrays
def top_hits(cls):
    d = reg[reg["class"] == cls]
    d = d[d.genome.isin(gen[gen.ok].genome)]
    d = d.sort_values(["score", "qcov"], ascending=False).groupby("genome").head(1)
    return d.set_index("genome")
def sumlen(g):
    c = clen.get(g)
    return np.sort(np.array(list(c.values()))) if c else None
sorted_len = {g: sumlen(g) for g in arr.genome.unique()}
def near_measure(a0, a1, x, y, L, D):
    """measure of window starts s in [0,L-span] with window within D of [x,y]"""
    span = a1 - a0; hi = L - span
    lo_s = max(0, x - D - span); hi_s = min(hi, y + D)
    return max(0.0, hi_s - lo_s), max(hi, 1.0)
def genome_denominator(g, span):
    ls = sorted_len.get(g)
    if ls is None: return np.nan
    return float(np.clip(ls - span, 0, None).sum())
def region_cols(arr, cls, Ds=(100_000, 250_000, 500_000)):
    th = top_hits(cls)
    res = {f"{cls}_obs{D//1000}": [] for D in Ds}
    res.update({f"{cls}_pc{D//1000}": [] for D in Ds}); res.update({f"{cls}_pg{D//1000}": [] for D in Ds})
    res[f"{cls}_samecontig"] = []; res[f"{cls}_has"] = []; res[f"{cls}_dist"] = []
    for r in arr.itertuples():
        has = r.genome in th.index
        res[f"{cls}_has"].append(has)
        if has:
            t = th.loc[r.genome]
            same = t.contig == r.contig
            d = max(0, max(r.start, t.start) - min(r.end, t.end)) if same else np.inf
        else:
            same, d = False, np.inf
        res[f"{cls}_samecontig"].append(same); res[f"{cls}_dist"].append(d)
        for D in Ds:
            res[f"{cls}_obs{D//1000}"].append(d <= D)
            if has and same:
                m, den = near_measure(r.start, r.end, t.start, t.end, r.contig_len, D)
                res[f"{cls}_pc{D//1000}"].append(m / den)
            else:
                res[f"{cls}_pc{D//1000}"].append(0.0)
            if has:
                # genome-position null: a random window anywhere in the genome; only the hit's contig can be near
                tl = clen[r.genome].get(t.contig, r.contig_len) if r.genome in clen else r.contig_len
                m, _ = near_measure(r.start, r.end, t.start, t.end, tl, D)
                res[f"{cls}_pg{D//1000}"].append(m / genome_denominator(r.genome, r.span))
            else:
                res[f"{cls}_pg{D//1000}"].append(0.0)
    for k, v in res.items(): arr[k] = v
for cls in ("STE20", "MIPBF"):
    region_cols(arr, cls)
# HD positions (pipeline called HD / bLocus loci) -- nearest, not top hit
hd = pd.read_csv("hd_loci.tsv.gz", sep="\t"); hd = hd[(hd.status == "called") & hd.genome.isin(arr.genome.unique())]
hdby = {k: v for k, v in hd.groupby("genome")}
arr["HD_has"] = arr.genome.isin(hdby.keys())
dHD = []
for r in arr.itertuples():
    h = hdby.get(r.genome)
    if h is None: dHD.append(np.inf); continue
    s = h[h.contig == r.contig]
    dHD.append(min([max(0, max(r.start, x.start) - min(r.end, x.end)) for x in s.itertuples()], default=np.inf))
arr["HD_dist"] = dHD
arr.to_csv("v_arrays.tsv.gz", sep="\t", index=False)

def enrich(cls, D, df):
    rows = []
    for o in ["ALL"] + ORDERS + ["other orders"]:
        s = df if o == "ALL" else df[df.order == o]
        s = s[s[f"{cls}_has"]]
        for lab, g in (("flagged", s[s.flag]), ("unflagged", s[~s.flag])):
            if len(g) == 0: continue
            O = g[f"{cls}_obs{D//1000}"].sum(); Ec = g[f"{cls}_pc{D//1000}"].sum(); Eg = g[f"{cls}_pg{D//1000}"].sum()
            rows.append(dict(order=o, group=lab, arrays=len(g), contig_len_median_Mb=round(g.contig_len.median() / 1e6, 2), same_contig=int(g[f"{cls}_samecontig"].sum()),
                             near_obs=int(O), near_pct=round(100 * O / len(g), 1), exp_genome_null=round(Eg, 1), exp_pct_genome=round(100 * Eg / len(g), 1),
                             obs_over_exp_genome=round(O / Eg, 1) if Eg > 0 else np.nan,
                             exp_same_contig_null=round(Ec, 1)))
    return pd.DataFrame(rows)
for cls in ("STE20", "MIPBF"):
    for D in (250_000,):
        T(enrich(cls, D, arr), f"v_{cls}_enrichment_{D//1000}kb.tsv")
# STE20 near-fraction by distance threshold (ALL)
rows = []
for cls in ("STE20", "MIPBF"):
    for D in (100_000, 250_000, 500_000):
        s = arr[arr[f"{cls}_has"]]
        for lab, g in (("flagged", s[s.flag]), ("unflagged", s[~s.flag])):
            rows.append(dict(gene=cls, D_kb=D // 1000, group=lab, arrays=len(g), near=int(g[f"{cls}_obs{D//1000}"].sum()), expected_genome_null=round(g[f"{cls}_pg{D//1000}"].sum(), 1),
                             expected_same_contig_null=round(g[f"{cls}_pc{D//1000}"].sum(), 1)))
T(pd.DataFrame(rows), "v_region_by_threshold.tsv")

# --- conditional on contig: arrays on the same contig as the STE20 top hit
s = arr[arr.STE20_samecontig]
rows = []
for D in (100_000, 250_000, 500_000):
    for lab, g in (("flagged", s[s.flag]), ("unflagged", s[~s.flag])):
        rows.append(dict(D_kb=D // 1000, group=lab, arrays_on_STE20_contig=len(g), near=int(g[f"STE20_obs{D//1000}"].sum()),
                         expected_if_position_random_on_contig=round(g[f"STE20_pc{D//1000}"].sum(), 1)))
T(pd.DataFrame(rows), "v_STE20_conditional_on_contig.tsv")

# --- within-genome comparison: CMH over genomes with >=1 flagged and >=1 unflagged array and a STE20 hit; and permutation of flag labels
from statsmodels.stats.contingency_tables import StratifiedTable
def within_genome(df, D=250_000, nperm=5000, label=""):
    d = df[df.STE20_has]
    tabs = []; keep = []
    for g, x in d.groupby("genome"):
        if x.flag.any() and (~x.flag).any():
            a = int(x[x.flag][f"STE20_obs{D//1000}"].sum()); b = int(x.flag.sum()) - a
            c = int(x[~x.flag][f"STE20_obs{D//1000}"].sum()); dd = int((~x.flag).sum()) - c
            tabs.append(np.array([[a, b], [c, dd]])); keep.append(g)
    st = StratifiedTable(tabs)
    # permutation: shuffle the flag labels among arrays inside each genome
    obs = sum(t[0, 0] for t in tabs)
    by = {g: (x[f"STE20_obs{D//1000}"].values.astype(int), int(x.flag.sum())) for g, x in d[d.genome.isin(keep)].groupby("genome")}
    sims = np.zeros(nperm)
    # vectorised: per genome sample k of n without replacement nperm times
    for g, (near, k) in by.items():
        n = len(near)
        idx = np.argsort(RNG.random((nperm, n)), axis=1)[:, :k]
        sims += near[idx].sum(1)
    pperm = (np.sum(sims >= obs) + 1) / (nperm + 1)
    # species-cluster bootstrap of the pooled MH odds ratio
    gsp = {g: meta.loc[g, "sp"] for g in keep}
    sp_keys = sorted(set(gsp.values())); bysp = defaultdict(list)
    for g, t in zip(keep, tabs): bysp[gsp[g]].append(t)
    ors = []
    for _ in range(1000):
        pick = RNG.integers(0, len(sp_keys), len(sp_keys))
        tt = [t for i in pick for t in bysp[sp_keys[i]]]
        try:
            ors.append(StratifiedTable(tt).oddsratio_pooled)
        except Exception:
            pass
    lo, hi = np.nanpercentile([o for o in ors if np.isfinite(o)], [2.5, 97.5])
    return dict(scope=label, D_kb=D // 1000, genomes_with_both=len(keep), flagged_near=int(obs), flagged_total=int(sum(t[0].sum() for t in tabs)),
                unflagged_near=int(sum(t[1, 0] for t in tabs)), unflagged_total=int(sum(t[1].sum() for t in tabs)),
                MH_odds_ratio=round(st.oddsratio_pooled, 2), MH_species_boot_lo=round(lo, 2), MH_species_boot_hi=round(hi, 2),
                CMH_p=float(st.test_null_odds().pvalue), perm_p=round(pperm, 4), mean_perm_flagged_near=round(sims.mean(), 1))
rows = [within_genome(arr, 250_000, label="ALL qpass genomes"), within_genome(arr, 100_000, label="ALL qpass genomes"), within_genome(arr, 500_000, label="ALL qpass genomes"),
        within_genome(arr[arr.one], 250_000, label="one genome per species")]
for o in ORDERS:
    try: rows.append(within_genome(arr[arr.order == o], 250_000, nperm=2000, label=o))
    except Exception as e: rows.append(dict(scope=o, note=str(e)[:60]))
# leave-one-order-out
for o in ORDERS:
    rows.append(within_genome(arr[arr.order != o], 250_000, nperm=2000, label="all except " + o))
WG = pd.DataFrame(rows); T(WG, "v_STE20_within_genome.tsv")
# same for MIPBF (negative-control gene)
def within_genome_cls(df, cls, D=250_000):
    d = df[df[f"{cls}_has"]]
    tabs = []
    for g, x in d.groupby("genome"):
        if x.flag.any() and (~x.flag).any():
            a = int(x[x.flag][f"{cls}_obs{D//1000}"].sum()); b = int(x.flag.sum()) - a
            c = int(x[~x.flag][f"{cls}_obs{D//1000}"].sum()); dd = int((~x.flag).sum()) - c
            tabs.append(np.array([[a, b], [c, dd]]))
    st = StratifiedTable(tabs)
    return dict(gene=cls, genomes_with_both=len(tabs), flagged_near=int(sum(t[0, 0] for t in tabs)), flagged_total=int(sum(t[0].sum() for t in tabs)),
                unflagged_near=int(sum(t[1, 0] for t in tabs)), unflagged_total=int(sum(t[1].sum() for t in tabs)),
                MH_odds_ratio=round(st.oddsratio_pooled, 2), CMH_p=float(st.test_null_odds().pvalue))
T(pd.DataFrame([within_genome_cls(arr, "MIPBF")]), "v_MIPBF_within_genome.tsv")
# contig length confounding: flagged vs unflagged contig lengths, array span
P("\nflagged vs unflagged arrays: median contig length (Mb), median array size, median span (kb)")
for lab, g in (("flagged", arr[arr.flag]), ("unflagged", arr[~arr.flag])):
    P(lab, len(g), round(g.contig_len.median() / 1e6, 2), g.n.median(), round(g.span.median() / 1e3, 1))
# species- and order-bootstrapped difference in near-fractions (flagged minus unflagged, genome-null-corrected excess)
def excess_stat(df, D=250_000, cls="STE20"):
    s = df[df[f"{cls}_has"]]
    f = s[s.flag]; u = s[~s.flag]
    return (f[f"{cls}_obs{D//1000}"].sum() - f[f"{cls}_pg{D//1000}"].sum()) / max(len(f), 1) - (u[f"{cls}_obs{D//1000}"].sum() - u[f"{cls}_pg{D//1000}"].sum()) / max(len(u), 1)
t = species_boot(arr[arr.STE20_has], excess_stat, nb=1000)
P(f"\nexcess near-STE20 fraction beyond genome-position null, flagged minus unflagged, species bootstrap: {fmt(t)}")
# order-level: per-order flagged near fraction vs unflagged
rows = []
for o in ORDERS + ["other orders"]:
    s = arr[(arr.order == o) & arr.STE20_has]
    if len(s) == 0: continue
    rows.append(dict(order=o, species=s.sp.nunique(), flagged=int(s.flag.sum()), flagged_near_pct=round(100 * s[s.flag].STE20_obs250.mean(), 1) if s.flag.any() else np.nan,
                     unflagged=int((~s.flag).sum()), unflagged_near_pct=round(100 * s[~s.flag].STE20_obs250.mean(), 1)))
T(pd.DataFrame(rows), "v_STE20_by_order.tsv")
# per-genome concentration: how many genomes/species hold the flagged-near arrays
fn = arr[arr.flag & arr.STE20_obs250]
P(f"\nflagged arrays near STE20 (250 kb): {len(fn)} in {fn.genome.nunique()} genomes, {fn.sp.nunique()} species, {fn.order.nunique()} order groups;",
  f"flagged arrays near STE20 per genome: {fn.groupby('genome').size().describe()[['mean','50%','max']].round(2).to_dict()}")
# HD: same contig
d = arr[arr.HD_has]
exp_same = []
for r in d.itertuples():
    c = clen.get(r.genome, {}); hdc = hdby[r.genome]
    tot = sum(c.values()) or np.nan
    # P(a random point shares a contig with at least one HD call)  -- contigs holding HD calls, length weighted
    exp_same.append(sum(c.get(x, 0) for x in set(hdc.contig)) / tot)
d = d.assign(exp_same=exp_same, same=np.isfinite(d.HD_dist))
P(f"\nHD calls: arrays in genomes with an HD call {len(d)}; on a contig that also holds an HD call {int(d.same.sum())}; expected if arrays placed at random by contig length {d.exp_same.sum():.1f}")
for lab, g in (("flagged", d[d.flag]), ("unflagged", d[~d.flag])):
    P(f"  {lab}: {len(g)} arrays, same contig as HD {int(g.same.sum())}, expected {g.exp_same.sum():.1f}")

# ================================================================== E  scenarios and calls
pr = pd.read_csv("pr_calls.tsv.gz", sep="\t")
ok_g = set(gen[gen.ok].genome)
called = pr[(pr.status == "called") & (pr.family == "PR") & pr.genome.isin(ok_g)].copy().reset_index(drop=True)
called["unverified"] = called.verification == "unverified"
L = L50.reset_index(drop=True); lby = {k: v for k, v in L.groupby(["genome", "contig"])}
cov = []
for r in called.itertuples():
    g = lby.get((r.genome, r.contig))
    cov.append([] if g is None else g[(g.start <= r.end) & (g.end >= r.start)].array.tolist())
called["arrays"] = cov
asize = arr.set_index("array").n
called["n_loci_covered"] = [sum(1 for _ in x) for x in cov]
called["array_ids"] = [sorted(set(x)) for x in cov]
called["order"] = called.genome.map(meta["ord"]); called["sp"] = called.genome.map(meta.sp)
arr["covered"] = arr.array.isin({a for x in called.array_ids for a in x})
arr["covered_unver"] = arr.array.isin({a for r in called[called.unverified].itertuples() for a in r.array_ids})
called["array_n_max"] = [max([asize[a] for a in x], default=0) for x in called.array_ids]
_hx = arr.set_index("array").hx; _fl = arr.set_index("array").flag
called["hx_arr"] = [any(_hx[a] for a in x) if x else False for x in called.array_ids]
called["any_flag_arr"] = [any(_fl[a] for a in x) if x else False for x in called.array_ids]
n_unv = int(called.unverified.sum())
S0 = len(called)
narr = len({a for x in called.array_ids for a in x}) + int((called.n_loci_covered == 0).sum())
unc = arr[arr.flag & ~arr.covered]
E_unc = arr[~arr.covered].p10.sum()          # all uncovered arrays, chance of being flagged
unc_hx = unc[unc.hx]
S = [("S0 current PR calls", S0, f"{n_unv} unverified, {S0 - n_unv} other evidence"),
     ("S1 one call per array", narr, f"{S0} -> {narr}"),
     ("S2 add every uncovered CAAX-flagged array", len(unc), f"in {unc.genome.nunique()} genomes; expected flagged by chance among all uncovered arrays {E_unc:.0f} (within-contig null); excess {len(unc) - E_unc:.0f}"),
     ("S2b add uncovered flagged arrays with >=2 loci", int((unc.n >= 2).sum()), f"chance for size>=2 uncovered: {arr[~arr.covered & (arr.n >= 2)].p10.sum():.0f}"),
     ("S2c add uncovered flagged arrays with precursor homology (Hx<=10kb)", len(unc_hx), "flag AND tblastn precursor homology"),
     ("S4 unverified calls in an array with precursor homology", int(called[called.unverified & called.hx_arr].shape[0]), f"of {n_unv}"),
     ("S5 unverified calls whose array is a single locus", int(((called.array_n_max == 1) & called.unverified).sum()), f"of {n_unv}"),
     ("S6 unverified calls in arrays of >=2 loci", int(((called.array_n_max >= 2) & called.unverified).sum()), f"of {n_unv}")]
SC = pd.DataFrame(S, columns=["scenario", "n", "note"]); T(SC, "v_scenarios.tsv")
rows = []
for o in ORDERS + ["other orders", "ALL"]:
    c = called if o == "ALL" else called[called.order == o]
    a = arr if o == "ALL" else arr[arr.order == o]
    u = a[a.flag & ~a.covered]
    gset = gen[gen.ok & ((gen.ord == o) | (o == "ALL"))]
    rows.append(dict(order=o, genomes=len(gset), calls=len(c), unverified=int(c.unverified.sum()), one_array_per_call=len({x for l in c.array_ids for x in l}) + int((c.n_loci_covered == 0).sum()),
                     unverified_single_locus_array=int(((c.array_n_max == 1) & c.unverified).sum()), unverified_array_ge2=int(((c.array_n_max >= 2) & c.unverified).sum()),
                     unverified_with_Hx_array=int((c.unverified & c.hx_arr).sum()),
                     flagged_arrays=int(a.flag.sum()), flagged_uncovered=len(u), uncovered_expected_by_chance=round(a[~a.covered].p10.sum(), 1),
                     uncovered_flagged_with_Hx=int(u.hx.sum()), uncovered_flagged_size_ge2=int((u.n >= 2).sum()), genomes_with_uncovered_flagged=u.genome.nunique(),
                     genomes_with_call=c.genome.nunique()))
CS = pd.DataFrame(rows); T(CS, "v_calls_and_uncovered_by_order.tsv")
# species bootstrap of the key per-genome rates
G = gen[gen.ok].copy()
pg = called.groupby("genome").size().rename("calls"); G["calls"] = G.genome.map(pg).fillna(0)
G["unc"] = G.genome.map(unc.groupby("genome").size()).fillna(0)
G["has_call"] = (G.calls > 0).astype(int); G["has_unc"] = (G.unc > 0).astype(int); G["one"] = 1
G["call_or_unc"] = ((G.calls > 0) | (G.unc > 0)).astype(int)
rows = []
for o in ORDERS + ["other orders", "ALL"]:
    s = G if o == "ALL" else G[G.ord == o]
    rows.append(dict(order=o, genomes=len(s), species=s.sp.nunique(),
                     pct_genomes_with_call=fmt(tuple(100 * x for x in ratio_ci(s, "has_call", "one")), 1),
                     pct_genomes_with_uncovered_flagged_array=fmt(tuple(100 * x for x in ratio_ci(s, "has_unc", "one")), 1),
                     pct_genomes_call_or_uncovered=fmt(tuple(100 * x for x in ratio_ci(s, "call_or_unc", "one")), 1),
                     calls_per_genome=fmt(ratio_ci(s, "calls", "one"), 2), uncovered_flagged_per_genome=fmt(ratio_ci(s, "unc", "one"), 2)))
T(pd.DataFrame(rows), "v_genome_rates_by_order.tsv")

# ================================================================== F  withheld loci inside arrays
w = pr[(pr.status != "called") & (pr.family == "PR") & pr.genome.isin(ok_g)].copy()
wc = []
for r in w.itertuples():
    g = lby.get((r.genome, r.contig))
    wc.append([] if g is None else sorted(set(g[(g.start <= r.end) & (g.end >= r.start)].array)))
w["array_ids"] = wc
wa = defaultdict(set)
for r in w.itertuples():
    for a in r.array_ids: wa[a].add(r.status)
arr["withheld_bar"] = arr.array.map(lambda a: "withheld:modelled_gene_bar" in wa.get(a, ()))
arr["withheld_floor"] = arr.array.map(lambda a: "withheld:below_fraction_floor" in wa.get(a, ()))
arr["withheld_any"] = arr.withheld_bar | arr.withheld_floor
rows = []
for lab, m in (("flagged, covered by a call", arr.flag & arr.covered), ("flagged, no call", arr.flag & ~arr.covered), ("unflagged, no call", ~arr.flag & ~arr.covered)):
    s = arr[m]
    rows.append(dict(arrays=lab, n=len(s), with_withheld_PR_locus=int(s.withheld_any.sum()), with_modelled_gene_bar=int(s.withheld_bar.sum()),
                     with_fraction_floor=int(s.withheld_floor.sum()), size_ge2=int((s.n >= 2).sum())))
T(pd.DataFrame(rows), "v_withheld_in_arrays.tsv")
# loci in called arrays but not inside any call (siblings)
arr["n_in_call"] = 0
cl = defaultdict(int)
for r in called.itertuples():
    g = lby.get((r.genome, r.contig))
    if g is None: continue
    for a in g[(g.start <= r.end) & (g.end >= r.start)].array: cl[a] += 1
arr["n_in_call"] = arr.array.map(cl).fillna(0).astype(int)
cv = arr[arr.covered]
P(f"\ncalled arrays {len(cv)}; loci {int(cv.n.sum())}; loci inside a call {int(cv.n_in_call.sum())}; sibling loci not inside any call {int((cv.n - cv.n_in_call).sum())};",
  f"arrays with all members inside a call {int((cv.n == cv.n_in_call).sum())}")
arr.to_csv("v_arrays.tsv.gz", sep="\t", index=False)
called.drop(columns=["arrays"]).assign(array_ids=called.array_ids.map(lambda a: ",".join(map(str, a)))).to_csv("v_calls.tsv.gz", sep="\t", index=False)
P("\ndone")
