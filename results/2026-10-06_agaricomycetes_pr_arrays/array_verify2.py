#!/usr/bin/env python3
"""Second pass of the array verification (reads v_arrays.tsv.gz, v_calls.tsv.gz written by array_verify.py).
 1 chance among uncovered arrays (species bootstrap), by order, size, precursor homology; the old per-order-share shortcut
 2 precursor-homology chance (tblastn Hx within 10 kb) and flag AND Hx
 3 identity of the receptor in withheld loci inside uncovered flagged arrays vs called unverified calls
 4 STE20: stratified by array size and by contig length
 5 per-order organisation with species intervals
 6 options table (before/after call counts)"""
import glob, os
from collections import defaultdict
import numpy as np, pandas as pd
exec(open("array_verify.py").read().split("# ================================================================== A  sample")[0].split("# ------------------------------------------------------------------ load")[0].split('RNG = ')[0]) if False else None
RNG = np.random.default_rng(11); NB = 2000
OUT = open("v2_summary.txt", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.write(s + "\n")
def T(df, name):
    df.to_csv(name, sep="\t", index=False); P(f"\n== {name}"); P(df.to_string(index=False))
ORDERS = ["Agaricales", "Boletales", "Polyporales", "Cantharellales", "Russulales", "Hymenochaetales", "Auriculariales", "other orders"]
arr = pd.read_csv("v_arrays.tsv.gz", sep="\t")
calls = pd.read_csv("v_calls.tsv.gz", sep="\t")
gen = pd.read_csv("agari_genomes.tsv", sep="\t")
def ratio_ci(df, num, den, key="sp", nb=NB):
    g = df.groupby(key)[[num, den]].sum(); v = g.values; n = len(v)
    pt = v[:, 0].sum() / v[:, 1].sum()
    idx = RNG.integers(0, n, size=(nb, n)); r = v[idx, 0].sum(1) / np.maximum(v[idx, 1].sum(1), 1e-9)
    return pt, *np.percentile(r, [2.5, 97.5])
def f3(t): return f"{t[0]:.2f} ({t[1]:.2f}-{t[2]:.2f})"

# ---- 2: Hx chance.  Hx positions from cand.tsv (tblastn E<=1), same within-contig analytic null.
cand = pd.concat([pd.read_csv(f, sep="\t").assign(genome=os.path.basename(f)[:-9]) for f in glob.glob("scan_out/*.cand.tsv") if os.path.getsize(f) > 30])
hx = cand[cand.kind == "Hx"]
Hp = defaultdict(dict)
for g, d in hx[hx.genome.isin(arr.genome.unique())].groupby("genome"):
    for c, e in d.groupby("contig"): Hp[g][c] = np.sort(e.pos.values)
def union_len(iv, lo, hi):
    if len(iv) == 0: return 0.0
    iv = iv[np.argsort(iv[:, 0])]; tot = 0.0; cs = ce = None
    for a, b in iv:
        a = max(a, lo); b = min(b, hi)
        if b <= a: continue
        if cs is None: cs, ce = a, b
        elif a <= ce: ce = max(ce, b)
        else: tot += ce - cs; cs, ce = a, b
    return tot + (ce - cs if cs is not None else 0)
ph = []
for r in arr.itertuples():
    t = Hp.get(r.genome, {}).get(r.contig)
    if t is None: ph.append(0.0); continue
    hi = max(r.contig_len - r.span, 1.0); W = 10_000
    iv = np.stack([t - r.span - W, t + W], 1)
    ph.append(union_len(iv, 0, hi) / hi)
arr["p_hx"] = ph
arr["p_flag_hx"] = arr.p10 * arr.p_hx            # independence between a strict-CAAX ORF and a tblastn hit (assumed)
arr["flag_hx"] = arr.flag & arr.hx
arr.to_csv("v_arrays.tsv.gz", sep="\t", index=False)
P("Hx (tblastn E<=1 within 10 kb): observed arrays", int(arr.hx.sum()), "expected by chance", round(arr.p_hx.sum(), 1))

# ---- 1: uncovered
unc = arr[~arr.covered].copy()
unc["one"] = 1
for c, name in (("flag", "O_flag"), ("p10", "E_flag")): pass
unc["O"] = unc.flag.astype(int); unc["E"] = unc.p10
rows = []
def addrow(label, d):
    if d.O.sum() == 0: return
    pt, lo, hi = ratio_ci(d, "E", "O")
    rows.append(dict(set=label, uncovered_arrays=len(d), flagged=int(d.O.sum()), expected_by_chance=round(d.E.sum(), 1), chance_share=f3((pt, lo, hi)),
                     excess_flagged=round(d.O.sum() - d.E.sum(), 0)))
addrow("ALL", unc)
for o in ORDERS: addrow(o, unc[unc.order == o])
addrow("size 1", unc[unc.n == 1]); addrow("size >=2", unc[unc.n >= 2]); addrow("size >=3", unc[unc.n >= 3])
addrow("one genome per species", unc[unc.one_x if "one_x" in unc else unc["one"] == 1])
T(pd.DataFrame(rows), "v_uncovered_chance.tsv")
# flagged AND Hx: observed vs expected under independence
d = unc.assign(O=unc.flag_hx.astype(int), E=unc.p_flag_hx)
rows = []
for lab, s in (("ALL", d), ("Agaricales", d[d.order == "Agaricales"]), ("Boletales", d[d.order == "Boletales"]), ("Russulales", d[d.order == "Russulales"]), ("size >=2", d[d.n >= 2])):
    if s.O.sum() == 0: continue
    pt, lo, hi = ratio_ci(s, "E", "O")
    rows.append(dict(set="flagged AND Hx, " + lab, observed=int(s.O.sum()), expected_by_chance=round(s.E.sum(), 1), chance_share=f3((pt, lo, hi))))
T(pd.DataFrame(rows), "v_uncovered_flag_and_Hx_chance.tsv")
# the old shortcut: per-order locus-level chance share x uncovered flagged
os_ = pd.read_csv("order_summary.tsv", sep="\t").set_index("order").chance_share
ofl = unc[unc.flag].groupby("order").size()
old = sum(os_.get(o if o != "other orders" else [x for x in os_.index if x.startswith("other")][0], np.nan) * n for o, n in ofl.items())
P(f"\nold shortcut (locus-level order chance share x uncovered flagged arrays) = {old:.0f}; array-level expected among uncovered = {unc.E.sum():.0f}; observed {int(unc.O.sum())}")

# ---- 3: identity of withheld receptor locus in uncovered flagged arrays vs called calls
pr = pd.read_csv("pr_calls.tsv.gz", sep="\t")
ok = set(arr.genome)
w = pr[(pr.status == "withheld:below_fraction_floor") & pr.genome.isin(ok)]
ls = pd.read_csv("loci_all.tsv.gz", sep="\t"); ls = ls[ls.genome.isin(ok)]
# array id of the new partition = v_arrays.array; recompute by overlap on contig coordinates
wid = {}
byc = {k: v for k, v in arr.groupby(["genome", "contig"])}
for r in w.itertuples():
    g = byc.get((r.genome, r.contig))
    if g is None: continue
    m = g[(g.start <= r.end) & (g.end >= r.start)]
    for a in m.array: wid[a] = max(wid.get(a, 0), r.rec_identity)
arr["floor_identity"] = arr.array.map(wid)
x = arr[arr.flag & ~arr.covered & arr.floor_identity.notna()]
c = calls[calls.unverified]
P(f"\nreceptor identity of the withheld-below-floor locus in uncovered flagged arrays: n {len(x)}, median {x.floor_identity.median():.1f}, share >=50 {100*(x.floor_identity>=50).mean():.0f}%, >=60 {100*(x.floor_identity>=60).mean():.0f}%")
P(f"receptor identity of called unverified PR calls: n {len(c)}, median {c.rec_identity.median():.1f}, share >=50 {100*(c.rec_identity>=50).mean():.0f}%, >=60 {100*(c.rec_identity>=60).mean():.0f}%")
P("called unverified: calls with receptor+caax only (no curated pheromone):", int(c.genes_found.eq("pheromone_receptor|caax_precursor").sum()), "of", len(c))
P("uncovered flagged arrays with a floor-withheld locus >=50% identity and size>=2:", int(((x.floor_identity >= 50) & (x.n >= 2)).sum()), " >=50% any size:", int((x.floor_identity >= 50).sum()))
arr.to_csv("v_arrays.tsv.gz", sep="\t", index=False)

# ---- 4: STE20 strata
from statsmodels.stats.contingency_tables import StratifiedTable
def cmh(d, strata_cols, D=250):
    tabs = []
    for k, x in d.groupby(strata_cols):
        if x.flag.any() and (~x.flag).any():
            a = int(x[x.flag][f"STE20_obs{D}"].sum()); b = int(x.flag.sum()) - a; cc = int(x[~x.flag][f"STE20_obs{D}"].sum()); dd = int((~x.flag).sum()) - cc
            tabs.append(np.array([[a, b], [cc, dd]]))
    st = StratifiedTable(tabs)
    return len(tabs), int(sum(t[0, 0] for t in tabs)), int(sum(t[0].sum() for t in tabs)), int(sum(t[1, 0] for t in tabs)), int(sum(t[1].sum() for t in tabs)), st.oddsratio_pooled, st.test_null_odds().pvalue
d = arr[arr.STE20_has].copy()
d["size_bin"] = pd.cut(d.n, [0, 1, 2, 3, 100], labels=["1", "2", "3", "4+"]).astype(str)
d["len_bin"] = pd.qcut(d.contig_len, 4, labels=False, duplicates="drop")
rows = []
for lab, cols in (("genome", ["genome"]), ("genome x array size", ["genome", "size_bin"]), ("genome x contig-length quartile", ["genome", "len_bin"]), ("genome x size x length quartile", ["genome", "size_bin", "len_bin"])):
    n, fn, ft, un, ut, orr, p = cmh(d, cols)
    rows.append(dict(strata=lab, informative_strata=n, flagged_near=fn, flagged_total=ft, unflagged_near=un, unflagged_total=ut, MH_OR=round(orr, 1), p=f"{p:.1e}"))
T(pd.DataFrame(rows), "v_STE20_strata.tsv")
# by array size, simple table
rows = []
for sb in ["1", "2", "3", "4+"]:
    s = d[d.size_bin == sb]
    for lab, g in (("flagged", s[s.flag]), ("unflagged", s[~s.flag])):
        rows.append(dict(array_size=sb, group=lab, arrays=len(g), near_STE20_250kb=int(g.STE20_obs250.sum()), near_pct=round(100 * g.STE20_obs250.mean(), 1) if len(g) else np.nan,
                         expected_genome_null=round(g.STE20_pg250.sum(), 1)))
T(pd.DataFrame(rows), "v_STE20_by_size.tsv")
# excess near STE20 among flagged: how much would a chance-flagged fraction dilute?  flagged near rate vs flagged genome-chance share
# genomes: how many genomes have a flagged array near STE20, and of those, is it one array
g1 = arr[arr.flag & (arr.STE20_obs250)].groupby("genome").size()
P(f"genomes with a flagged array within 250 kb of STE20: {len(g1)} of {arr.genome.nunique()} ({100*len(g1)/arr.genome.nunique():.0f}%)")
# does STE20 proximity separate covered vs uncovered flagged arrays? (a check that calls are not simply the STE20 ones)
f = arr[arr.flag & arr.STE20_has]
for lab, g in (("flagged, covered by a call", f[f.covered]), ("flagged, uncovered", f[~f.covered]), ("flagged, uncovered, size>=2", f[~f.covered & (f.n >= 2)]), ("flagged, uncovered, size 1", f[~f.covered & (f.n == 1)])):
    P(f"  {lab}: {len(g)} arrays, near STE20 250 kb {int(g.STE20_obs250.sum())} ({100*g.STE20_obs250.mean():.1f}%)")
u = arr[~arr.flag & ~arr.covered & arr.STE20_has]
P(f"  unflagged uncovered: {len(u)}, near STE20 {int(u.STE20_obs250.sum())} ({100*u.STE20_obs250.mean():.1f}%)")

# ---- 5: per-order organisation with species intervals
gt = pd.read_csv("genome_table.tsv", sep="\t")
gt = gt[gt.scanned & gt.qpass & (gt.contigs <= 5000)].copy()
gt["sp"] = gt.species.fillna("").str.split().str[:2].str.join(" ")
gt["ord"] = gt.order.where(gt.order.isin(ORDERS), "other orders")
gt = gt[["genome", "species", "order", "sp", "ord"]]
ag = arr.groupby("genome").agg(n_arrays=("array", "size"), nl=("n", "sum"), largest=("n", "max"), n_flag_arr=("flag", "sum"), n_unc_flag=("flag", lambda s: 0)).reset_index()
ag["n_unc_flag"] = ag.genome.map(arr[arr.flag & ~arr.covered].groupby("genome").size()).fillna(0)
ag["n_cov"] = ag.genome.map(arr[arr.covered].groupby("genome").size()).fillna(0)
ag["n_big"] = ag.genome.map(arr[arr.n >= 3].groupby("genome").size()).fillna(0)
gt = gt.merge(ag, on="genome", how="left").fillna({"n_arrays": 0, "nl": 0, "largest": 0, "n_flag_arr": 0, "n_unc_flag": 0, "n_cov": 0, "n_big": 0})
gt["has2"] = (gt.largest >= 2).astype(int); gt["hasflag"] = (gt.n_flag_arr > 0).astype(int); gt["ge5"] = (gt.nl >= 5).astype(int); gt["one"] = 1
gt["has_call"] = gt.genome.isin(calls.genome).astype(int)
rows = []
for o in ORDERS + ["ALL"]:
    s = gt if o == "ALL" else gt[gt.ord == o]
    rows.append(dict(order=o, genomes=len(s), species=s.sp.nunique(), STE3_loci_per_genome_median=s.nl.median(), q1_q3=f"{s.nl.quantile(.25):.0f}-{s.nl.quantile(.75):.0f}", max=int(s.nl.max()),
                     pct_ge5_loci=f3(tuple(100 * x for x in ratio_ci(s, "ge5", "one"))),
                     arrays_per_genome=f3(ratio_ci(s, "n_arrays", "one")), pct_genomes_array_ge2=f3(tuple(100 * x for x in ratio_ci(s, "has2", "one"))),
                     pct_genomes_flagged_array=f3(tuple(100 * x for x in ratio_ci(s, "hasflag", "one"))),
                     flagged_arrays_per_genome=f3(ratio_ci(s, "n_flag_arr", "one")), pct_genomes_with_PR_call=f3(tuple(100 * x for x in ratio_ci(s, "has_call", "one")))))
T(pd.DataFrame(rows), "v2_order_organisation.tsv")
# genome-level: number of arrays >=3 loci
P("\ngenomes with at least one array of >=3 loci:", int((gt.n_big > 0).sum()), "of", len(gt))
# singletons / 2 / 3+ by order (arrays)
rows = []
for o in ORDERS + ["ALL"]:
    s = arr if o == "ALL" else arr[arr.order == o]
    rows.append(dict(order=o, arrays=len(s), singleton=int((s.n == 1).sum()), size2=int((s.n == 2).sum()), size3=int((s.n == 3).sum()), size4plus=int((s.n >= 4).sum()), max=int(s.n.max()),
                     span_kb_median_ge2=round(s[s.n >= 2].span.median() / 1000, 1), flagged_pct_size1=round(100 * s[s.n == 1].flag.mean(), 1), flagged_pct_size2plus=round(100 * s[s.n >= 2].flag.mean(), 1)))
T(pd.DataFrame(rows), "v2_array_sizes_by_order.tsv")
P("done")
