#!/usr/bin/env python3
"""Agaricomycetes receptor arrays: organisation per genome, pipeline calls vs arrays, rule scenarios.
Inputs: agari_genomes.tsv, scan_out/ (array_scan.py), pr_calls.tsv + hd_loci.tsv (parse_reports.py),
../2026-10-05_caax_receptor_test/panel_loci.tsv.  Outputs: *.tsv and summary.txt in this directory.
Quality filter for headline tables: BUSCO complete >= 70 and N50 >= 20 kb (column qpass).
Intervals: bootstrap over species (resample species with replacement within order), 2000 draws, seed 1."""
import glob, os, sys
from collections import defaultdict
import numpy as np
import pandas as pd

GAP = int(os.environ.get("ARRAY_GAP", 50000))
RNG = np.random.default_rng(1)
NB = 2000
gen = pd.read_csv("agari_genomes.tsv", sep="\t")
gen["sp"] = gen.species.fillna("").str.split().str[:2].str.join(" ")
scanned = {os.path.basename(f)[:-11] for f in glob.glob("scan_out/*.chance.tsv")}
gen["scanned"] = gen.genome.isin(scanned)

def rd(suffix):
    fs = glob.glob(f"scan_out/*.{suffix}.tsv")
    dfs = []
    for f in fs:
        d = pd.read_csv(f, sep="\t")
        d["genome"] = os.path.basename(f)[:-len(suffix) - 5]
        dfs.append(d)
    return pd.concat(dfs, ignore_index=True)

loci = rd("loci").drop(columns=["asm"], errors="ignore")
reg = rd("region")
chance = rd("chance").drop(columns=["asm"], errors="ignore")
pr = pd.read_csv("pr_calls.tsv.gz", sep="\t")
hd = pd.read_csv("hd_loci.tsv.gz", sep="\t")
loci["flag"] = loci.nT_10kb > 0
loci["hx"] = loci.nHx_10kb > 0

# ---- arrays: single linkage, same contig, gap <= GAP
def make_arrays(df, gap):
    df = df.sort_values(["genome", "contig", "start"]).reset_index(drop=True)
    ids, aid, last, cur_end = [], 0, None, -1
    for g, c, st, en in zip(df.genome, df.contig, df.start, df.end):
        if (g, c) != last or st - cur_end > gap:
            aid += 1; cur_end = en
        else:
            cur_end = max(cur_end, en)
        last = (g, c); ids.append(aid)
    df["array"] = ids
    return df
loci = make_arrays(loci, GAP)
arr = loci.groupby("array").agg(genome=("genome", "first"), contig=("contig", "first"), start=("start", "min"), end=("end", "max"),
                                n=("start", "size"), n_flag=("flag", "sum"), n_hx=("hx", "sum"), contig_len=("contig_len", "first")).reset_index()
arr["span"] = arr.end - arr.start
arr["has_flag"] = arr.n_flag > 0

# ---- calls
bsub = pr[(pr.status == "called") & (pr.family != "PR")]
pr = pr[pr.family == "PR"]
called = pr[pr.status == "called"].copy()
withheld = pr[pr.status != "called"].copy()
def overlap_loci(df_calls, with_flag=False):
    rows = []
    L = {k: v for k, v in loci.groupby(["genome", "contig"])}
    for r in df_calls.itertuples():
        g = L.get((r.genome, r.contig))
        if g is None:
            rows.append((r.Index, [], [])); continue
        m = g[(g.start <= r.end) & (g.end >= r.start)]
        rows.append((r.Index, m.index.tolist(), m.array.tolist()))
    return rows
cov = overlap_loci(called)
called["n_loci_covered"] = [len(x[1]) for x in cov]
called["arrays"] = [sorted(set(x[2])) for x in cov]
called["n_flag_covered"] = [int(loci.loc[x[1], "flag"].sum()) if x[1] else 0 for x in cov]
called["unverified"] = called.verification == "unverified"
loci["in_call"] = False
loci.loc[[i for x in cov for i in x[1]], "in_call"] = True
wcov = overlap_loci(withheld)
loci["in_withheld"] = False
loci.loc[[i for x in wcov for i in x[1]], "in_withheld"] = True

# ---- HD positions (pipeline called HD / bLocus loci) and STE20 candidate (top-scoring hit)
hdc = hd[hd.status == "called"]
hd_by = {k: v for k, v in hdc.groupby("genome")}
ste = reg[reg["class"] == "STE20"].copy()
top = ste.sort_values("score", ascending=False).groupby("genome").head(1).set_index("genome") if len(ste) else pd.DataFrame()
def dist(a0, a1, b0, b1):
    return max(0, max(a0, b0) - min(a1, b1))
d_hd, same_hd, d_ste, nst = [], [], [], []
for r in arr.itertuples():
    h = hd_by.get(r.genome)
    best = np.inf
    if h is not None:
        s = h[h.contig == r.contig]
        for x in s.itertuples():
            best = min(best, dist(r.start, r.end, x.start, x.end))
    d_hd.append(best)
    best = np.inf
    if r.genome in top.index:
        t = top.loc[r.genome]
        if t.contig == r.contig:
            best = dist(r.start, r.end, t.start, t.end)
    d_ste.append(best)
arr["d_HD"] = d_hd; arr["d_STE20"] = d_ste
arr["has_HDcall"] = arr.genome.isin(hd_by.keys())
arr["has_STE20"] = arr.genome.isin(top.index)
arr["covered"] = arr.array.isin({a for x in cov for a in x[2]}) if "array" in arr else False
arr["covered"] = arr["array"].isin({a for x in cov for a in x[2]})
unver_arrays = {a for r, x in zip(called.itertuples(), cov) if r.unverified for a in x[2]}
arr["covered_unverified"] = arr["array"].isin(unver_arrays)
loci_cov = defaultdict(int)
for r, x in zip(called.itertuples(), cov):
    pass

# ---- per genome table
gl = loci.groupby("genome").agg(n_ste3=("start", "size"), n_flag=("flag", "sum"), n_in_call=("in_call", "sum"),
                                n_flag_in_call=("in_call", lambda s: 0)).reset_index()
gl["n_flag_in_call"] = loci[loci.flag & loci.in_call].groupby("genome").size().reindex(gl.genome).fillna(0).astype(int).values
ga = arr.groupby("genome").agg(n_arrays=("array", "size"), largest_array=("n", "max"), n_flag_arrays=("has_flag", "sum"),
                               n_arrays_called=("covered", "sum"), max_span=("span", "max")).reset_index()
gc = called.groupby("genome").agg(n_pr_calls=("start", "size"), n_unverified=("unverified", "sum")).reset_index()
G = gen.merge(chance[["genome", "p_random"]], on="genome", how="left").merge(gl, on="genome", how="left").merge(ga, on="genome", how="left").merge(gc, on="genome", how="left")
G["n_hd_calls"] = G.genome.map(hdc.groupby("genome").size()).fillna(0).astype(int)
for c in ["n_ste3", "n_flag", "n_in_call", "n_flag_in_call", "n_arrays", "largest_array", "n_flag_arrays", "n_arrays_called", "n_pr_calls", "n_unverified"]:
    G[c] = G[c].fillna(0).astype(int)
G.loc[~G.scanned, ["n_ste3", "n_flag"]] = np.nan
G.to_csv("genome_table.tsv", sep="\t", index=False)
arr.to_csv("arrays.tsv", sep="\t", index=False)
called.drop(columns=["arrays"]).assign(arrays=called.arrays.map(lambda a: ",".join(map(str, a)))).to_csv("pr_calls_vs_loci.tsv", sep="\t", index=False)
loci.to_csv("loci_all.tsv.gz", sep="\t", index=False)

# ---- bootstrap helpers (species within order)
def boot_ratio(df, num, den, nb=NB):
    """ratio sum(num)/sum(den) with species-level bootstrap."""
    sp = df.groupby("sp")[[num, den]].sum()
    n = len(sp)
    if n == 0 or sp[den].sum() == 0:
        return np.nan, np.nan, np.nan
    pt = sp[num].sum() / sp[den].sum()
    v = sp.values
    idx = RNG.integers(0, n, size=(nb, n))
    nums, dens = v[idx, 0].sum(1), v[idx, 1].sum(1)
    ok = dens > 0
    r = nums[ok] / dens[ok]
    return pt, *np.percentile(r, [2.5, 97.5])
def boot_stat(vals_by_sp, f, nb=NB):
    keys = list(vals_by_sp); n = len(keys)
    pt = f(np.concatenate([vals_by_sp[k] for k in keys]))
    out = []
    for _ in range(nb):
        pick = RNG.integers(0, n, n)
        out.append(f(np.concatenate([vals_by_sp[keys[i]] for i in pick])))
    return pt, *np.percentile(out, [2.5, 97.5])
def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)

Q = G[G.scanned & G.qpass].copy()
Q["exp_flag"] = Q.n_ste3 * Q.p_random
Q["has_flag_arr"] = (Q.n_flag > 0).astype(int); Q["has_call"] = (Q.n_pr_calls > 0).astype(int)
Q["has_unver"] = (Q.n_unverified > 0).astype(int); Q["one"] = 1
Q["has_ste3"] = (Q.n_ste3 > 0).astype(int)
Q["two_plus"] = (Q.n_ste3 >= 2).astype(int)
Q["no_call_flag"] = ((Q.n_pr_calls == 0) & (Q.n_flag > 0)).astype(int)
rows = []
orders = ["Agaricales", "Boletales", "Polyporales", "Cantharellales", "Russulales", "Hymenochaetales", "Auriculariales"]
rest = [o for o in Q.order.unique() if o not in orders]
for name, sub in [(o, Q[Q.order == o]) for o in orders] + [("other orders (n=%d)" % len(rest), Q[Q.order.isin(rest)]), ("ALL", Q)]:
    allg = G[(G.order == name) | ((G.order.isin(rest)) & name.startswith("other")) | (name == "ALL")]
    n_all = len(allg); n_q = len(sub)
    ste = sub.n_ste3
    ex, exlo, exhi = boot_ratio(sub, "exp_flag", "n_flag")  # chance share of flagged loci
    rows.append(dict(order=name, genomes_all=n_all, genomes_scanned_qpass=n_q, species=sub.sp.nunique(),
                     ste3_median=ste.median(), ste3_q1=ste.quantile(.25), ste3_q3=ste.quantile(.75), ste3_mean=round(ste.mean(), 2),
                     pct_0=round(100 * (ste == 0).mean(), 1), pct_ge5=round(100 * (ste >= 5).mean(), 1), ste3_max=ste.max(),
                     pct_genomes_cluster2=round(100 * (sub.largest_array >= 2).mean(), 1),
                     largest_array_median=sub.largest_array.median(), pct_genomes_flagarray=round(100 * sub.has_flag_arr.mean(), 1),
                     flagged_loci=int(sub.n_flag.sum()), expected_flagged_chance=round(sub.exp_flag.sum(), 1),
                     chance_share=round(ex, 3), chance_lo=round(exlo, 3), chance_hi=round(exhi, 3),
                     pct_genomes_PRcall=round(100 * sub.has_call.mean(), 1), pct_genomes_unverified_call=round(100 * sub.has_unver.mean(), 1),
                     genomes_flag_no_call=int(sub.no_call_flag.sum()), PRcalls=int(sub.n_pr_calls.sum()), unverified_calls=int(sub.n_unverified.sum())))
O = pd.DataFrame(rows); O.to_csv("order_summary.tsv", sep="\t", index=False)
print(O.T.to_string()); print("B-sublocus calls (Balpha/Bbeta, Schizophyllum record):", len(bsub), bsub.verification.value_counts().to_dict())

# =====================================================================================
# 2. array organisation (quality-pass, scanned genomes)
# =====================================================================================
A = arr.merge(Q[["genome", "order", "sp"]], on="genome")          # arrays of quality-pass genomes only
out = []
groups = [(o, A[A.order == o]) for o in orders] + [("other orders", A[A.order.isin(rest)]), ("ALL", A)]
for name, a in groups:
    f = a[a.has_flag]
    g_ok = Q if name == "ALL" else (Q[Q.order.isin(rest)] if name == "other orders" else Q[Q.order == name])
    hdok = f[f.has_HDcall]; stok = f[f.has_STE20]
    out.append(dict(order=name, genomes=len(g_ok), arrays=len(a), arrays_per_genome=round(len(a) / max(1, len(g_ok)), 2),
                    singletons=int((a.n == 1).sum()), size2_3=int(a.n.between(2, 3).sum()), size4plus=int((a.n >= 4).sum()),
                    array_size_max=int(a.n.max()) if len(a) else 0, span_median_kb=round(a[a.n > 1].span.median() / 1000, 1) if (a.n > 1).any() else np.nan,
                    flagged_arrays=len(f), flagged_arrays_per_genome=round(len(f) / max(1, len(g_ok)), 2),
                    flagged_size_median=f.n.median() if len(f) else np.nan, flagged_nflag_median=f.n_flag.median() if len(f) else np.nan,
                    flagged_with_Hx=int((f.n_hx > 0).sum()),
                    flagged_arrays_in_HDcall_genome=len(hdok),
                    HD_same_contig=int((hdok.d_HD < np.inf).sum()), HD_within_100kb=int((hdok.d_HD <= 100000).sum()),
                    flagged_arrays_in_STE20_genome=len(stok), STE20_same_contig=int((stok.d_STE20 < np.inf).sum()),
                    STE20_within_100kb=int((stok.d_STE20 <= 100000).sum()), STE20_within_250kb=int((stok.d_STE20 <= 250000).sum())))
AO = pd.DataFrame(out); AO.to_csv("array_organisation.tsv", sep="\t", index=False)
print(AO.T.to_string())

# unflagged arrays for contrast (descriptive only; flag != label)
ctr = []
for name, a in groups:
    for lab, s in (("flagged", a[a.has_flag]), ("unflagged", a[~a.has_flag])):
        sh = s[s.has_STE20]; hh = s[s.has_HDcall]
        ctr.append(dict(order=name, group=lab, arrays=len(s), size_median=s.n.median(), STE20_genome_arrays=len(sh),
                        STE20_within_250kb=int((sh.d_STE20 <= 250000).sum()), HD_genome_arrays=len(hh), HD_same_contig=int((hh.d_HD < np.inf).sum())))
pd.DataFrame(ctr).to_csv("flagged_vs_unflagged_arrays.tsv", sep="\t", index=False)

# =====================================================================================
# 3. what the calls cover  (qpass and all scanned genomes)
# =====================================================================================
from collections import Counter
ca = arr.set_index("array")
def coverage_table(GS, tag):
    cq = called.merge(GS[["genome", "order", "sp"]], on="genome")
    cq["array_size"] = cq.arrays.map(lambda a: max([ca.loc[x, "n"] for x in a], default=0))
    cq["frac_array_covered"] = cq.n_loci_covered / cq.array_size.replace(0, np.nan)
    cq.drop(columns=["arrays"]).to_csv(f"calls_{tag}.tsv", sep="\t", index=False)
    A_ = arr.merge(GS[["genome", "order", "sp"]], on="genome")
    ncall = Counter(x for l in cq.arrays for x in l)
    rows = []
    for name, sub in [(o, cq[cq.order == o]) for o in orders] + [("other orders", cq[cq.order.isin(rest)]), ("ALL", cq)]:
        gq = GS if name == "ALL" else (GS[GS.order.isin(rest)] if name == "other orders" else GS[GS.order == name])
        a = A_ if name == "ALL" else (A_[A_.order.isin(rest)] if name == "other orders" else A_[A_.order == name])
        cv = a[a.covered]
        wh = withheld.merge(gq[["genome"]], on="genome")
        rows.append(dict(order=name, genomes=len(gq), PR_genomes=int(gq.has_call.sum()) if "has_call" in gq else int((gq.n_pr_calls > 0).sum()),
                         calls=len(sub), unverified=int(sub.unverified.sum()), verified_other=int((~sub.unverified).sum()),
                         calls_no_scan_locus=int((sub.n_loci_covered == 0).sum()),
                         loci_per_call_median=sub.n_loci_covered.median(), calls_1locus=int((sub.n_loci_covered == 1).sum()),
                         calls_2plus_loci=int((sub.n_loci_covered >= 2).sum()),
                         array_size_of_called_median=sub.array_size.median(), frac_array_covered_median=round(sub.frac_array_covered.median(), 2),
                         calls_in_array_with_uncovered_loci=int((sub.frac_array_covered < 1).sum()),
                         arrays_with_gt1_call=int(sum(1 for x in cv.array if ncall.get(x, 0) > 1)),
                         arrays_total=len(a), arrays_flagged=int(a.has_flag.sum()), flagged_arrays_covered=int((a.has_flag & a.covered).sum()),
                         flagged_arrays_uncovered=int((a.has_flag & ~a.covered).sum()),
                         unflagged_arrays_covered=int((~a.has_flag & a.covered).sum()),
                         loci_total=int(a.n.sum()), loci_in_called_arrays=int(cv.n.sum()), loci_inside_calls=int(sub.n_loci_covered.sum()),
                         withheld_PR_loci=len(wh), withheld_modelled_gene_bar=int((wh.status == "withheld:modelled_gene_bar").sum()),
                         withheld_below_fraction_floor=int((wh.status == "withheld:below_fraction_floor").sum())))
    T_ = pd.DataFrame(rows); T_.to_csv(f"call_coverage_{tag}.tsv", sep="\t", index=False)
    return cq, T_
cq, CO = coverage_table(Q, "qpass")
print(CO.T.to_string())
cq_all, CO_all = coverage_table(G[G.scanned], "allscanned")

# =====================================================================================
# 4. curated panel (6 Agaricomycete genomes of the PR #33 study): where do the curated B receptors sit?
#    Independent labels (genetically mapped): Coprinopsis B43 and Schizophyllum bar3/bbr2 only; Trametes, Grifola,
#    Russula, Heterobasidion records were chosen from CAAX positional evidence (not independent).
# =====================================================================================
pnl = pd.read_csv("../2026-10-05_caax_receptor_test/panel_loci.tsv", sep="\t")
pnl = pnl[pnl.agarico.astype(str) == "True"]
prow = []
for r in pnl.itertuples():
    g = loci[(loci.genome == r.asm) & (loci.contig == r.contig) & (loci.start <= r.end) & (loci.end >= r.start)]
    if g.empty:
        prow.append(dict(asm=r.asm, mating=r.mating, independent=r.independent, found=False)); continue
    x = g.iloc[0]; a = arr[arr["array"] == x["array"]].iloc[0]
    prow.append(dict(asm=r.asm, mating=r.mating, independent=r.independent, found=True, contig=r.contig, start=r.start, array=x["array"],
                     array_n=a.n, array_flagged=bool(a.has_flag), array_n_flag=int(a.n_flag), locus_flag=bool(x.flag), in_call=bool(x.in_call),
                     array_covered=bool(a.covered), d_HD=a.d_HD, d_STE20=a.d_STE20))
PN = pd.DataFrame(prow); PN.to_csv("panel_agaricomycetes.tsv", sep="\t", index=False)
if len(PN) and PN.found.any():
    pg = PN[PN.found].groupby(["asm", "mating"]).agg(loci=("array", "size"), arrays=("array", "nunique"), max_array=("array_n", "max"),
                                                     flagged_loci=("locus_flag", "sum"), in_flagged_array=("array_flagged", "sum"),
                                                     in_call=("in_call", "sum"), array_called=("array_covered", "sum")).reset_index()
    pg.to_csv("panel_per_genome.tsv", sep="\t", index=False); print(pg.to_string())
    # largest array of each panel genome: does it hold the curated receptors?
    lar = []
    for asm in PN.asm.unique():
        a = arr[arr.genome == asm].sort_values(["n", "n_flag"], ascending=False)
        mating_arrays = set(PN[(PN.asm == asm) & (PN.mating == "mating") & PN.found].array)
        top_a = a.iloc[0]["array"] if len(a) else None
        fl = a[a.has_flag]
        lar.append(dict(asm=asm, arrays=len(a), largest_array_n=int(a.n.max()), largest_holds_mating=top_a in mating_arrays,
                        flagged_arrays=len(fl), flagged_hold_mating=int(fl.array.isin(mating_arrays).sum()),
                        mating_arrays=len(mating_arrays), mating_array_ids_are_flagged=int(a[a.array.isin(mating_arrays)].has_flag.sum())))
    pd.DataFrame(lar).to_csv("panel_array_rank.tsv", sep="\t", index=False); print(pd.DataFrame(lar).to_string())

# =====================================================================================
# 5. rule scenarios, before/after (call and genome counts)
# =====================================================================================
chance_share = {r.order: r.chance_share for r in O.itertuples()}
hi_share = {r.order: r.chance_hi for r in O.itertuples()}
S = []
def scen(label, scope, calls, extra=""):
    S.append(dict(scenario=label, scope=scope, n=calls, note=extra))
for tag, GS, cqx in (("qpass", Q, cq), ("allscanned", G[G.scanned], cq_all)):
    nloc = cqx.groupby("order").size()
    n_unv = int(cqx.unverified.sum())
    scen("S0 current: PR calls (loci)", tag, len(cqx), f"{n_unv} unverified, {len(cqx) - n_unv} verified by other evidence")
    # S1: one call per array
    arr_ids = {x for l in cqx.arrays for x in l}
    scen("S1 one call per receptor array (merge calls sharing an array; calls overlapping no scan locus kept)", tag,
         len(arr_ids) + int((cqx.n_loci_covered == 0).sum()), f"{len(cqx)} -> {len(arr_ids) + int((cqx.n_loci_covered == 0).sum())}")
    # S2: admit flagged arrays not covered by a call
    A_ = arr.merge(GS[["genome", "order"]], on="genome")
    unc = A_[A_.has_flag & ~A_.covered]
    exp_fp = sum(chance_share.get(o, np.nan) * n for o, n in unc.groupby("order").size().items() if not np.isnan(chance_share.get(o, np.nan)))
    scen("S2 admit CAAX-flagged arrays that no call covers (new arrays)", tag, len(unc),
         f"in {unc.genome.nunique()} genomes; approx chance share (per-order flagged-locus level) {exp_fp:.0f}")
    # S3: drop `unverified` where chance-share upper bound <= 0.30 (excess lower bound >= 0.7)
    ok_orders = [o for o in orders if hi_share.get(o, 1) <= 0.30]
    n3 = int(cqx[cqx.unverified & cqx.order.isin(ok_orders)].shape[0])
    scen("S3 lift `unverified` in orders with chance-share upper 95% <= 0.30", tag, n3, "orders: " + ", ".join(ok_orders))
    # S4: T and Hx tier: flagged arrays with precursor homology, unverified calls in them
    hx_arr = set(arr[(arr.n_hx > 0) & arr.has_flag].array)
    n4 = int(sum(1 for r in cqx.itertuples() if r.unverified and any(x in hx_arr for x in r.arrays)))
    scen("S4 calls whose array also has tblastn precursor homology (T and Hx tier)", tag, n4, f"of {n_unv} unverified calls")
    # S5: single-locus arrays with CAAX: weakest tier
    sing = int(sum(1 for r in cqx.itertuples() if r.unverified and all(ca.loc[x, "n"] == 1 for x in r.arrays) and r.arrays))
    scen("S5 unverified calls whose array has one STE3 locus (no array support)", tag, sing, f"of {n_unv}")
SC = pd.DataFrame(S); SC.to_csv("scenarios.tsv", sep="\t", index=False); print(SC.to_string())

# =====================================================================================
# 6. call composition (all 1,276 PR calls of Agaricomycetes, all genomes, no quality filter)
# =====================================================================================
cc = called.merge(gen[["genome", "order", "qpass"]], on="genome")
cc["has_caax"] = cc.genes_found.str.contains("caax_precursor")
cc["has_curated_pher"] = cc.genes_found.str.contains("pheromone_B4|fungal_mating_type_pheromone")
cc["hd_in_genome"] = cc.genome.isin(hd_by.keys())
comp = cc.groupby("order").agg(calls=("start", "size"), unverified=("unverified", "sum"), caax_in_call=("has_caax", "sum"),
                               median_receptor_rows=("n_receptor_rows", "median"), calls_ge2_receptor_rows=("n_receptor_rows", lambda s: int((s >= 2).sum())),
                               median_loci_overlapped=("n_loci_covered", "median"), qpass_calls=("qpass", "sum"),
                               genomes=("genome", "nunique"), genomes_with_HD_call=("hd_in_genome", lambda s: 0)).reset_index()
hdg = cc.drop_duplicates("genome").groupby("order").hd_in_genome.sum()
comp["genomes_with_HD_call"] = comp.order.map(hdg)
comp.sort_values("calls", ascending=False).to_csv("call_composition.tsv", sep="\t", index=False)
print(comp.sort_values("calls", ascending=False).to_string())
vo = cc[~cc.unverified]
print("verified-other genes_found:", vo.genes_found.value_counts().head(8).to_dict())
print("unverified genes_found:", cc[cc.unverified].genes_found.value_counts().head(5).to_dict())
print("caax motifs unverified:", cc[cc.unverified].caax_motif.str.extract(r"([A-Z]{4})")[0].value_counts().head(8).to_dict())
# genome-level: PR-only (no HD call) vs both
gg = G[G.n_pr_calls > 0]
print("PR-call genomes:", len(gg), "with HD call:", int((gg.n_hd_calls > 0).sum()))
