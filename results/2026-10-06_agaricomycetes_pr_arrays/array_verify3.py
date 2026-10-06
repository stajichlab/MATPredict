#!/usr/bin/env python3
"""Third pass: scenario sensitivity to the array gap, call tiers for the options table, per-genome flagged-array structure."""
import glob, os
import numpy as np, pandas as pd
OUT = open("v3_summary.txt", "w")
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.write(s + "\n")
def T(df, name):
    df.to_csv(name, sep="\t", index=False); P(f"\n== {name}"); P(df.to_string(index=False))
arr = pd.read_csv("v_arrays.tsv.gz", sep="\t")
calls = pd.read_csv("v_calls.tsv.gz", sep="\t")
gen = pd.read_csv("agari_genomes.tsv", sep="\t")
ok = set(arr.genome)
def rd(suffix):
    out = []
    for f in glob.glob(f"scan_out/*.{suffix}.tsv"):
        d = pd.read_csv(f, sep="\t")
        if len(d): d["genome"] = os.path.basename(f)[:-len(suffix) - 5]; out.append(d)
    return pd.concat(out, ignore_index=True)
loci = rd("loci"); loci = loci[loci.genome.isin(ok)].copy()
def make_arrays(df, gap):
    df = df.sort_values(["genome", "contig", "start"]).copy()
    key = df.genome + "|" + df.contig
    cm = df.groupby(key, sort=False).end.cummax(); prev = cm.groupby(key, sort=False).shift()
    df["array"] = ((key != key.shift()) | ((df.start - prev) > gap)).cumsum().values
    return df
# ---- scenario counts by gap
rows = []
calls = calls[["genome", "contig", "start", "end", "unverified"]].copy()
for gap in (10_000, 25_000, 50_000, 100_000, 200_000):
    L = make_arrays(loci, gap)
    A = L.groupby("array").agg(genome=("genome", "first"), contig=("contig", "first"), n=("start", "size"), flag=("d_T", lambda s: bool((s <= 10_000).any())))
    lby = {k: v for k, v in L.groupby(["genome", "contig"])}
    cov = []
    for r in calls.itertuples():
        g = lby.get((r.genome, r.contig)); cov.append([] if g is None else sorted(set(g[(g.start <= r.end) & (g.end >= r.start)].array)))
    calls["arrs"] = cov
    cset = {a for x in cov for a in x}
    s1 = len(cset) + sum(1 for x in cov if not x)
    amax = [max([A.n[a] for a in x], default=0) for x in cov]
    unv = calls.unverified.values
    unc = A[A.flag & ~A.index.isin(cset)]
    rows.append(dict(gap_kb=gap // 1000, arrays=len(A), S0_calls=len(calls), S1_one_per_array=s1, S2_uncovered_flagged=len(unc), S2_size_ge2=int((unc.n >= 2).sum()),
                     unverified_single_locus=int(sum(1 for m, u in zip(amax, unv) if u and m == 1)), unverified_array_ge2=int(sum(1 for m, u in zip(amax, unv) if u and m >= 2)),
                     arrays_with_gt1_call=int(sum(1 for a, c in pd.Series([a for x in cov for a in x]).value_counts().items() if c > 1))))
T(pd.DataFrame(rows), "v_scenarios_by_gap.tsv")

# ---- tiers of called arrays (50 kb)
a = arr.set_index("array")
c = pd.read_csv("v_calls.tsv.gz", sep="\t"); c["ids"] = c.array_ids.fillna("").map(lambda s: [int(x) for x in str(s).split(",") if x != ""])
c["nT2"] = [max([a.nT[x] for x in ids], default=0) >= 2 for ids in c.ids]
c["size"] = [max([a.n[x] for x in ids], default=0) for ids in c.ids]
c["hx"] = [any(a.hx[x] for x in ids) for ids in c.ids]
c["support"] = (c["size"] >= 2) | c.hx | c.nT2
u = c[c.unverified]
rows = [dict(tier="unverified calls", n=len(u)),
        dict(tier="  array >= 2 STE3 loci", n=int((u["size"] >= 2).sum())),
        dict(tier="  precursor homology (tblastn) within 10 kb", n=int(u.hx.sum())),
        dict(tier="  >= 2 strict-CAAX ORFs within 10 kb", n=int(u.nT2.sum())),
        dict(tier="  any of the three (array-supported)", n=int(u.support.sum())),
        dict(tier="  none (single locus, one CAAX ORF, no homology)", n=int((~u.support).sum())),
        dict(tier="  none, receptor + CAAX only gene content", n=int((~u.support & u.genes_found.eq("pheromone_receptor|caax_precursor")).sum()))]
T(pd.DataFrame(rows), "v_call_tiers.tsv")
by = u.groupby("order").agg(unverified=("order", "size"), supported=("support", "sum")).assign(unsupported=lambda d: d.unverified - d.supported)
T(by.reset_index(), "v_call_tiers_by_order.tsv")
# chance share by size among flagged arrays (any call status)
P("\nchance share (expected/observed flagged) by array size and CAAX multiplicity, all arrays")
for lab, m in (("size 1", arr.n == 1), ("size 2", arr.n == 2), ("size>=3", arr.n >= 3)):
    s = arr[m]; P(f"  {lab}: flagged {int(s.flag.sum())}, expected {s.p10.sum():.0f}, share {s.p10.sum()/s.flag.sum():.2f}")
s = arr[arr.n == 1]; P(f"  size 1 with >=2 CAAX ORFs: flagged {int(s.flag10_k2.sum())}, expected {s.p10_k2.sum():.1f}")

# ---- per-genome structure of flagged arrays (candidate B locus uniqueness)
g = arr.groupby("genome").agg(flagged=("flag", "sum"), flagged_ge2=("flag", lambda s: 0)).reset_index()
g["flagged_ge2"] = g.genome.map(arr[arr.flag & (arr.n >= 2)].groupby("genome").size()).fillna(0)
g["flagged_ge3"] = g.genome.map(arr[arr.flag & (arr.n >= 3)].groupby("genome").size()).fillna(0)
g["called"] = g.genome.isin(c.genome)
rows = []
for lab, col in (("flagged arrays", "flagged"), ("flagged arrays of >=2 loci", "flagged_ge2"), ("flagged arrays of >=3 loci", "flagged_ge3")):
    vc = g[col].clip(upper=3).value_counts().sort_index()
    rows.append(dict(count_per_genome=lab, zero=int(vc.get(0, 0)), one=int(vc.get(1, 0)), two=int(vc.get(2, 0)), three_plus=int(vc.get(3, 0)), genomes=len(g)))
T(pd.DataFrame(rows), "v_flagged_arrays_per_genome.tsv")
# genomes that would gain a first PR call under S2 variants
unc = arr[arr.flag & ~arr.covered]
nocall = set(g[~g.called].genome)
rows = [dict(rule="S2 all uncovered flagged arrays", arrays=len(unc), genomes_new_first_call=unc[unc.genome.isin(nocall)].genome.nunique(), genomes_extra_call=unc[~unc.genome.isin(nocall)].genome.nunique()),
        dict(rule="S2b size >= 2", arrays=int((unc.n >= 2).sum()), genomes_new_first_call=unc[(unc.n >= 2) & unc.genome.isin(nocall)].genome.nunique(), genomes_extra_call=unc[(unc.n >= 2) & ~unc.genome.isin(nocall)].genome.nunique()),
        dict(rule="S2c flagged + precursor homology", arrays=int(unc.hx.sum()), genomes_new_first_call=unc[unc.hx & unc.genome.isin(nocall)].genome.nunique(), genomes_extra_call=unc[unc.hx & ~unc.genome.isin(nocall)].genome.nunique()),
        dict(rule="S2d size >= 2 AND Hx", arrays=int((unc.hx & (unc.n >= 2)).sum()), genomes_new_first_call=unc[unc.hx & (unc.n >= 2) & unc.genome.isin(nocall)].genome.nunique(), genomes_extra_call=unc[unc.hx & (unc.n >= 2) & ~unc.genome.isin(nocall)].genome.nunique())]
T(pd.DataFrame(rows), "v_S2_genome_effect.tsv")
P("genomes without a PR call:", len(nocall), "of", len(g), "; with >=1 flagged array (any):", int(g[~g.called].flagged.gt(0).sum()))
# by order: S2 size>=2 and Hx
rows = []
for o in sorted(arr.order.unique()):
    s = unc[unc.order == o]
    rows.append(dict(order=o, S2=len(s), S2b_size_ge2=int((s.n >= 2).sum()), S2c_Hx=int(s.hx.sum()), S2d_both=int((s.hx & (s.n >= 2)).sum())))
T(pd.DataFrame(rows), "v_S2_by_order.tsv")
P("done")
