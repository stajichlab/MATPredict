"""Tables for the note: python tables.py DIR   (after analyze.py and, optionally, the scan/ chunks)."""
import sys, json, math
from pathlib import Path
import pandas as pd

D = Path(sys.argv[1]); T = D / "tables"; T.mkdir(exist_ok=True)
L = pd.read_csv(D / "loci_audit_flagged.tsv", sep="\t", low_memory=False)
G = pd.read_csv(D / "genome_audit.tsv", sep="\t", low_memory=False)
M = pd.read_csv(D / "genome_list.tsv", sep="\t", dtype=str)
v06 = D / "v060" / "genomes_dothideomycetes_v060.tsv"
v06l = D / "v060" / "loci_dothideomycetes_v060.tsv"
rec = D / "v060" / "genomes_records_run_20261005.tsv"
for f in (v06, v06l, rec):
    if not f.exists():
        print("missing", f)

# ---- scan (genome-wide second-gene search), if present
scan = None
sc = sorted((D / "scan").glob("scan_*.tsv")) if (D / "scan").exists() else []
if sc:
    scan = pd.concat([pd.read_csv(f, sep="\t") for f in sc]).drop_duplicates("genome")
    scan.to_csv(D / "scan_all.tsv", sep="\t", index=False)
    G = G.merge(scan.add_prefix("gw_").rename(columns={"gw_genome": "genome"}), on="genome", how="left")

# ---- v0.6.0 and records-run
G["v060_idiomorphs"] = ""
if v06l.exists():
    l6 = pd.read_csv(v06l, sep="\t", dtype=str)
    G["v060_idiomorphs"] = G["genome"].map(l6[l6.idiomorph.isin(["MAT1-1", "MAT1-2"])].groupby("genome")["idiomorph"].agg(lambda s: "+".join(sorted(set(s)))))
    g6 = pd.read_csv(v06, sep="\t", dtype=str)
    G["in_v060"] = G["genome"].isin(set(g6["genome"]))
if rec.exists():
    r5 = pd.read_csv(rec, sep="\t", dtype=str).set_index("genome")
    G["records_run_idiomorphs"] = G["genome"].map(r5["idiomorphs"])

# ---- locus-level "both domains" test for homothallic-like loci
for k in ("pf_alpha_on_a", "pf_hmg_on_h"):
    L[k] = L[k].astype(str).eq("True")
L["hmg_score"] = pd.to_numeric(L["prof_hmg_on_hmg"], errors="coerce")
HMG_STRONG = 150.0   # bits on the HMG profile: lower end of single-HMG MAT1-2 loci in orders with a reference (156-250)
L["locus_both_domains"] = L["called"] & L["pf_alpha_on_a"] & L["pf_hmg_on_h"]
L["both_strong"] = L["locus_both_domains"] & ((L["idiomorph"] == "MAT1-2") | (L["hmg_score"] >= HMG_STRONG))
L["hmg_weak_beside_alpha"] = L["locus_both_domains"] & (L["idiomorph"] == "MAT1-1") & ~(L["hmg_score"] >= HMG_STRONG)
c = L[L["called"]].copy()
prim = (c.assign(rc=c["confidence"].map({"high": 0, "medium": 1, "low": 2}).fillna(3),
                 rl=c["locus_class"].map({"mat_locus": 0, "partial_locus": 1, "idiomorph_gene_only": 2}).fillna(3))
        .sort_values(["rc", "rl", "win_bits"], ascending=[True, True, False]).drop_duplicates("genome").set_index("genome"))
G["primary_both_domains"] = G["genome"].map(prim["locus_both_domains"])
G["any_both_domains"] = G["genome"].map(c.groupby("genome")["both_strong"].any())
G["any_hmg_weak"] = G["genome"].map(c.groupby("genome")["hmg_weak_beside_alpha"].any())
G["primary_core_lost_to_optional"] = G["genome"].map(prim["core_lost_to_optional"])
G["primary_own_id"] = G["genome"].map(prim["own_id"]); G["primary_own_cov"] = G["genome"].map(prim["own_cov"])
G["primary_win_gene"] = G["genome"].map(prim["win_gene"])
G["primary_win_optional"] = G["genome"].map(prim["win_optional"])
G["primary_flank"] = G["genome"].map(prim["flank_found"])


def verdict(r):
    if not isinstance(r["primary_idiomorph"], str):
        return "no_call"
    both_idio = isinstance(r["idiomorphs_called"], str) and "+" in r["idiomorphs_called"]
    bad, soft = [], []
    if r["primary_vs_profile"] == "disagree":
        bad.append("profile_disagrees")
    if str(r["any_both_domains"]) == "True":
        bad.append("both_genes_strong")
    if both_idio:
        bad.append("two_idiomorph_loci")
    if str(r["primary_core_lost_to_optional"]) == "True":
        bad.append("core_set_aside_for_optional")
    if bad:
        return "conflicting:" + "+".join(bad)
    if r["primary_domain"] == "no":
        soft.append("domain_absent")
    if str(r["any_hmg_weak"]) == "True":
        soft.append("hmg_beside_alpha_below_MAT1-2-1_range")
    if r["primary_vs_profile"] in ("profile_undetermined", "profile_not_scored"):
        soft.append("profile_" + r["primary_vs_profile"].split("_", 1)[1])
    if soft:
        return "unverifiable:" + "+".join(soft)
    if r["primary_vs_profile"] == "agree" and r["primary_domain"] == "yes":
        return "confirmed"
    return "unverifiable:other"


G["verdict"] = G.apply(verdict, axis=1)
G["verdict_class"] = G["verdict"].str.split(":").str[0]
G.to_csv(D / "genome_audit_verdict.tsv", sep="\t", index=False)
L.to_csv(D / "loci_audit_flagged.tsv", sep="\t", index=False)

# ---- counts
def counts(df, by):
    g = df.groupby(by)
    t = pd.DataFrame({
        "genomes": g.size(),
        "report": g.apply(lambda x: (x["status"] != "no_report").sum()),
        "called_single": g.apply(lambda x: x["primary_idiomorph"].notna().sum()),
        "MAT1-1": g.apply(lambda x: (x["primary_idiomorph"] == "MAT1-1").sum()),
        "MAT1-2": g.apply(lambda x: (x["primary_idiomorph"] == "MAT1-2").sum()),
        "both_idio_loci": g.apply(lambda x: x["idiomorphs_called"].astype(str).str.contains(r"\+").sum()),
        "undetermined_only": g.apply(lambda x: ((x["n_detected"] > 0) & x["primary_idiomorph"].isna()).sum()),
        "not_detected": g.apply(lambda x: (x["status"] == "not_detected").sum()),
        "v060_called": g.apply(lambda x: (x["v060_idiomorphs"].astype(str).str.len() > 0).sum() if "v060_idiomorphs" in x else 0),
    })
    for k in ("confirmed", "conflicting", "unverifiable"):
        t[k] = g.apply(lambda x, k=k: (x["verdict_class"] == k).sum())
    return t


G["order"] = G["order"].fillna("(blank)")
counts(G, "order").sort_values("genomes", ascending=False).to_csv(T / "by_order.tsv", sep="\t")
G["og"] = G["order"] + " / " + G["genus"].fillna("?")
counts(G, ["order", "genus"]).sort_values("genomes", ascending=False).to_csv(T / "by_genus.tsv", sep="\t")
counts(G, ["order", "family"]).sort_values("genomes", ascending=False).to_csv(T / "by_family.tsv", sep="\t")
print(open(T / "by_order.tsv").read())

# call classes
c = L[L["called"]]
pd.crosstab(c["order"], [c["locus_class"]]).to_csv(T / "locus_class_by_order.tsv", sep="\t")
pd.crosstab(c["order"], [c["confidence"]]).to_csv(T / "confidence_by_order.tsv", sep="\t")
pd.crosstab([c["order"], c["idiomorph"]], c["vote_vs_profile"]).to_csv(T / "vote_vs_profile_by_order.tsv", sep="\t")
pd.crosstab([c["genus"], c["idiomorph"]], c["vote_vs_profile"]).to_csv(T / "vote_vs_profile_by_genus.tsv", sep="\t")
pd.crosstab([c["order"], c["idiomorph"]], c["domain_support"]).to_csv(T / "domain_by_order.tsv", sep="\t")
pd.crosstab(G["order"], G["verdict_class"]).to_csv(T / "verdict_by_order.tsv", sep="\t")

# species ratios
sp = G[G["primary_idiomorph"].notna()].groupby("species")["primary_idiomorph"].agg(
    n="size", m1=lambda s: (s == "MAT1-1").sum(), m2=lambda s: (s == "MAT1-2").sum())
def binom_p(k, n):
    # two-sided exact binomial p at 0.5
    if n == 0:
        return float("nan")
    pk = [math.comb(n, i) * 0.5 ** n for i in range(n + 1)]
    return min(1.0, sum(p for p in pk if p <= pk[k] + 1e-12))
sp["p_binom_1to1"] = [binom_p(int(r.m1), int(r.n)) for r in sp.itertuples()]
sp.sort_values("n", ascending=False).to_csv(T / "species_ratio.tsv", sep="\t")
