"""Aggregate loci_audit.tsv / genome_status.tsv (audit_extract.py) with taxonomy, v0.6.0 calls, BioSample truth and the
curated-record genomes into the tables of analysis/2026-10-18_dothideomycetes-mat-audit.md.   python analyze.py DIR
Writes tables/*.tsv under DIR.  Pure aggregation; no thresholds other than the ones named in FLAGS below."""
import sys, json, re
from pathlib import Path
import pandas as pd

D = Path(sys.argv[1])
T = D / "tables"; T.mkdir(exist_ok=True)
L = pd.read_csv(D / "loci_audit.tsv", sep="\t", dtype={"genome": str})
G = pd.read_csv(D / "genome_status.tsv", sep="\t", dtype={"genome": str})
M = pd.read_csv(D / "genome_list.tsv", sep="\t", dtype=str)
rk = D / "taxid_ranks.tsv"
fam = {}
if rk.exists():
    for line in open(rk):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 5:
            fam[p[0]] = (p[2], p[3], p[4])
M["family"] = M["NCBI_TAXID"].map(lambda t: fam.get(t, ("", "", ""))[1])
M["order_tk"] = M["NCBI_TAXID"].map(lambda t: fam.get(t, ("", "", ""))[0])
M["order_final"] = M["ORDER"].where(M["ORDER"].fillna("") != "", M["order_tk"])
meta = M.set_index("ASM_ACCESSION")
for df in (L, G):
    df["family"] = df["genome"].map(meta["family"])
    df["order"] = df["genome"].map(meta["order_final"]).fillna(df.get("order"))
    df["genus"] = df["genome"].map(meta["GENUS"])
    df["species"] = df["genome"].map(meta["SPECIES"])

# ---------------------------------------------------------------- flags (per locus)
single = L["idiomorph"].isin(["MAT1-1", "MAT1-2"])
L["called"] = single
L["prof_agree"] = L["prof_call"].map(lambda x: x if x in ("MAT1-1", "MAT1-2") else "undet_or_none")
L["vote_vs_profile"] = "n/a"
c = L["called"]
L.loc[c & (L["prof_call"] == L["idiomorph"]), "vote_vs_profile"] = "agree"
L.loc[c & L["prof_call"].isin(["MAT1-1", "MAT1-2"]) & (L["prof_call"] != L["idiomorph"]), "vote_vs_profile"] = "disagree"
L.loc[c & (L["prof_call"] == "undetermined"), "vote_vs_profile"] = "profile_undetermined"
L.loc[c & L["prof_call"].isin(["no_protein"]) | (c & L["prof_call"].isna()), "vote_vs_profile"] = "profile_not_scored"
for k in ("own_cov", "own_id", "own_bits", "oth_cov", "own_edge", "win_bits"):
    L[k] = pd.to_numeric(L[k], errors="coerce")
L["own_core_gene"] = L["idiomorph"].map({"MAT1-1": "MAT1-1-1", "MAT1-2": "MAT1-2-1"})
L["F_optional_wins"] = c & L["win_optional"].astype(str).eq("True")
L["F_own_core_missing"] = c & ~L["own_present"].astype(str).eq("True")
L["F_own_partial_edge"] = c & L["own_present"].astype(str).eq("True") & (L["own_cov"] < 80) & (L["own_edge"] < 500)
L["F_own_partial_interior"] = c & L["own_present"].astype(str).eq("True") & (L["own_cov"] < 80) & ~(L["own_edge"] < 500)
L["F_low_identity"] = c & (L["own_id"] < 40)
L["F_profile_disagree"] = L["vote_vs_profile"].eq("disagree")


BITS = {"MAT1-1-1": "a_bits", "MAT1-2-1": "h_bits"}


def pf(row, gene, names):
    """True/False when the gene's protein was scored (an empty Pfam string means scored, no domain above GA);
    None when the gene is not at the locus or its protein could not be built."""
    s = row.get("pfam_" + gene)
    if isinstance(s, str):
        return any(n in s for n in names)
    present = pd.notna(row.get(BITS[gene])) and pd.notna(row.get("prof_call"))
    return False if present else None


ALPHA = ("MATalpha_HMGbox",)
HMG = ("HMG_box", "HMG_box_2", "HMG_box_5", "HMG-box")
L["pf_alpha_on_a"] = L.apply(lambda r: pf(r, "MAT1-1-1", ALPHA), axis=1)
L["pf_hmg_on_h"] = L.apply(lambda r: pf(r, "MAT1-2-1", HMG), axis=1)
L["pf_alpha_on_h"] = L.apply(lambda r: pf(r, "MAT1-2-1", ALPHA), axis=1)   # alpha box in the HMG gene: paralog / mislabelled
L["pf_hmg_on_a"] = L.apply(lambda r: pf(r, "MAT1-1-1", HMG), axis=1)
# domain support for the called idiomorph: MAT1-1 needs an alpha-box protein at the locus; MAT1-2 needs an HMG protein
# at the locus AND no alpha-box protein (an HMG protein beside MAT1-1-1 is the MAT1-1-3 pattern, not MAT1-2-1)
def dom_support(r):
    if r["idiomorph"] == "MAT1-1":
        v = r["pf_alpha_on_a"]
        return "yes" if v is True else "no" if v is False else "not_scored"
    if r["idiomorph"] == "MAT1-2":
        v = r["pf_hmg_on_h"]
        if v is None or pd.isna(v):
            return "not_scored"
        if v is False:
            return "no"
        return "yes_but_alpha_gene_also_present" if r["pf_alpha_on_a"] is True else "yes"
    return "n/a"
L["domain_support"] = L.apply(dom_support, axis=1)
L.to_csv(D / "loci_audit_flagged.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- genome-level
called = L[L["called"]]
rank_c = {"high": 0, "medium": 1, "low": 2}
rank_l = {"mat_locus": 0, "partial_locus": 1, "idiomorph_gene_only": 2}
prim = (called.assign(rc=called["confidence"].map(rank_c).fillna(3), rl=called["locus_class"].map(rank_l).fillna(3))
        .sort_values(["rc", "rl", "win_bits"], ascending=[True, True, False])
        .drop_duplicates("genome").set_index("genome"))
G["n_called_loci"] = G["genome"].map(called.groupby("genome").size()).fillna(0).astype(int)
G["idiomorphs_called"] = G["genome"].map(called.groupby("genome")["idiomorph"].agg(lambda s: "+".join(sorted(set(s)))))
G["primary_idiomorph"] = G["genome"].map(prim["idiomorph"])
G["primary_confidence"] = G["genome"].map(prim["confidence"])
G["primary_class"] = G["genome"].map(prim["locus_class"])
G["primary_profile"] = G["genome"].map(prim["prof_call"])
G["primary_vs_profile"] = G["genome"].map(prim["vote_vs_profile"])
G["primary_domain"] = G["genome"].map(prim["domain_support"])
G.to_csv(D / "genome_audit.tsv", sep="\t", index=False)
print("loci", len(L), "called", int(L["called"].sum()), "genomes", len(G), "with called locus", int((G["n_called_loci"] > 0).sum()))
