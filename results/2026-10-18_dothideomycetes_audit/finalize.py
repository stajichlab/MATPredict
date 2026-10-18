"""Final per-genome verdict (adds the genome-wide scan) and group classification.  python finalize.py DIR
Inputs: genome_audit_verdict.tsv (tables.py + scan_integrate.py), literature.tsv.  Outputs: genome_audit_final.tsv,
tables/group_genus.tsv, group_family.tsv, group_order.tsv, group_literature_check.tsv.

Verdict per genome (primary locus = highest confidence, then locus class, then bitscore):
  no_call      no single-idiomorph call
  conflicting  the profile disagrees with the vote, OR a strong second-gene signal sits at the locus (a homothallic-like
               locus called as one idiomorph), OR two loci with different idiomorphs, OR a core gene was set aside for an
               optional gene
  unverifiable weak second-gene signal (HMG ORF beside MAT1-1-1 below the MAT1-2-1 range), domain absent on the called gene,
               profile undetermined or not scored
  confirmed    profile and Pfam domain agree with the vote and nothing else flags
Group label (genus / family / order), over genomes with a call: confirmed if >= 90% confirmed and none conflicting;
conflicting if >= 20% conflicting; unverifiable if < 50% confirmed and < 20% conflicting; mixed otherwise."""
import sys
from pathlib import Path
import pandas as pd

D = Path(sys.argv[1]); T = D / "tables"
G = pd.read_csv(D / "genome_audit_verdict.tsv", sep="\t", low_memory=False)
have_scan = "gw_class" in G.columns and G["gw_class"].notna().any()
G["second_gene_scan"] = "not_scanned"
if have_scan:
    def sg(r):
        cl = r.get("gw_class")
        if not isinstance(cl, str):
            return "not_scanned"
        if cl.startswith("MAT1-1: HMG ORF at the locus"):
            return "strong" if float(r.get("gw_hmg_best_at_locus") or 0) >= 90 else "weak"
        if cl.startswith("MAT1-2: alpha-box ORF at the locus"):
            return "strong" if float(r.get("gw_alpha_best_at_locus") or 0) >= 90 else "weak"
        return "none"
    G["second_gene_scan"] = G.apply(sg, axis=1)

# domain evidence for the called gene: Pfam on the pipeline's gene model; where that model is a fragment with no Pfam hit, an
# annotation-free scan ORF at the locus (own-idiomorph profile >= 60 bits, within 20 kb) counts as the domain being there
def dom_final(r):
    d = r["primary_domain"]
    if d == "yes" or not have_scan or not isinstance(r.get("gw_class"), str):
        return d
    own = r.get("gw_alpha_best_at_locus") if r["primary_idiomorph"] == "MAT1-1" else r.get("gw_hmg_best_at_locus")
    own = pd.to_numeric(own, errors="coerce")
    if d in ("no", "not_scored") and pd.notna(own) and own >= 60:
        return "yes_by_scan"
    return d


G["domain_final"] = G.apply(dom_final, axis=1)
G["second_gene_pfam"] = "none"
G.loc[G["any_hmg_weak"].astype(str) == "True", "second_gene_pfam"] = "weak"
G.loc[G["any_both_domains"].astype(str) == "True", "second_gene_pfam"] = "strong"
rank = {"none": 0, "not_scanned": 0, "weak": 1, "strong": 2}
G["second_gene"] = [a if rank[a] >= rank[b] else b for a, b in zip(G["second_gene_pfam"], G["second_gene_scan"])]


def verdict(r):
    if not isinstance(r["primary_idiomorph"], str):
        return "no_call"
    both_idio = isinstance(r["idiomorphs_called"], str) and "+" in r["idiomorphs_called"]
    bad, soft = [], []
    if r["primary_vs_profile"] == "disagree":
        bad.append("profile_disagrees")
    if r["second_gene"] == "strong":
        bad.append("second_gene_at_locus")
    if both_idio:
        bad.append("two_idiomorph_loci")
    if str(r["primary_core_lost_to_optional"]) == "True":
        bad.append("core_set_aside_for_optional")
    if bad:
        return "conflicting:" + "+".join(bad)
    if r["second_gene"] == "weak":
        soft.append("weak_second_gene_at_locus")
    if r["domain_final"] == "no":
        soft.append("domain_absent")
    if r["primary_vs_profile"] in ("profile_undetermined", "profile_not_scored"):
        soft.append(r["primary_vs_profile"])
    if soft:
        return "unverifiable:" + "+".join(soft)
    if r["primary_vs_profile"] == "agree" and r["domain_final"] in ("yes", "yes_by_scan"):
        return "confirmed"
    return "unverifiable:other"


G["verdict"] = G.apply(verdict, axis=1)
G["verdict_class"] = G["verdict"].str.split(":").str[0]
G.to_csv(D / "genome_audit_final.tsv", sep="\t", index=False)


def label(row):
    n = row["called"]
    if n == 0:
        return "no_call"
    fc, cc = row["confirmed"] / n, row["conflicting"] / n
    if cc >= 0.20:
        return "conflicting"
    if fc >= 0.90 and row["conflicting"] == 0:
        return "confirmed"
    if fc < 0.50:
        return "unverifiable"
    return "mixed"


def group(by, name):
    g = G.groupby(by)
    t = pd.DataFrame({
        "genomes": g.size(),
        "called": g.apply(lambda x: (x["verdict_class"] != "no_call").sum()),
        "MAT1-1": g.apply(lambda x: (x["primary_idiomorph"] == "MAT1-1").sum()),
        "MAT1-2": g.apply(lambda x: (x["primary_idiomorph"] == "MAT1-2").sum()),
        "confirmed": g.apply(lambda x: (x["verdict_class"] == "confirmed").sum()),
        "conflicting": g.apply(lambda x: (x["verdict_class"] == "conflicting").sum()),
        "unverifiable": g.apply(lambda x: (x["verdict_class"] == "unverifiable").sum()),
    })
    # most common reason among the non-confirmed
    def reasons(x):
        v = x.loc[x["verdict_class"].isin(["conflicting", "unverifiable"]), "verdict"]
        return "; ".join(f"{k} ({n})" for k, n in v.value_counts().head(3).items())
    t["main_reasons"] = g.apply(reasons)
    t["label"] = t.apply(label, axis=1)
    t["call_rate"] = (t["called"] / t["genomes"]).round(2)
    t.loc[(t["called"] > 0) & (t["call_rate"] < 0.5), "label"] += " (low call rate)"
    t.sort_values("genomes", ascending=False).to_csv(T / f"group_{name}.tsv", sep="\t")
    return t


G["family"] = G["family"].fillna("(blank)")
go = group("order", "order"); gg = group(["order", "genus"], "genus"); gf = group(["order", "family"], "family")
print(go[["genomes", "called", "confirmed", "conflicting", "unverifiable", "label"]].sort_values("genomes", ascending=False).head(12))
print(gg["label"].value_counts()); print(G["verdict_class"].value_counts())
print("genomes in groups labelled conflicting:", int(gg.loc[gg["label"] == "conflicting", "genomes"].sum()))
