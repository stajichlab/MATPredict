"""Join the genome-wide alpha/HMG ORF scan (scan/ chunks) to the primary call of each genome.  python scan_integrate.py DIR
Writes tables/genomewide_second_gene.tsv, genomewide_by_genus.tsv, genomewide_by_order.tsv and adds gw_* columns to
genome_audit_verdict.tsv.  Independent of the pipeline: six-frame ORFs scored with Pfam and the curated profiles."""
import sys
from pathlib import Path
import pandas as pd

D = Path(sys.argv[1]); T = D / "tables"
sc = sorted((D / "scan").glob("scan_*.tsv"))
scan = pd.concat([pd.read_csv(f, sep="\t") for f in sc]).drop_duplicates("genome")
scan.to_csv(D / "scan_all.tsv", sep="\t", index=False)
L = pd.read_csv(D / "loci_audit_flagged.tsv", sep="\t", low_memory=False)
G = pd.read_csv(D / "genome_audit_verdict.tsv", sep="\t", low_memory=False)
G = G[[c for c in G.columns if not c.startswith("gw_")]]
c = L[L["called"]]
prim = (c.assign(rc=c["confidence"].map({"high": 0, "medium": 1, "low": 2}).fillna(3),
                 rl=c["locus_class"].map({"mat_locus": 0, "partial_locus": 1, "idiomorph_gene_only": 2}).fillna(3))
        .sort_values(["rc", "rl", "win_bits"], ascending=[True, True, False]).drop_duplicates("genome").set_index("genome"))
S = scan.set_index("genome")


def parse_hits(s):
    out = []
    if not isinstance(s, str) or not s:
        return out
    for h in s.split(";"):
        try:
            loc, sc_, ln = h.rsplit(":", 2)
            contig, rng = loc.rsplit(":", 1)
            a, b = rng[:-1].split("-")
            out.append((contig, int(a), int(b), float(sc_), int(ln)))
        except Exception:  # noqa: BLE001
            continue
    return out


def at_locus(h, contig, start, end, pad=20000):
    return h[0] == contig and h[2] >= start - pad and h[1] <= end + pad


rows = []
for g in G["genome"]:
    row = dict(genome=g, gw_scanned=g in S.index, gw_class=None)
    if g in S.index:
        s = S.loc[g]
        a_hits, h_hits = parse_hits(s["alpha_hits40"]), parse_hits(s["hmg_hits40"])
        row.update(gw_alpha_best=s["alpha_prof_best"], gw_hmg_best=s["hmg_prof_best"],
                   gw_alpha_n100=sum(1 for h in a_hits if h[3] >= 100), gw_hmg_n150=sum(1 for h in h_hits if h[3] >= 150),
                   gw_alpha_pfam_best=s["alpha_pfam_best"])
        if g in prim.index:
            p = prim.loc[g]
            c0, s0, e0, idio = p["contig"], int(p["start"]), int(p["end"]), p["idiomorph"]
            a_loc = [h for h in a_hits if at_locus(h, c0, s0, e0)]
            h_loc = [h for h in h_hits if at_locus(h, c0, s0, e0)]
            row["gw_alpha_best_at_locus"] = max([h[3] for h in a_loc], default=None)
            row["gw_hmg_best_at_locus"] = max([h[3] for h in h_loc], default=None)
            a_else = [h for h in a_hits if h[3] >= 100 and not at_locus(h, c0, s0, e0)]
            h_else = [h for h in h_hits if h[3] >= 150 and not at_locus(h, c0, s0, e0)]
            # Baseline: in heterothallic genera (Parastagonospora, Alternaria, Corynespora, Ascochyta, Fulvia, Hortaea) the scan
            # finds no HMG ORF >= 40 bits within 20 kb of a MAT1-1 call and no alpha ORF >= 40 bits at a MAT1-2 call.
            if h_loc:
                bh = max(h_loc, key=lambda h: h[3]); row["gw_hmg_orf_len_at_locus"] = bh[4]
            if a_loc:
                ba = max(a_loc, key=lambda h: h[3]); row["gw_alpha_orf_len_at_locus"] = ba[4]
            if idio == "MAT1-1":
                row["gw_class"] = ("MAT1-1: HMG ORF at the locus" if h_loc else
                                   "MAT1-1: strong HMG ORF elsewhere" if h_else else
                                   "MAT1-1: no HMG ORF at or beyond the locus")
            else:
                row["gw_class"] = ("MAT1-2: alpha-box ORF at the locus" if a_loc else
                                   "MAT1-2: alpha-box ORF elsewhere" if a_else else
                                   "MAT1-2: no alpha-box ORF in the genome")
    rows.append(row)
GW = pd.DataFrame(rows)
GW.to_csv(T / "genomewide_second_gene.tsv", sep="\t", index=False)
G2 = G.merge(GW, on="genome", how="left")
G2.to_csv(D / "genome_audit_verdict.tsv", sep="\t", index=False)
x = G2[G2["gw_class"].notna()]
pd.crosstab([x["order"], x["genus"]], x["gw_class"]).to_csv(T / "genomewide_by_genus.tsv", sep="\t")
pd.crosstab(x["order"], x["gw_class"]).to_csv(T / "genomewide_by_order.tsv", sep="\t")
print(len(scan), "scanned;", len(x), "with a primary call")
print(x["gw_class"].value_counts())
