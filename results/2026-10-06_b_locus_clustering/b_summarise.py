#!/usr/bin/env python3
"""Part 2/3: B array vs other arrays in the 6 genomes with curated B records; cassette statistics over the 1,287 genomes."""
import os, glob, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
L = pd.read_csv(HERE + "/b_loci_cassette.tsv.gz", sep="\t", low_memory=False)
A = pd.read_csv(HERE + "/b_arrays_features.tsv.gz", sep="\t")
lab = pd.read_csv(HERE + "/labelled_loci_from_pr-receptor-ml-investigation.tsv", sep="\t")
C = pd.read_csv(HERE + "/curated_agaricomycete_B_records.tsv", sep="\t")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
agar = lab[lab.order.isin(["Agaricales", "Polyporales", "Russulales"]) & (lab.label == "mating")]
# map labelled mating loci to inventory loci by overlap
Bm = []
for _, r in agar.iterrows():
    x = L[(L.genome == r.asm) & (L.contig == r.contig) & (L.start < r.end) & (L.end > r.start)]
    for _, y in x.iterrows(): Bm.append((r.species, r.grade, y.array_id, y.locus_id))
Bm = pd.DataFrame(Bm, columns=["species", "grade", "array_id", "locus_id"]).drop_duplicates()
Bm.to_csv(HERE + "/b_labelled_loci_mapped.tsv", sep="\t", index=False)
bar = Bm.groupby("array_id").agg(species=("species", "first"), grade=("grade", "first"), n_B_loci=("locus_id", "size")).reset_index()
G = A[A.genome.isin(L[L.locus_id.isin(Bm.locus_id)].genome.unique())].merge(bar[["array_id", "n_B_loci"]], on="array_id", how="left")
G["is_B"] = G.n_B_loci.notna()
# held-out precursor homology (curated own-species records excluded) inside array +-10 kb
def hx_ho(r):
    f = HERE + f"/hx_heldout/{r.genome}.hx_heldout.tsv"
    h = pd.read_csv(f, sep="\t"); h = h[h.contig == r.contig]
    p = np.sort(h.pos.values); p = p[(p >= r.start - 10000) & (p <= r.end + 10000)]
    if len(p) == 0: return 0
    return 1 + int((np.diff(p) > 300).sum())
G["aH10_heldout"] = G.apply(hx_ho, axis=1)
G["n_loci_in_array_nonB"] = G.n_loci - G.n_B_loci.fillna(0)
cols = ["species", "array_id", "is_B", "n_loci", "n_B_loci", "span", "aT10", "aH10", "aH10_heldout", "aTH10", "any_casA", "d_STE20", "in_call", "withheld"]
G = G.sort_values(["species", "is_B", "n_loci"], ascending=[True, False, False])
G[cols].to_csv(HERE + "/b_vs_other_arrays_6genomes.tsv", sep="\t", index=False)
print(G[cols].to_string())
print("\nB arrays n=", G.is_B.sum(), " other arrays n=", (~G.is_B).sum(), " other arrays with >=2 loci:", ((~G.is_B) & (G.n_loci >= 2)).sum())
# curated receptors vs inventory loci
for _, c in C.iterrows():
    pass
# ---- full set
A["order2"] = A.order.where(A.order.isin(["Agaricales", "Boletales", "Polyporales", "Cantharellales", "Russulales", "Hymenochaetales", "Auriculariales"]), "other")
print("\nFULL SET arrays", len(A), "genomes", A.genome.nunique())
print("loci with cassette A (>=2 strict CAAX ORFs within 5 kb):", L.casA.sum(), "of", len(L), " expected by contig-uniform null:", round(L.null_casA.sum(), 1))
print("casB (>=2 precursor cands incl >=1 CAAX+homology):", L.casB.sum(), " casC (>=2 CAAX ORFs each with homology):", L.casC.sum())
# per order per class
res = []
for o, g in L.groupby("order"):
    res.append(dict(order=o, loci=len(g), genomes=g.genome.nunique(), casA=g.casA.sum(), casA_exp=round(g.null_casA.sum(), 1), casB=g.casB.sum(), casC=g.casC.sum(),
                    genomes_casA=g[g.casA].genome.nunique(), genomes_casB=g[g.casB].genome.nunique(), genomes_casC=g[g.casC].genome.nunique(),
                    called=g.pipeline_call.sum(), casA_called=(g.casA & g.pipeline_call).sum(), casA_notcalled=(g.casA & ~g.pipeline_call).sum(),
                    casB_notcalled=(g.casB & ~g.pipeline_call).sum(), casC_notcalled=(g.casC & ~g.pipeline_call).sum(),
                    casC_notcalled_genomes_without_call=0))
R = pd.DataFrame(res); R.loc[len(R)] = ["ALL"] + [R[c].sum() if R[c].dtype != object else "" for c in R.columns[1:]]
gc = L.groupby("genome").pipeline_call.any()
R.to_csv(HERE + "/cassette_by_order.tsv", sep="\t", index=False); print(R.to_string())
# genomes recoverable: no call, but cassette
for k in ["casA", "casB", "casC"]:
    gk = set(L[L[k]].genome); nocall = set(gc[~gc].index)
    gk_nc = set(L[L[k] & ~L.pipeline_call].genome)
    print(k, "genomes with a cassette locus", len(gk), "; with cassette but no PR-called locus anywhere", len(gk & nocall), "; with an uncalled cassette locus", len(gk_nc))
# array-level: number of loci in cassette per array, size classes
A["cas_cls"] = np.where(A.any_casC, "C", np.where(A.any_casB, "B", np.where(A.any_casA, "A", "none")))
tab = pd.crosstab(A.n_loci.clip(upper=5), A.cas_cls); print(tab); tab.to_csv(HERE + "/cassette_by_array_size.tsv", sep="\t")
# per-genome number of cassette arrays (casA)
pg = A.groupby("genome").any_casA.sum(); print("cassetteA arrays per genome: dist", pg.value_counts().sort_index().to_dict())
pg2 = A.groupby("genome").any_casC.sum(); print("cassetteC arrays per genome: dist", pg2.value_counts().sort_index().to_dict())
# STE20 near as semi-independent check by cassette status among arrays on STE20-bearing genomes
for k in ["none", "A", "B", "C"]:
    x = A[(A.cas_cls == k) & A.has_STE20.fillna(False)]
    print("STE20 within 250 kb, cassette class", k, len(x), round(100 * x.ste20_250.mean(), 1))
for lab_, m in [("size1 none", (A.n_loci == 1) & (A.cas_cls == "none")), ("size1 A+", (A.n_loci == 1) & (A.cas_cls != "none")),
                ("size>=2 none", (A.n_loci >= 2) & (A.cas_cls == "none")), ("size>=2 A+", (A.n_loci >= 2) & (A.cas_cls != "none")),
                ("size>=2 C", (A.n_loci >= 2) & (A.cas_cls == "C")), ("size>=3 A+", (A.n_loci >= 3) & (A.cas_cls != "none"))]:
    x = A[m & A.has_STE20.fillna(False)]; print("STE20<=250kb", lab_, len(x), round(100 * x.ste20_250.mean(), 1))
# cassette not called, by order and size
nc = A[(A.cas_cls != "none") & ~A.in_call]
print("arrays with cassette, no PR call in array:", len(nc), nc.groupby("order2").size().to_dict())
print("  with >=2 loci:", (nc.n_loci >= 2).sum(), " class C:", (nc.cas_cls == "C").sum(), " C & >=2 loci:", ((nc.cas_cls == "C") & (nc.n_loci >= 2)).sum())
nc.to_csv(HERE + "/cassette_arrays_without_call.tsv.gz", sep="\t", index=False)
# arrays multiple cassette-bearing loci (size>=2 with >=2 casA loci)
print("arrays with >=2 casA loci:", (A.n_casA >= 2).sum(), " of arrays >=2:", (A.n_loci >= 2).sum())
print("unique genomes: arrays>=3 loci", A[A.n_loci >= 3].genome.nunique(), "; with casA and >=3", A[(A.n_loci >= 3) & A.any_casA].genome.nunique())
# array-property percentile of B arrays in the full set
for _, b in G[G.is_B].iterrows():
    print(b.species, "size", b.n_loci, "span", b.span, "arrays with size>= :", (A.n_loci >= b.n_loci).sum(), "span>= :", (A.span >= b.span).sum(), "aT10>= :", (A.aT10 >= b.aT10).sum(), "(>=2 loci & aT10>=): ", ((A.n_loci >= b.n_loci) & (A.aT10 >= b.aT10)).sum())
