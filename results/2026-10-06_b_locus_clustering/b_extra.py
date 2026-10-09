import os, numpy as np, pandas as pd
from scipy.stats import fisher_exact
HERE = os.path.dirname(os.path.abspath(__file__))
L = pd.read_csv(HERE + "/b_loci_cassette.tsv.gz", sep="\t", low_memory=False); A = pd.read_csv(HERE + "/b_arrays_features.tsv.gz", sep="\t")
Bm = pd.read_csv(HERE + "/b_labelled_loci_mapped.tsv", sep="\t"); print(Bm.groupby("species").size())
print("null casA among uncalled loci:", round(L[~L.pipeline_call].null_casA.sum(), 1), "obs", int((L.casA & ~L.pipeline_call).sum()), "; among called:", round(L[L.pipeline_call].null_casA.sum(), 1), "obs", int((L.casA & L.pipeline_call).sum()))
print("null casA uncalled, size>=2 arrays:")
n = L[L.array_size >= 2]; print(" obs", int((n.casA & ~n.pipeline_call).sum()), "exp", round(n[~n.pipeline_call].null_casA.sum(), 1))
A["cas"] = A.any_casA
S = A[A.has_STE20.fillna(False)]
for nm, m in [("size>=2 flagged, not cassette", (S.n_loci >= 2) & S.caax_flag & ~S.cas), ("size>=2 flagged, cassette", (S.n_loci >= 2) & S.caax_flag & S.cas), ("size>=2 unflagged", (S.n_loci >= 2) & ~S.caax_flag),
              ("size1 flagged not cassette", (S.n_loci == 1) & S.caax_flag & ~S.cas), ("size1 cassette", (S.n_loci == 1) & S.cas)]:
    x = S[m]; print("STE20<=250kb", nm, len(x), round(100 * x.ste20_250.mean(), 1))
nc = A[A.cas & ~A.in_call]; hascall = A.groupby("genome").in_call.any()
for nm, m in [("cassette A uncalled", nc), ("A uncalled >=2 loci", nc[nc.n_loci >= 2]), ("C uncalled", nc[nc.any_casC]), ("C or (A and >=2 loci)", nc[nc.any_casC | (nc.n_loci >= 2)])]:
    g = set(m.genome); print(nm, "arrays", len(m), "genomes", len(g), "genomes with no call at all", sum(not hascall[x] for x in g))
G = pd.read_csv(HERE + "/b_vs_other_arrays_6genomes.tsv", sep="\t")
t = [[G[G.is_B].any_casA.sum(), (~G[G.is_B].any_casA).sum()], [G[~G.is_B].any_casA.sum(), (~G[~G.is_B].any_casA).sum()]]
print("cassetteA B vs other (6 genomes)", t, fisher_exact(t))
ind = G[G.species.isin(["Schizophyllum commune", "Coprinopsis cinerea"])]
t = [[ind[ind.is_B].any_casA.sum(), (~ind[ind.is_B].any_casA).sum()], [ind[~ind.is_B].any_casA.sum(), (~ind[~ind.is_B].any_casA).sum()]]
print("independent-label genomes only", t, fisher_exact(t))
# cassette array per genome by order: genomes with >=1 and exactly 1
pg = A.groupby(["order", "genome"]).cas.sum().reset_index()
print(pg.groupby("order").apply(lambda d: pd.Series(dict(genomes=len(d), ge1=(d.cas >= 1).sum(), exactly1=(d.cas == 1).sum()))).query("ge1>0"))
