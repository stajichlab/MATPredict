#!/usr/bin/env python3
"""Review corrections (headline set H, LOCO).
1. ESM-2 150M plus CAAX refit with the bug fixed: the embedding block and the CAAX block are
   standardised separately and the CAAX block is weighted by sqrt(n_emb / n_caax) so the 6 CAAX
   columns are not drowned among 640 embedding columns (the first version multiplied by 3
   before a StandardScaler, which cancelled the weight).
2. Baseline 'CAAX within 10 kb OR flank gene within 50 kb' (no learning).
3. Bootstrap of dAUC vs the rule with genome, species and order as the resampling unit.
Reuses saved embeddings and out-of-fold scores; CPU only, seconds.
Usage: review_fixes.py  (run in the results directory)
"""
import os

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261006)
L = pd.read_csv(f"{HERE}/labelled_loci.tsv", sep="\t")
L["genome"] = L.asm.str[:15]
H = L[L.headline & L.label.isin(["mating", "other"])].reset_index(drop=True)
# species unit: strains of one species count once (5 R. toruloides, JEC20/JEC21 C. neoformans)
H["species"] = np.where(H.species.str.startswith(("Rhodotorula", "Cryptococcus")), H.species.str.split().str[:2].str.join(" "), H.species)
y = (H.label == "mating").astype(int).values
oof = pd.read_csv(f"{HERE}/oof_H_LOCO.tsv", sep="\t").set_index("locus_id").loc[H.locus_id]
ext = pd.read_csv(f"{HERE}/oof_H_extra_LOCO.tsv", sep="\t").set_index("locus_id").loc[H.locus_id]
emb = dict(zip([l.strip() for l in open(f"{HERE}/emb/emb.ids")], np.load(f"{HERE}/emb/emb.t30150M.npy")))
E = np.vstack([emb[i] for i in H.locus_id])
BIG = 1e7
prec = np.column_stack([
    (H.nT_10kb > 0).astype(float), np.log10(1 + H.d_T.clip(upper=BIG)), H.nT_20kb, H.nT_50kb,
    H.nHx_10kb_xo, np.log10(1 + H.d_Hx_xo.clip(upper=BIG))])
score = np.full(len(H), np.nan)
for o in H.order.unique():
    te = (H.order == o).values
    tr = ~te
    se, sp = StandardScaler().fit(E[tr]), StandardScaler().fit(prec[tr])
    w = np.sqrt(E.shape[1] / prec.shape[1])
    X = lambda m: np.hstack([se.transform(E[m]), w * sp.transform(prec[m])])
    m = LogisticRegression(C=0.01, class_weight="balanced", max_iter=5000).fit(X(tr), y[tr])
    score[te] = m.predict_proba(X(te))[:, 1]
or_rule = ((H.nT_10kb > 0) | (H.n_FLANK_50kb_xo > 0)).astype(float).values
S = {"rule_caax10": oof.rule_caax10.values, "CAAX_or_flank": or_rule, "ESM150_plus_CAAX_fixed": score,
     "ESM150_plus_CAAX_first_version": oof.ESMprec_LR_t30150M.values,
     "LR_flank": ext.LR_flank.values, "LR_all": oof.LR_all.values}


def sens95(yy, s):
    f, t, _ = roc_curve(yy, s)
    return t[f <= 0.05].max()


def units(col):
    u = H[col].values
    ug = np.unique(u)
    mem = {g: np.where(u == g)[0] for g in ug}
    return ug, mem


rows = []
for name, s in S.items():
    r = dict(model=name, n_pos=int(y.sum()), n_neg=int((1 - y).sum()), AUC=roc_auc_score(y, s), AP=average_precision_score(y, s),
             sens95=sens95(y, s))
    if name == "CAAX_or_flank":
        r["flagged_mating"] = int(s[y == 1].sum())
        r["flagged_other"] = int(s[y == 0].sum())
    for col in ("genome", "species", "order"):
        ug, mem = units(col)
        d, a = [], []
        for _ in range(2000):
            ix = np.concatenate([mem[g] for g in rng.choice(ug, len(ug))])
            if 0 < y[ix].sum() < len(ix):
                d.append(roc_auc_score(y[ix], s[ix]) - roc_auc_score(y[ix], S["rule_caax10"][ix]))
                a.append(roc_auc_score(y[ix], s[ix]))
        r[f"AUC_ci_{col}"] = f"{np.percentile(a, 2.5):.3f}-{np.percentile(a, 97.5):.3f}"
        r[f"dAUC_vs_rule_ci_{col}"] = f"{np.percentile(d, 2.5):.3f} to {np.percentile(d, 97.5):.3f}"
    rows.append(r)
out = pd.DataFrame(rows).round(3)
out.to_csv(f"{HERE}/review_fixes.tsv", sep="\t", index=False)
print(out.T.to_string())
print("species:", H.species.nunique(), "genomes:", H.genome.nunique(), "orders:", H.order.nunique())
print("mating loci with CAAX within 50kb:", int(((H.nT_50kb > 0) & (y == 1)).sum()), "within 10kb:", int(((H.nT_10kb > 0) & (y == 1)).sum()))
