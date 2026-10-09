#!/usr/bin/env python3
"""Learning curve: AUC on held-out genomes as a function of the number of training genomes
(headline set H; handcrafted logistic regression, ESM-2 150M logistic regression, and the CAAX
rule as reference). Training genomes are drawn at random; the test set is every genome not drawn.
Writes learning_curve.tsv."""
import importlib.util
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("an", os.path.join(HERE, "analysis.py"))
an = importlib.util.module_from_spec(spec)
spec.loader.exec_module(an)
H = an.L[an.L.headline].reset_index(drop=True)
y = H.y.values
g = H.genome.values
ug = np.unique(g)
rng = np.random.default_rng(7)
Xh = H[an.FEAT_PREC + an.FEAT_MAT + an.FEAT_PROT].values.astype(float)
Xe = np.vstack([an.EMB["t30150M"][i] for i in H.locus_id])
rows = []
for k in (2, 4, 6, 8, 10, 12):
    res = {"hand": [], "esm": [], "rule": []}
    for _ in range(150):
        tr_g = rng.choice(ug, k, replace=False)
        tr = np.isin(g, tr_g)
        te = ~tr
        if y[tr].sum() < 2 or (1 - y[tr]).sum() < 2 or y[te].sum() < 2 or (1 - y[te]).sum() < 2:
            continue
        m = an.lr(0.3).fit(Xh[tr], y[tr])
        res["hand"].append(an.auc(y[te], m.predict_proba(Xh[te])[:, 1]))
        m = an.lr(0.01).fit(Xe[tr], y[tr])
        res["esm"].append(an.auc(y[te], m.predict_proba(Xe[te])[:, 1]))
        res["rule"].append(an.auc(y[te], H.caax10.values[te].astype(float)))
    for name, v in res.items():
        rows.append(dict(train_genomes=k, model=name, reps=len(v), AUC_mean=round(np.mean(v), 3),
                         AUC_p5=round(np.percentile(v, 5), 3), AUC_p95=round(np.percentile(v, 95), 3)))
out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, "learning_curve.tsv"), sep="\t", index=False)
print(out.to_string())
