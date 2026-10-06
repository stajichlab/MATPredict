#!/usr/bin/env python3
"""Univariate separation and held-out-clade / held-out-genome classifier comparison.

Inputs (same directory): labelled_loci.tsv, protein_features.tsv, labelled.faa,
emb/*.npy + emb.ids (ESM-2 mean-pooled), Outputs: univariate_<set>.tsv,
cv_summary_<set>_<cv>.tsv, cv_per_clade_<set>_<cv>.tsv, oof_<set>_<cv>.tsv, caax_distance_profile.tsv.

Sets: H (headline) = genomes whose mating label does not come from CAAX evidence and is not
putative (grades A and B); ALL = every labelled genome.
CV: LOCO = leave one order out; LOSO = leave one genome out (leaks close homologs; secondary).
Uncertainty: cluster bootstrap over genomes (2000), permutation of labels within genome (2000).
"""
import itertools
import os
import sys
import warnings

import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from scipy.stats import mannwhitneyu
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(sys.argv[0]))
RNG = np.random.default_rng(20261005)
NBOOT = int(os.environ.get("NBOOT", 2000))

# ---------------------------------------------------------------- data
L = pd.read_csv(f"{HERE}/labelled_loci.tsv", sep="\t")
P = pd.read_csv(f"{HERE}/protein_features.tsv", sep="\t")
L = L.merge(P.drop(columns=["prot_len"]), on="locus_id", how="inner")
L = L[L.label.isin(["mating", "other"]) & (L.has_prot == 1)].reset_index(drop=True)
L["y"] = (L.label == "mating").astype(int)
L["genome"] = L.asm.str[:15]
seq = {}
k = None
for line in open(f"{HERE}/labelled.faa"):
    if line.startswith(">"):
        k = line[1:].strip(); seq[k] = ""
    else:
        seq[k] += line.strip()
L["seq"] = L.locus_id.map(seq)
ids = [x.strip() for x in open(f"{HERE}/emb/emb.ids")]
EMB = {}
for f in sorted(os.listdir(f"{HERE}/emb")):
    if f.endswith(".npy"):
        M = np.load(f"{HERE}/emb/{f}")
        d = dict(zip(ids, M))
        EMB[f.split(".")[1]] = d
BIG = 1e7
L["caax10"] = (L.nT_10kb > 0).astype(int)
L["log_dT"] = np.log10(1 + L.d_T.clip(upper=BIG))
L["log_dHx"] = np.log10(1 + L.d_Hx.clip(upper=BIG))
L["log_dHx_xo"] = np.log10(1 + L.d_Hx_xo.clip(upper=BIG))
L["log_dHD_xo"] = np.log10(1 + L.d_HD_xo.clip(upper=BIG))
L["log_dFLANK_xo"] = np.log10(1 + L.d_FLANK_xo.clip(upper=BIG))
L["log_dHD"] = np.log10(1 + L.d_HD.clip(upper=BIG))
L["log_dFLANK"] = np.log10(1 + L.d_FLANK.clip(upper=BIG))
L["introns"] = (L.n_cds - 1).clip(lower=0)
L["neg_dT"] = -L.log_dT
FEAT_PREC = ["caax10", "log_dT", "nT_20kb", "nT_50kb", "nHx_10kb_xo", "log_dHx_xo"]  # nT2 omitted: motif widened using labelled positives
FEAT_MAT = ["log_dHD_xo", "n_HD_20kb_xo", "log_dFLANK_xo", "n_FLANK_50kb_xo"]
FEAT_PROT = ["prot_len", "introns", "tm", "hydrophobic_frac", "pf02076_score", "pf02076_cov"]
UNI = FEAT_PREC + FEAT_MAT + FEAT_PROT + ["nL_10kb", "nT2_10kb", "nHx_10kb", "log_dHx", "log_dHD", "n_HD_20kb", "log_dFLANK", "n_FLANK_50kb"]


# ---------------------------------------------------------------- metrics
def auc(y, s):
    return roc_auc_score(y, s) if 0 < y.sum() < len(y) else np.nan


def ap(y, s):
    return average_precision_score(y, s) if 0 < y.sum() < len(y) else np.nan


def sens95(y, s):
    if not 0 < y.sum() < len(y):
        return np.nan
    fpr, tpr, _ = roc_curve(y, s)
    return tpr[fpr <= 0.05].max()


def within_auc(y, s, g):
    out = []
    for gg in np.unique(g):
        m = g == gg
        if 0 < y[m].sum() < m.sum():
            out.append(roc_auc_score(y[m], s[m]))
    return np.mean(out) if out else np.nan, len(out)


def boot_idx(g, B):
    ug = np.unique(g)
    members = {u: np.where(g == u)[0] for u in ug}
    for _ in range(B):
        pick = RNG.choice(ug, len(ug), replace=True)
        yield np.concatenate([members[u] for u in pick])


def boot_ci(y, s, g, fn, B=NBOOT):
    vals = []
    for ix in boot_idx(g, B):
        v = fn(y[ix], s[ix])
        if not np.isnan(v):
            vals.append(v)
    return (np.percentile(vals, 2.5), np.percentile(vals, 97.5)) if len(vals) > 50 else (np.nan, np.nan)


def perm_p(y, s, g, B=NBOOT):
    obs = auc(y, s)
    ge = 0
    members = [np.where(g == u)[0] for u in np.unique(g)]
    for _ in range(B):
        yp = y.copy()
        for m in members:
            yp[m] = RNG.permutation(y[m])
        if auc(yp, s) >= obs:
            ge += 1
    return (ge + 1) / (B + 1)


# ---------------------------------------------------------------- univariate
def univariate(D, tag):
    y, g = D.y.values, D.genome.values
    rows = []
    for f in UNI:
        s = D[f].values.astype(float)
        a = auc(y, s)
        lo, hi = boot_ci(y, s, g, auc, 1000)
        wa, wn = within_auc(y, s, g)
        try:
            p = mannwhitneyu(s[y == 1], s[y == 0]).pvalue
        except ValueError:
            p = np.nan
        rows.append(dict(feature=f, n_pos=int(y.sum()), n_neg=int((1 - y).sum()), n_genomes=len(np.unique(g)),
                         mean_mating=s[y == 1].mean(), mean_other=s[y == 0].mean(), AUC=a, AUC_lo=lo, AUC_hi=hi,
                         cliffs_delta=2 * a - 1, within_genome_AUC=wa, n_genomes_both=wn, MW_p=p))
    pd.DataFrame(rows).round(3).to_csv(f"{HERE}/univariate_{tag}.tsv", sep="\t", index=False)


# ---------------------------------------------------------------- models
def lr(C=0.1):
    return make_pipeline(StandardScaler(), LogisticRegression(C=C, class_weight="balanced", max_iter=5000))


def kmer_matrix(seqs, k):
    aa = "ACDEFGHIKLMNPQRSTVWY"
    idx = {"".join(p): i for i, p in enumerate(itertools.product(aa, repeat=k))}
    X = np.zeros((len(seqs), len(idx)))
    for i, s in enumerate(seqs):
        for j in range(len(s) - k + 1):
            t = idx.get(s[j:j + k])
            if t is not None:
                X[i, t] += 1
        X[i] /= max(1, X[i].sum())
    return X


aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score, aligner.extend_gap_score, aligner.mode = -11, -1, "local"


def align_matrix(seqs):
    n = len(seqs)
    self_ = np.array([aligner.score(s, s) for s in seqs])
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            v = aligner.score(seqs[i], seqs[j]) / min(self_[i], self_[j])
            M[i, j] = M[j, i] = v
    return M


def fit_predict(name, tr, te, D, ctx):
    y = D.y.values
    if name == "rule_caax10":
        return D.caax10.values[te].astype(float)
    if name == "rule_dist":
        return D.neg_dT.values[te]
    if name == "pf02076_only":
        return D.pf02076_score.values[te]
    if name.startswith("LR_") or name.startswith("GBM_") or name.startswith("RF_"):
        cols = {"prec": FEAT_PREC, "all": FEAT_PREC + FEAT_MAT + FEAT_PROT, "noprec": FEAT_MAT + FEAT_PROT,
                "prot": FEAT_PROT, "mat": FEAT_MAT}[name.split("_", 1)[1]]
        X = D[cols].values.astype(float)
        if name.startswith("LR_"):
            m = lr(0.3)
        elif name.startswith("GBM_"):
            m = HistGradientBoostingClassifier(max_depth=3, max_iter=100, learning_rate=0.1, min_samples_leaf=3, class_weight="balanced")
        else:
            m = RandomForestClassifier(300, max_depth=4, min_samples_leaf=2, class_weight="balanced", random_state=1)
        m.fit(X[tr], y[tr])
        return m.predict_proba(X[te])[:, 1]
    if name.startswith("ESM_LR_") or name.startswith("ESM_kNN_") or name.startswith("ESMprec_LR_"):
        mod = name.split("_")[-1]
        X = np.vstack([EMB[mod][i] for i in D.locus_id])
        if name.startswith("ESMprec"):
            X = np.hstack([X, D[FEAT_PREC].values.astype(float) * 3.0])
        if name.startswith("ESM_kNN"):
            Z = X / np.linalg.norm(X, axis=1, keepdims=True)
            S = Z[te] @ Z[tr].T
            pos = S[:, y[tr] == 1]
            neg = S[:, y[tr] == 0]
            kp, kn = min(3, pos.shape[1]), min(3, neg.shape[1])
            return np.sort(pos, 1)[:, -kp:].mean(1) - np.sort(neg, 1)[:, -kn:].mean(1)
        m = lr(0.01)
        m.fit(X[tr], y[tr])
        return m.predict_proba(X[te])[:, 1]
    if name.startswith("kmer"):
        k_ = int(name[4])
        X = ctx["kmer"][k_]
        m = lr(0.01)
        m.fit(X[tr], y[tr])
        return m.predict_proba(X[te])[:, 1]
    if name == "alignNN":
        S = ctx["align"][np.ix_(te, tr)]
        pos, neg = S[:, y[tr] == 1], S[:, y[tr] == 0]
        kp, kn = min(3, pos.shape[1]), min(3, neg.shape[1])
        return np.sort(pos, 1)[:, -kp:].mean(1) - np.sort(neg, 1)[:, -kn:].mean(1)
    raise ValueError(name)


def run_cv(D, tag, cvname, models):
    D = D.reset_index(drop=True)
    y = D.y.values
    gcol = "order" if cvname == "LOCO" else "genome"
    groups = D[gcol].values
    ix = [GIDX[i] for i in D.locus_id]
    ctx = {"kmer": {2: GK[2][ix], 3: GK[3][ix]}, "align": GA[np.ix_(ix, ix)]}
    oof = {m: np.full(len(D), np.nan) for m in models}
    for g in np.unique(groups):
        te = np.where(groups == g)[0]
        tr = np.where(groups != g)[0]
        if len(np.unique(y[tr])) < 2:
            continue
        for m in models:
            oof[m][te] = fit_predict(m, tr, te, D, ctx)
    gen = D.genome.values
    keep = ~np.isnan(oof[models[0]])
    out, per = [], []
    rule = oof["rule_caax10"]
    for m in models:
        s = oof[m]
        a, ap_, s95 = auc(y[keep], s[keep]), ap(y[keep], s[keep]), sens95(y[keep], s[keep])
        alo, ahi = boot_ci(y[keep], s[keep], gen[keep], auc)
        plo, phi = boot_ci(y[keep], s[keep], gen[keep], ap)
        slo, shi = boot_ci(y[keep], s[keep], gen[keep], sens95)
        wa, wn = within_auc(y[keep], s[keep], gen[keep])
        pp = perm_p(y[keep], s[keep], gen[keep], 1000)
        # paired bootstrap of AUC and AP difference to the rule
        d_auc, d_ap = [], []
        for ix in boot_idx(gen[keep], NBOOT):
            yy = y[keep][ix]
            if 0 < yy.sum() < len(yy):
                d_auc.append(auc(yy, s[keep][ix]) - auc(yy, rule[keep][ix]))
                d_ap.append(ap(yy, s[keep][ix]) - ap(yy, rule[keep][ix]))
        out.append(dict(model=m, cv=cvname, set=tag, n=int(keep.sum()), n_pos=int(y[keep].sum()), n_neg=int((1 - y[keep]).sum()),
                        n_genomes=len(np.unique(gen[keep])), AUC=a, AUC_lo=alo, AUC_hi=ahi, AP=ap_, AP_lo=plo, AP_hi=phi,
                        sens_at_95spec=s95, sens_lo=slo, sens_hi=shi, within_genome_AUC=wa, n_genomes_both=wn, perm_p=pp,
                        dAUC_vs_rule=np.mean(d_auc), dAUC_lo=np.percentile(d_auc, 2.5), dAUC_hi=np.percentile(d_auc, 97.5),
                        dAP_vs_rule=np.mean(d_ap), dAP_lo=np.percentile(d_ap, 2.5), dAP_hi=np.percentile(d_ap, 97.5)))
        for g in np.unique(groups):
            mm = (groups == g) & keep
            per.append(dict(model=m, cv=cvname, set=tag, heldout=g, n_pos=int(y[mm].sum()), n_neg=int((1 - y[mm]).sum()),
                            AUC=auc(y[mm], s[mm]), AP=ap(y[mm], s[mm])))
    pd.DataFrame(out).round(3).to_csv(f"{HERE}/cv_summary_{tag}_{cvname}.tsv", sep="\t", index=False)
    pd.DataFrame(per).round(3).to_csv(f"{HERE}/cv_per_clade_{tag}_{cvname}.tsv", sep="\t", index=False)
    o = D[["locus_id", "order", "genome", "label"]].copy()
    for m in models:
        o[m] = oof[m]
    o.round(4).to_csv(f"{HERE}/oof_{tag}_{cvname}.tsv", sep="\t", index=False)


GIDX = {i: n for n, i in enumerate(L.locus_id)}
GA = align_matrix(list(L.seq))
GK = {2: kmer_matrix(L.seq, 2), 3: kmer_matrix(L.seq, 3)}


def add_conservation():
    A = GA
    ms, ng = [], []
    for i in range(len(L)):
        o = np.where((L.order.values == L.order.values[i]) & (L.genome.values != L.genome.values[i]))[0]
        if len(o) == 0:
            ms.append(np.nan); ng.append(np.nan); continue
        ms.append(A[i, o].max())
        ng.append(len({L.genome.values[j] for j in o if A[i, j] >= 0.5}))
    L["cons_maxsim_order"], L["cons_ngenomes_order"] = ms, ng


MODELS = ["rule_caax10", "rule_dist", "pf02076_only", "LR_prec", "LR_all", "LR_noprec", "LR_prot", "GBM_all", "RF_all",
          "kmer2", "kmer3", "alignNN"] + [f"ESM_LR_{m}" for m in EMB] + [f"ESM_kNN_{m}" for m in EMB] + [f"ESMprec_LR_{m}" for m in EMB]

if __name__ == "__main__":
    add_conservation()
    UNI += ["cons_maxsim_order", "cons_ngenomes_order"]
    SETS = {"H": L[L.headline].reset_index(drop=True), "ALL": L.reset_index(drop=True)}
    # sensitivity set: headline without Coprinopsis "other" copies (many receptors in that genome)
    SETS["Hnocc"] = L[L.headline & ~((L.species == "Coprinopsis cinerea") & (L.y == 0))].reset_index(drop=True)
    # complete-model subset: protein >= 200 aa and PF02076 covering >= 50 percent of the HMM
    SETS["Hq"] = L[L.headline & (L.prot_len >= 200) & (L.pf02076_cov >= 0.5)].reset_index(drop=True)
    only = os.environ.get("ONLY")
    for tag, D in SETS.items():
        if only and tag != only:
            continue
        print(tag, len(D), D.y.sum(), D.genome.nunique(), flush=True)
        univariate(D, tag)
        for cv in ("LOCO", "LOSO"):
            run_cv(D, tag, cv, MODELS)
            print("done", tag, cv, flush=True)
    if only and only != "PROFILE":
        sys.exit(0)
    # CAAX distance profile
    rows = []
    for tag, D in SETS.items():
        for lab in (1, 0):
            d = D[D.y == lab]
            for w in (1000, 2000, 5000, 10000, 20000, 50000, 100000):
                rows.append(dict(set=tag, label="mating" if lab else "other", window_bp=w, n=len(d), n_within=int((d.d_T <= w).sum()),
                                 frac=round((d.d_T <= w).mean(), 3)))
    pd.DataFrame(rows).to_csv(f"{HERE}/caax_distance_profile.tsv", sep="\t", index=False)
