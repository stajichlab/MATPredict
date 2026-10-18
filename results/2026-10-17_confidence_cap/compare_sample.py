#!/usr/bin/env python
"""Compare --confidence-cap-policy off vs partial_weak runs on the stratified sample.

    compare_sample.py RUNS_DIR sample.tsv features.tsv truth.tsv out_prefix

Per genome: (1) the two reports must be identical apart from `confidence` / `confidence_cap` of the detected loci (anything else is a
FAILURE of the "cap changes only confidence" property); (2) loci whose confidence differs are listed with the note; (3) each dropped
locus is checked against the v0.6.0 replay prediction (features.tsv, rule PE/WK) for the same genome/family/contig/start.
Extrapolation: weight = N_h / n_h per stratum (sample.tsv).
"""
import csv, copy, sys
from collections import Counter, defaultdict
from pathlib import Path
import yaml

L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
runs, sample, feats, truthp, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
S = {r["genome"]: r for r in csv.DictReader(open(sample), delimiter="\t")}
truth = {r["genome"]: r for r in csv.DictReader(open(truthp), delimiter="\t")}
rep = {(r["genome"], r["family"], r["contig"], r["start"]): r for r in csv.DictReader(open(feats), delimiter="\t")}


def load(v, g):
    p = runs / v / g / "detection_report.yaml"
    return yaml.load(p.open(), Loader=L) if p.exists() else None


def strip(doc):
    d = copy.deepcopy(doc)
    d.pop("run", None)  # parameters (the flag itself), start time and wall time differ by construction
    for c in d.get("detected") or []:
        c.pop("confidence", None); c.pop("confidence_cap", None)
    return d


def f(x):
    return float(x) if x not in ("", None) else None


done, paired, missing, changed_other, drops, hi_total = 0, 0, [], [], [], 0
per_genome = []
for g, s in S.items():
    b = load("partial_weak", g)
    if b is None:
        missing.append(g); continue
    done += 1
    lb = b.get("detected") or []
    noted = [c for c in lb if c.get("confidence_cap")]
    a = load("off", g)
    if a is not None:
        paired += 1
        if strip(a) != strip(b):
            changed_other.append(g)
        la = a.get("detected") or []
        for ca, cb in zip(la, lb):
            if (ca["confidence"] != cb["confidence"]) != bool(cb.get("confidence_cap")):
                changed_other.append(g + ":cap-note-mismatch")
    for cb in noted:
        ca = copy.deepcopy(cb); ca["confidence"] = "high"   # the `off` confidence the note records (from: high)
        drops.append((g, ca, cb))
    n_high_off = sum(1 for c in lb if c["confidence"] == "high") + len(noted)
    hi_total += n_high_off
    per_genome.append((g, s, n_high_off, len(noted)))

print(f"{done} genomes with a partial_weak run, {paired} of them also with an off run; {len(missing)} missing; "
      f"paired genomes that differ in something other than confidence/confidence_cap: {changed_other[:8]}")
bad = [d for d in drops if d[2]["confidence"] != "medium"]
print(f"{len(drops)} loci carry a confidence_cap note; {len(bad)} of them are not at medium")
W = lambda g: float(S[g]["weight"])
hi_genomes = [x for x in per_genome if x[2] > 0]
est = lambda sel: sum(W(x[0]) for x in sel)
print(f"high loci in the sample (off): {hi_total}; weighted high loci: {sum(W(x[0]) * x[2] for x in per_genome):.0f}")
print(f"dropped loci: {len(drops)} raw, weighted {sum(W(d[0]) for d in drops):.0f}; genomes affected: {len({d[0] for d in drops})}, weighted {est([x for x in per_genome if x[3] > 0]):.0f}")
crit = Counter(); fam = Counter(); fam_w = Counter(); order = Counter(); order_w = Counter(); cls = Counter(); cls_w = Counter()
for g, ca, cb in drops:
    k = "+".join(cb["confidence_cap"]["criteria"])
    crit[k] += 1
    fam[ca["family"]] += 1; fam_w[ca["family"]] += W(g)
    o = S[g]["order"] or "(none)"; order[o] += 1; order_w[o] += W(g)
    cls[ca["locus_class"]] += 1; cls_w[ca["locus_class"]] += W(g)
print("criteria:", dict(crit)); print("family:", dict(fam)); print("locus_class:", dict(cls))
print("order (raw/weighted):", [(k, v, round(order_w[k])) for k, v in order.most_common(12)])

# agreement with the replay: match the v0.6.0 locus by genome, family, contig and start within 10 kb; did the replay rule predict it?
repg = defaultdict(list)
for r in rep.values():
    repg[(r["genome"], r["family"], r["contig"])].append(r)
def predicted(r):
    cov, mx, ident, edge, nfl = f(r["best_cov"]), f(r["max_cov"]), f(r["best_ident"]), f(r["best_gene_edge"]), int(r["n_flank"])
    pe = cov is not None and cov < 80 and mx is not None and mx < 90 and edge is not None and edge < 500
    wk = ident is not None and ident < 45 and nfl == 0
    return pe, wk
agree = Counter()
for g, ca, cb in drops:
    near = [r for r in repg[(g, ca["family"], ca["contig"])] if abs(int(r["start"]) - int(ca["start"])) < 10000]
    if not near:
        agree["no matching locus in v0.6.0 report"] += 1; continue
    r = near[0]
    pe, wk = predicted(r)
    agree[f"v0.6.0 {r['confidence']}, replay PE={pe} WK={wk}"] += 1
print("drops vs v0.6.0 replay:", dict(agree))
# truth
tr = [(g, ca, cb) for g, ca, cb in drops if g in truth]
print(f"dropped loci in truth genomes: {len(tr)}", [(truth[g]['organism'], truth[g]['stated'], ca['idiomorph']) for g, ca, cb in tr])
# stratified estimate of the number of dropped loci / high loci over all called genomes, with the usual stratified SE
def strat(valfn):
    by = defaultdict(list)
    for g, s_, nh, nd in per_genome:
        by[(s_["class_"], s_["order"])].append(valfn(nh, nd))
    tot = var = 0.0
    for k, v in by.items():
        n = len(v); N = float(S[next(g for g, s_, _, _ in per_genome if (s_["class_"], s_["order"]) == k)]["N_h"])
        m = sum(v) / n; tot += N * m
        if n > 1:
            s2 = sum((x - m) ** 2 for x in v) / (n - 1); var += N * N * (1 - n / N) * s2 / n
    return tot, var ** 0.5
for lab, fn in (("dropped loci", lambda nh, nd: nd), ("high loci (off)", lambda nh, nd: nh + nd * 0), ("genomes with a drop", lambda nh, nd: 1 if nd else 0)):
    t, se = strat(fn)
    print(f"estimated {lab} over all called genomes: {t:.0f} +- {se:.0f} (SE, stratified)")
with open(out + "_dropped.tsv", "w") as fh:
    fh.write("genome\torder\tfamily\tidiomorph\tlocus_class\tcriteria\tbest_gene\tbest_cov\tbest_ident\tedge\tweight\ttruth_stated\n")
    for g, ca, cb in drops:
        n = cb["confidence_cap"]
        fh.write(f"{g}\t{S[g]['order']}\t{ca['family']}\t{ca['idiomorph']}\t{ca['locus_class']}\t{'+'.join(n['criteria'])}\t{n['best_core_gene']}\t{n['best_core_coverage']}\t{n['best_core_identity']}\t{n['best_core_edge_distance']}\t{S[g]['weight']}\t{truth.get(g, {}).get('stated', '')}\n")
with open(out + "_per_genome.tsv", "w") as fh:
    fh.write("genome\tclass_\torder\tweight\tn_high_off\tn_dropped\n")
    for g, s, nh, nd in per_genome:
        fh.write(f"{g}\t{s['class_']}\t{s['order']}\t{s['weight']}\t{nh}\t{nd}\n")
