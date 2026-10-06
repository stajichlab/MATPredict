"""Replay stricter MAT-gene-gate variants on existing detect reports.

Read-only. Reproduces the gate in src/MATPredict/detect/mat_gene_gate.py
(remote polish-scope-cuts 076afe4) from report fields, checks that every
reported classified call passes it (reproduction), then applies variants:

  V1: flank route needs >=2 distinct flanks modelled at >=40% AND at least one
      that is not rnhA.
  V2(X): a partial_locus call whose best model score is below the absolute
      threshold is withheld unless margin >= X.
  V2b(X): ANY partial_locus call is withheld unless margin >= X (added: the
      target calls pass the absolute route, so V2 cannot reach them).
  V3: V1 with the absolute route unchanged -- identical to V1 by definition.

Usage: python3 replay.py > replay_output.txt
"""
import collections
import csv
import glob
import os
import sys

import yaml

R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results"
SETS = {
    "LCG": f"{R}/2026-09-28_lcg_holdout/runs2/*/detection_report.yaml",
    "Mucoro293": f"{R}/2026-09-28_next_fixes/Mucoromycota_a9fb69c/runs/*/detection_report.yaml",
    "Jena": f"{R}/2026-09-28_mucor_jena_holdout/runs_scaffolds/*/detection_report.yaml",
    "Zygo23": f"{R}/2026-09-28_next_fixes/zygo23_a9fb69c/scaffold/runs/*/detection_report.yaml",
}
TARGETS = {"Circinella_angarensis_RSA_198_Plus", "Circinella_umbellata_RSA_505_Plus",
           "Thamnostylum_repens_RSA_459_Plus", "Backusella_lamprospora_NRRL_6044_Plus",
           "Pilaira_anomala_RSA_1997_Plus"}
DISPUTED = {"Ellisomyces_anomalus_RSA_581-", "Gilbertella_persicaria_CBS_442.64-",
            "Pirella_circinans_RSA_622-"}
MIN_SCORE, MIN_ID, MIN_GENES = 100.0, 40.0, 2
MODELLED = {"polished_agree", "polished_disagree", "polished_single"}
LICHT = ("Circinella", "Lichtheimia", "Rhizomucor", "Thermomucor", "Dichotomocladium",
         "Absidia_corymbifera", "Absidia_ramosa", "Zychaea")


def flanks(call):
    return sorted({e["gene"] for e in call.get("gene_evidence", [])
                   if e["role"].startswith("flanking")
                   and (e.get("status") in MODELLED or e.get("method") == "diamond_proteome")
                   and e.get("identity") is not None and e["identity"] >= MIN_ID})


def gate_route(call):
    """'score' | 'flank' | 'fail' under the current gate; None if not gated."""
    clf = call.get("idiomorph_classifier") or {}
    if not clf or call.get("split_locus"):
        return None
    best = max((clf.get("scores") or {}).values(), default=0.0)
    if clf.get("classifier_input") == "model" and best >= MIN_SCORE:
        return "score"
    return "flank" if len(flanks(call)) >= MIN_GENES else "fail"


def variant_withholds(call, v):
    route = gate_route(call)
    if route is None:
        return False
    clf = call["idiomorph_classifier"]
    best = max((clf.get("scores") or {}).values(), default=0.0)
    margin = clf.get("margin") or 0.0
    partial = call.get("locus_class") == "partial_locus"
    if v in ("V1", "V3"):
        if route == "flank":
            fl = flanks(call)
            return not (len(fl) >= 2 and any(g != "rnhA" for g in fl))
        return False
    if v.startswith("V2b"):
        return partial and margin < float(v.split("_")[1])
    if v.startswith("V2"):
        return partial and route != "score" and margin < float(v.split("_")[1])
    raise ValueError(v)


def load_labels():
    lab = {}
    p = f"{R}/2026-09-29_strain_labels_and_absidia"
    for f in glob.glob(f"{p}/*.tsv"):
        with open(f) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            cols = rd.fieldnames or []
            gcol = next((c for c in cols if c in ("genome", "organism", "org", "name")), None)
            lcol = next((c for c in cols if "label" in c.lower()), None)
            tcol = next((c for c in cols if c in ("tier", "leak_tier", "leakage")), None)
            if not gcol or not lcol:
                continue
            for r in rd:
                if r[lcol] in ("Plus", "Minus"):
                    lab[r[gcol]] = (r[lcol], (r.get(tcol) or "") if tcol else "")
    return lab


def main():
    labels = load_labels()
    variants = ["V1", "V2_75", "V2_100", "V2b_75", "V2b_100"]
    calls = []
    for s, pat in SETS.items():
        for f in sorted(glob.glob(pat)):
            g = f.split("/")[-2]
            r = yaml.safe_load(open(f)) or {}
            for c in r.get("detected") or []:
                calls.append((s, g, c))
    repro = collections.Counter()
    for s, g, c in calls:
        repro[(s, gate_route(c))] += 1
    print("reproduction (current gate route of reported calls):")
    for k, v in sorted(repro.items(), key=str):
        print("  ", k, v)
    print("  a 'fail' route here would mean the replay does not reproduce the gate\n")
    for v in variants:
        print(f"=== {v}")
        wh = [(s, g, c) for s, g, c in calls if variant_withholds(c, v)]
        per = collections.Counter(s for s, g, c in wh)
        print("  withheld per set:", dict(per))
        print("  targets withheld:", sorted(g for s, g, c in wh if g in TARGETS))
        for s, g, c in wh:
            clf = c["idiomorph_classifier"]
            best = max((clf.get("scores") or {}).values(), default=0.0)
            tag = []
            if g in TARGETS:
                tag.append("TARGET")
            if s == "Zygo23":
                tag.append("ZYGO")
            if g.startswith(LICHT):
                tag.append("licht")
            if g in labels:
                lab = labels[g][0]
                tag.append(f"label={lab}{'(agree)' if lab == c['idiomorph'] else '(disagree)'}")
            if g in DISPUTED:
                tag.append("disputed")
            print(f"   {s:9} {g[:45]:45} {c['idiomorph']:12} {c['confidence']:6} {c['locus_class']:13} "
                  f"in={clf.get('classifier_input')} best={best:.1f} m={clf.get('margin')} "
                  f"fl={','.join(flanks(c))} {' '.join(tag)}")
        # label agreement on the clean labelled set (LCG only), excluding disputed labels
        whset = {(s2, g2, id(c2)) for s2, g2, c2 in wh}
        agree = dis = unc = 0
        byg = collections.defaultdict(list)
        for s2, g2, c2 in calls:
            if s2 == "LCG" and (s2, g2, id(c2)) not in whset:
                byg[g2].append(c2["idiomorph"])
        for g2, (lab, tier) in labels.items():
            if g2 in DISPUTED or tier != "clean":
                continue
            got = [x for x in byg.get(g2, []) if x in ("Plus", "Minus")]
            if not got:
                unc += 1
            elif lab in got and len(set(got)) == 1:
                agree += 1
            else:
                dis += 1
        print(f"  clean labelled (disputed excluded): agree {agree} disagree {dis} uncalled {unc}")
        print()


if __name__ == "__main__":
    main()
