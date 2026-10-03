"""Compare R4 (41bd471) runs with their baselines and score strain labels.

Arms (reports read directly):
  Mucoromycota 293: baseline results/2026-09-29_gate_threshold/Mucoromycota_3f755db
                    (PR #9 code = 9a458de, the R4 parent's code) vs Mucoromycota_41bd471
  LCG Mucoromycotina 621: baseline results/2026-09-28_lcg_holdout/runs2 (076afe4) vs lcg_runs
  Jena 64 scaffolds: baseline results/2026-09-28_mucor_jena_holdout/runs_scaffolds (076afe4)
                    vs jena/runs_scaffolds
LCG/Jena baselines predate the per-build gate threshold (99.9 vs 100 bits) and the
report-only two_idiomorphs field; every changed call is attributed by its
withheld_reason so R4 (paralog_class) effects are separated.
"""
import collections
import csv
import sys
from pathlib import Path

import yaml

R = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/results")
HERE = R / "2026-09-29_r4_paralog"
sys.path.insert(0, str(R / "2026-09-29_strain_labels_and_absidia"))
from parse_labels import label  # noqa: E402

ARMS = {
    "Mucoromycota293": (R / "2026-09-29_gate_threshold/Mucoromycota_3f755db/runs",
                        HERE / "Mucoromycota_41bd471/runs"),
    "LCG": (R / "2026-09-28_lcg_holdout/runs2", HERE / "lcg_runs"),
    "Jena": (R / "2026-09-28_mucor_jena_holdout/runs_scaffolds", HERE / "jena/runs_scaffolds"),
}


def load(d: Path) -> dict:
    out = {}
    for f in d.glob("*/detection_report.yaml"):
        try:
            out[f.parent.name] = yaml.safe_load(f.read_text()) or {}
        except Exception as exc:  # noqa: BLE001
            print("unreadable", f, exc, file=sys.stderr)
    return out


def calls(rep: dict) -> list[tuple]:
    return sorted(
        (c.get("family"), c.get("contig"), c.get("start"), c.get("end"), c.get("idiomorph"),
         c.get("confidence"))
        for c in (rep.get("detected") or []) if str(c.get("family", "")).startswith("Mucoromycota"))


def paralog_withheld(rep: dict) -> list[dict]:
    return [s for s in (rep.get("suppressed_loci") or [])
            if s.get("withheld_reason") == "paralog_class"]


def overlaps(c, s):
    return c[1] == s.get("contig") and not (c[3] < s.get("start", 0) or c[2] > s.get("end", 0))


rows, summary = [], []
for arm, (base_dir, new_dir) in ARMS.items():
    base, new = load(base_dir), load(new_dir)
    genomes = sorted(set(base) & set(new))
    n_gain = n_lost = n_changed = n_par = 0
    called_b = sum(1 for g in genomes if calls(base[g]))
    called_n = sum(1 for g in genomes if calls(new[g]))
    for g in genomes:
        b, n = calls(base[g]), calls(new[g])
        par = paralog_withheld(new[g])
        n_par += len(par)
        for s in par:
            clf = s.get("idiomorph_classifier") or {}
            rows.append(dict(arm=arm, genome=g, change="withheld_paralog_class",
                             contig=s.get("contig"), start=s.get("start"), end=s.get("end"),
                             idiomorph=s.get("idiomorph"),
                             mat_scores=clf.get("scores"), paralog_scores=clf.get("paralog_scores"),
                             genes=",".join(s.get("genes_found") or []),
                             only_call_lost=bool(b) and not n))
        if b == n:
            continue
        for c in b:
            if c in n:
                continue
            same_place = [x for x in n if x[1] == c[1] and not (x[3] < c[2] or x[2] > c[3])]
            if same_place:
                n_changed += 1
                rows.append(dict(arm=arm, genome=g, change="changed", contig=c[1], start=c[2],
                                 end=c[3], idiomorph=f"{c[4]}/{c[5]} -> "
                                 f"{same_place[0][4]}/{same_place[0][5]}", mat_scores="",
                                 paralog_scores="", genes="", only_call_lost=False))
            else:
                n_lost += 1
                why = "paralog_class" if any(overlaps(c, s) for s in par) else "other"
                rows.append(dict(arm=arm, genome=g, change=f"lost ({why})", contig=c[1],
                                 start=c[2], end=c[3], idiomorph=f"{c[4]}/{c[5]}",
                                 mat_scores="", paralog_scores="", genes="",
                                 only_call_lost=not n))
        for c in n:
            if c not in b and not any(x[1] == c[1] and not (x[3] < c[2] or x[2] > c[3]) for x in b):
                n_gain += 1
                rows.append(dict(arm=arm, genome=g, change="gained", contig=c[1], start=c[2],
                                 end=c[3], idiomorph=f"{c[4]}/{c[5]}", mat_scores="",
                                 paralog_scores="", genes="", only_call_lost=False))
    summary.append(f"{arm}: genomes compared {len(genomes)} (baseline {len(base)}, R4 {len(new)}); "
                   f"called {called_b} -> {called_n}; paralog_class withheld {n_par}; "
                   f"calls lost {n_lost}, gained {n_gain}, changed {n_changed}")

with open(HERE / "changes.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["arm"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# Strain labels (LCG), excluding disputed, misidentified and leaked strains.
SL = R / "2026-09-29_strain_labels_and_absidia"
excl = {r["genome"] for f in ("disputed_labels.tsv", "misidentified_strains.tsv")
        for r in csv.DictReader(open(SL / f), delimiter="\t")}
leak = {r["org"]: r["leakage"] for r in csv.DictReader(
    open(R / "2026-09-28_lcg_holdout/per_genome.tsv"), delimiter="\t")}
base, new = load(ARMS["LCG"][0]), load(ARMS["LCG"][1])


def verdict(rep, lab):
    idios = {c[4] for c in calls(rep)}
    if not idios:
        return "uncalled"
    if idios == {lab}:
        return "agree"
    if {"Plus", "Minus"} <= idios:
        return "two_idiomorph"
    if idios == {"undetermined"}:
        return "undetermined"
    return "agree_plus_other" if lab in idios else "disagree"


lab_rows = []
for g in sorted(set(base) & set(new)):
    lab, _ = label(g)
    if not lab or g in excl or leak.get(g) != "clean":
        continue
    lab_rows.append((g, lab, verdict(base[g], lab), verdict(new[g], lab)))
cb = collections.Counter(x[2] for x in lab_rows)
cn = collections.Counter(x[3] for x in lab_rows)
summary.append(f"LCG clean strain labels (n={len(lab_rows)}; disputed/misidentified/leaked "
               f"excluded): baseline {dict(cb)} -> R4 {dict(cn)}")
for g, lab, vb, vn in lab_rows:
    if vb != vn:
        summary.append(f"  label change {g} ({lab}): {vb} -> {vn}")

zy = (HERE / "zygo23_41bd471" / "score.txt")
summary.append("Zygo 23: " + (zy.read_text().strip().replace("\n", " | ") if zy.exists()
                              else "score file missing (see log)"))
(HERE / "compare_output.txt").write_text("\n".join(summary) + "\n")
print("\n".join(summary))
