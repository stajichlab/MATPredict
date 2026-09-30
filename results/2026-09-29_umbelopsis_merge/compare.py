"""Compare curation-umbelopsis rebuilt deterministically (f59353c) with the
shipped PR #9 classifier runs (aligner_default candidate = 8eff88e HMMs).
"""
import collections
import csv
import sys
from pathlib import Path

import yaml

R = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/results")
HERE = R / "2026-09-29_umbelopsis_merge"
R4 = R / "2026-09-29_aligner_default"
sys.path.insert(0, str(R / "2026-09-29_strain_labels_and_absidia"))
from parse_labels import label  # noqa: E402

ARMS = {
    "Mucoromycota293": (R4 / "Mucoromycota_4f29df3/runs", HERE / "Mucoromycota_f59353c/runs"),
    "LCG": (R4 / "lcg_runs", HERE / "lcg_runs"),
    "Jena": (R4 / "jena/runs_scaffolds", HERE / "jena/runs_scaffolds"),
}
SAMPLES = {r["ASMID"]: r for r in csv.DictReader(
    open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}


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


def detail(rep: dict) -> list[str]:
    out = []
    for c in rep.get("detected") or []:
        if not str(c.get("family", "")).startswith("Mucoromycota"):
            continue
        clf = c.get("idiomorph_classifier") or {}
        out.append(f"{c.get('idiomorph')}/{c.get('confidence')} {c.get('locus_class')} "
                   f"input={clf.get('input') or c.get('classifier_input')} "
                   f"margin={clf.get('margin')} {c.get('contig')}:{c.get('start')}-{c.get('end')}")
    return out


def paralog_withheld(rep: dict) -> list[dict]:
    return [s for s in (rep.get("suppressed_loci") or [])
            if s.get("withheld_reason") == "paralog_class"]


def same_place(a, b):
    return a[1] == b[1] and not (a[3] < b[2] or a[2] > b[3])


rows, summary, loaded = [], [], {}
for arm, (base_dir, new_dir) in ARMS.items():
    base, new = load(base_dir), load(new_dir)
    loaded[arm] = (base, new)
    genomes = sorted(set(base) & set(new))
    n_gain = n_lost = n_changed = n_span = 0
    par_b = sum(len(paralog_withheld(base[g])) for g in genomes)
    par_n = sum(len(paralog_withheld(new[g])) for g in genomes)
    called_b = sum(1 for g in genomes if calls(base[g]))
    called_n = sum(1 for g in genomes if calls(new[g]))
    for g in genomes:
        b, n = calls(base[g]), calls(new[g])
        if b == n:
            continue
        for c in b:
            if c in n:
                continue
            sp = [x for x in n if same_place(c, x)]
            if sp and (sp[0][4], sp[0][5]) == (c[4], c[5]):
                n_span += 1
                continue
            if sp:
                n_changed += 1
                ch = f"{c[4]}/{c[5]} -> {sp[0][4]}/{sp[0][5]}"
                rows.append(dict(arm=arm, genome=g, change="changed", contig=c[1], start=c[2],
                                 end=c[3], detail=ch))
            else:
                n_lost += 1
                why = [s.get("withheld_reason") for s in (new[g].get("suppressed_loci") or [])
                       if s.get("contig") == c[1] and not (c[3] < s.get("start", 0)
                                                           or c[2] > s.get("end", 0))]
                rows.append(dict(arm=arm, genome=g, change="lost", contig=c[1], start=c[2],
                                 end=c[3], detail=f"{c[4]}/{c[5]}; now {why}; "
                                 f"only_call_lost={not n}"))
        for c in n:
            if not any(same_place(c, x) for x in b):
                n_gain += 1
                rows.append(dict(arm=arm, genome=g, change="gained", contig=c[1], start=c[2],
                                 end=c[3], detail=f"{c[4]}/{c[5]}"))
    summary.append(f"{arm}: compared {len(genomes)} (base {len(base)}, new {len(new)}); "
                   f"called {called_b} -> {called_n}; paralog_class withheld {par_b} -> {par_n}; "
                   f"lost {n_lost}, gained {n_gain}, changed label/confidence {n_changed}, "
                   f"span-only {n_span}")

with open(HERE / "changes.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["arm", "genome", "change", "contig", "start", "end",
                                       "detail"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# Umbelopsidaceae per genome (Mucoromycota 293 arm).
base, new = loaded["Mucoromycota293"]
urows = []
for g in sorted(new):
    if SAMPLES.get(g, {}).get("FAMILY") != "Umbelopsidaceae":
        continue
    urows.append(dict(genome=g, species=SAMPLES[g].get("SPECIES"),
                      pr9_r4=" | ".join(detail(base.get(g, {}))) or "uncalled",
                      branch=" | ".join(detail(new[g])) or "uncalled"))
with open(HERE / "umbelopsidaceae.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(urows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(urows)
summary.append(f"Umbelopsidaceae called: PR9-shipped {sum(1 for r in urows if r['pr9_r4'] != 'uncalled')}"
               f"/{len(urows)} -> branch {sum(1 for r in urows if r['branch'] != 'uncalled')}"
               f"/{len(urows)}")

# P1 withholdings in the new runs.
for arm, (_, new) in loaded.items():
    for g, rep in sorted(new.items()):
        for s in paralog_withheld(rep):
            clf = s.get("idiomorph_classifier") or {}
            summary.append(f"  P1 withheld [{arm}] {g} {s.get('contig')}:{s.get('start')}-"
                           f"{s.get('end')} mat={clf.get('scores')} paralog={clf.get('paralog_scores')}"
                           f" still_called={bool(calls(rep))}")

# Strain labels (LCG), excluding disputed, misidentified and leaked strains.
SL = R / "2026-09-29_strain_labels_and_absidia"
excl = {r["genome"] for f in ("disputed_labels.tsv", "misidentified_strains.tsv")
        for r in csv.DictReader(open(SL / f), delimiter="\t") if r.get("genome")}
leak = {r["org"]: r["leakage"] for r in csv.DictReader(
    open(R / "2026-09-28_lcg_holdout/per_genome.tsv"), delimiter="\t")}
base, new = loaded["LCG"]


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
summary.append(f"LCG clean strain labels (n={len(lab_rows)}): PR9-shipped {dict(cb)} -> branch {dict(cn)}")
for g, lab, vb, vn in lab_rows:
    if vb != vn:
        summary.append(f"  label change {g} ({lab}): {vb} -> {vn}")

zy = HERE / "zygo23_f59353c" / "score.txt"
summary.append("Zygo 23: " + (zy.read_text().strip().replace("\n", " | ") if zy.exists()
                              else "score file missing (see log)"))
(HERE / "compare_output.txt").write_text("\n".join(summary) + "\n")
print("\n".join(summary))
