"""Compare the candidate deterministic L-INS-i rebuild (run-algcand, 4f29df3)
with the shipped classifier (PR #9 R4 runs, 41bd471) on the same inputs.

Only the Mucoromycota classifier HMMs (sexM, sexP, P1) and the manifest's gate
threshold differ between the two arms. Every call, label and confidence change
is listed, with classifier margins where the report carries them.
"""
import collections
import csv
import statistics
import sys
from pathlib import Path

import yaml

R = Path("/bigdata/stajichlab/jstajich/projects/MATPredict/results")
HERE = R / "2026-09-29_aligner_default"
R4 = R / "2026-09-29_r4_paralog"
sys.path.insert(0, str(R / "2026-09-29_strain_labels_and_absidia"))
from parse_labels import label  # noqa: E402

ARMS = {
    "Mucoromycota293": (R4 / "Mucoromycota_41bd471/runs", HERE / "Mucoromycota_4f29df3/runs"),
    "LCG": (R4 / "lcg_runs", HERE / "lcg_runs"),
    "Jena": (R4 / "jena/runs_scaffolds", HERE / "jena/runs_scaffolds"),
}


def load(d: Path) -> dict:
    out = {}
    for f in d.glob("*/detection_report.yaml"):
        try:
            out[f.parent.name] = yaml.safe_load(f.read_text()) or {}
        except Exception as exc:  # noqa: BLE001
            print("unreadable", f, exc, file=sys.stderr)
    return out


def mcalls(rep: dict) -> list[dict]:
    return [c for c in (rep.get("detected") or [])
            if str(c.get("family", "")).startswith("Mucoromycota")]


def key(c):
    return (c.get("family"), c.get("contig"), c.get("start"), c.get("end"),
            c.get("idiomorph"), c.get("confidence"))


def margin(c):
    clf = c.get("idiomorph_classifier") or {}
    return clf.get("margin")


def overlap(a, b):
    return a.get("contig") == b.get("contig") and not (
        a.get("end", 0) < b.get("start", 0) or a.get("start", 0) > b.get("end", 0))


rows, summary, deltas = [], [], []
for arm, (bdir, ndir) in ARMS.items():
    base, new = load(bdir), load(ndir)
    genomes = sorted(set(base) & set(new))
    cb = sum(1 for g in genomes if mcalls(base[g]))
    cn = sum(1 for g in genomes if mcalls(new[g]))
    ng = nl = nc = 0
    for g in genomes:
        b, n = mcalls(base[g]), mcalls(new[g])
        for c in b:
            same = [x for x in n if overlap(x, c)]
            if same and margin(c) is not None and margin(same[0]) is not None:
                deltas.append(margin(same[0]) - margin(c))
            if any(key(x) == key(c) for x in n):
                continue
            if same:
                nc += 1
                x = same[0]
                rows.append(dict(arm=arm, genome=g, change="changed", contig=c.get("contig"),
                                 start=c.get("start"), end=c.get("end"),
                                 before=f"{c.get('idiomorph')}/{c.get('confidence')} m={margin(c)}",
                                 after=f"{x.get('idiomorph')}/{x.get('confidence')} m={margin(x)}"))
            else:
                nl += 1
                rows.append(dict(arm=arm, genome=g, change="lost", contig=c.get("contig"),
                                 start=c.get("start"), end=c.get("end"),
                                 before=f"{c.get('idiomorph')}/{c.get('confidence')} m={margin(c)}",
                                 after=""))
        for x in n:
            if not any(overlap(x, c) for c in b):
                ng += 1
                rows.append(dict(arm=arm, genome=g, change="gained", contig=x.get("contig"),
                                 start=x.get("start"), end=x.get("end"), before="",
                                 after=f"{x.get('idiomorph')}/{x.get('confidence')} m={margin(x)}"))
    summary.append(f"{arm}: genomes compared {len(genomes)} (shipped {len(base)}, candidate "
                   f"{len(new)}); called {cb} -> {cn}; calls lost {nl}, gained {ng}, changed {nc}")

if deltas:
    summary.append(f"classifier margin change on matched calls (n={len(deltas)}): min "
                   f"{min(deltas):.1f}, max {max(deltas):.1f}, mean |d| "
                   f"{statistics.mean(abs(d) for d in deltas):.2f} bits")

with open(HERE / "changes.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["arm", "genome", "change", "contig", "start", "end",
                                       "before", "after"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# Clean strain labels (LCG): exclude disputed, misidentified and leaked strains.
SL = R / "2026-09-29_strain_labels_and_absidia"
excl = {r["genome"] for f in ("disputed_labels.tsv", "misidentified_strains.tsv")
        for r in csv.DictReader(open(SL / f), delimiter="\t") if r.get("genome")}
leak = {r["org"]: r["leakage"] for r in csv.DictReader(
    open(R / "2026-09-28_lcg_holdout/per_genome.tsv"), delimiter="\t")}
base, new = load(ARMS["LCG"][0]), load(ARMS["LCG"][1])


def verdict(rep, lab):
    idios = {c.get("idiomorph") for c in mcalls(rep)}
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
summary.append(f"LCG clean strain labels (n={len(lab_rows)}): shipped "
               f"{dict(collections.Counter(x[2] for x in lab_rows))} -> candidate "
               f"{dict(collections.Counter(x[3] for x in lab_rows))}")
for g, lab, vb, vn in lab_rows:
    if vb != vn:
        summary.append(f"  label change {g} ({lab}): {vb} -> {vn}")

zlog = sorted(Path("/bigdata/stajichlab/jstajich/projects/MATPredict/logs").glob("zygo23_29222587.log"))
summary.append("Zygo 23: " + (" | ".join(l.strip() for l in zlog[0].read_text().splitlines()
                                           if "locus on truth" in l) if zlog else "log missing"))
(HERE / "compare_output.txt").write_text("\n".join(summary) + "\n")
print("\n".join(summary))
