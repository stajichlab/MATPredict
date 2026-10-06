"""Compare cap-protection variants (V1/V2/V3) against the a2fe1b4 baseline.

Baseline runs (same code/db as a2fe1b4): results/2026-09-29_umbelopsis_merge/
{Mucoromycota_f59353c/runs, lcg_runs, jena/runs_scaffolds, zygo23_f59353c}.
Controls: ctrl_base vs ctrl_V1 (Ascomycota/Basidiomycota samples; no classifier).
"""
import collections
import csv
import glob
import os
import re
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.dirname(HERE)
BASE = os.path.join(R, "2026-09-29_umbelopsis_merge")
SETS = {
    "mucoro": (os.path.join(BASE, "Mucoromycota_f59353c", "runs"), "mucoro/runs"),
    "lcg": (os.path.join(BASE, "lcg_runs"), "lcg"),
    "jena": (os.path.join(BASE, "jena", "runs_scaffolds"), "jena/runs_scaffolds"),
}


def calls(d):
    out = {}
    for f in glob.glob(os.path.join(d, "*", "detection_report.yaml")):
        g = f.split("/")[-2]
        try:
            r = yaml.safe_load(open(f)) or {}
        except Exception:
            continue
        cs = []
        for x in r.get("detected") or []:
            clf = x.get("idiomorph_classifier") or {}
            cs.append((x.get("family"), x.get("contig"), x.get("start"), x.get("end"),
                       x.get("idiomorph"), x.get("confidence"),
                       clf.get("margin"), clf.get("classifier_input"), clf.get("paralog_class")))
        wall = None
        wf = os.path.join(os.path.dirname(f), "wall_seconds")
        if os.path.exists(wf):
            try:
                wall = float(open(wf).read().strip())
            except ValueError:
                pass
        out[g] = (cs, wall)
    return out


def key(c):
    return (c[0], c[1])


def diff(base, var):
    lost, gained, changed = [], [], []
    for g in sorted(set(base) & set(var)):
        b = {key(c): c for c in base[g][0]}
        v = {key(c): c for c in var[g][0]}
        for k in b:
            if k not in v:
                # same family+contig moved? count as lost only if no overlapping same-family call
                lost.append((g, b[k]))
            elif (b[k][4], b[k][5]) != (v[k][4], v[k][5]):
                changed.append((g, b[k], v[k]))
        for k in v:
            if k not in b:
                gained.append((g, v[k]))
    return lost, gained, changed


def labels():
    lab = {}
    dis = set()
    for fn in ("disputed_labels.tsv", "misidentified_strains.tsv"):
        p = os.path.join(R, "2026-09-29_strain_labels_and_absidia", fn)
        for row in open(p):
            if row.startswith("genome") or not row.strip() or row.startswith("Also"):
                continue
            dis.add(row.split("\t")[0])
    for g in glob.glob(os.path.join(R, "2026-09-28_lcg_holdout", "curator_table.tsv")):
        pass
    return dis


def main():
    variants = [v for v in ("V1", "V2", "V3") if os.path.isdir(os.path.join(HERE, v))]
    dis = labels()
    rows = []
    for s, (bdir, vsub) in SETS.items():
        base = calls(bdir)
        for v in variants:
            var = calls(os.path.join(HERE, v, vsub))
            lost, gained, changed = diff(base, var)
            nb = sum(1 for g in set(base) & set(var) if base[g][0])
            nv = sum(1 for g in set(base) & set(var) if var[g][0])
            wb = sorted(base[g][1] for g in set(base) & set(var) if base[g][1] is not None)
            wv = sorted(var[g][1] for g in set(base) & set(var) if var[g][1] is not None)
            med = lambda a: a[len(a) // 2] if a else None
            print(f"\n== {s} {v}: genomes compared {len(set(base) & set(var))}; called {nb} -> {nv}; "
                  f"lost {len(lost)} gained {len(gained)} changed {len(changed)}; "
                  f"median wall {med(wb)} -> {med(wv)} s; total {sum(wb):.0f} -> {sum(wv):.0f} s")
            for g, c in lost:
                print("   LOST   ", g, c)
            for g, c in gained:
                print("   GAINED ", g, c)
            for g, b, c in changed:
                print("   CHANGED", g, b[4:6], "->", c[4:6], "margin", b[6], "->", c[6])
            rows.append((s, v, len(set(base) & set(var)), nb, nv, len(lost), len(gained), len(changed),
                         med(wb), med(wv), sum(wb), sum(wv)))
    # controls
    for L in ("asco50", "basidio50"):
        b = calls(os.path.join(HERE, "ctrl_base", L, "runs"))
        v = calls(os.path.join(HERE, "ctrl_V1", L, "runs"))
        lost, gained, changed = diff(b, v)
        print(f"\n== control {L}: genomes {len(set(b) & set(v))}; lost {len(lost)} gained {len(gained)} "
              f"changed {len(changed)}")
    # zygo
    for v in variants:
        p = os.path.join(HERE, v, "zygo23", "score.txt")
        if os.path.exists(p):
            print(f"\n== zygo {v}:", open(p).read().strip().replace("\n", " | "))
    # exemptions
    for v in variants:
        p = os.path.join(HERE, f"exempt_{v}.log")
        if os.path.exists(p):
            n = collections.Counter()
            for l in open(p):
                f = l.rstrip("\n").split("\t")
                if len(f) < 6:
                    continue
                n["strong"] += 1
                n["strong_in_top6"] += int(f[4])
                n["strong_allowed"] += int(f[5])
                n["newly_allowed"] += int(f[5]) and not int(f[4])
            print(f"\n== exemptions {v}:", dict(n))
    with open(os.path.join(HERE, "summary.tsv"), "w") as fo:
        fo.write("set\tvariant\tgenomes\tcalled_base\tcalled_var\tlost\tgained\tchanged\t"
                 "median_wall_base\tmedian_wall_var\ttotal_wall_base\ttotal_wall_var\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")


if __name__ == "__main__":
    main()
