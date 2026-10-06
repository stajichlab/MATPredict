"""Stress test: low caps (3, 2) plain vs V1/V2/V3, against a cap-OFF reference.

A reference call is a (family, contig, idiomorph) in the cap-off run.
- lost by plain: reference call absent from the plain low-cap run.
- recovered by Vx: a lost call that Vx reports with the same family+contig+label.
- false promotion by Vx: a Vx call absent from the plain run that is NOT a
  cap-off call (different contig, different label, or not in cap-off at all),
  or carries a paralog class.
Also runs the early-diverging (Mortierellomycota/Kickxellomycota) comparison.
"""
import collections
import glob
import os

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
S = os.path.join(HERE, "stress")
SUBS = {"mucoro": "mucoro/runs", "lcg": "lcg", "jena": "jena/runs_scaffolds"}


def load(d):
    out = {}
    for f in glob.glob(os.path.join(d, "*", "detection_report.yaml")):
        g = f.split("/")[-2]
        try:
            r = yaml.safe_load(open(f)) or {}
        except Exception:
            continue
        cs = set()
        meta = {}
        for x in r.get("detected") or []:
            clf = x.get("idiomorph_classifier") or {}
            k = (x.get("family"), x.get("contig"), x.get("idiomorph"))
            cs.add(k)
            meta[k] = (x.get("confidence"), clf.get("margin"), clf.get("paralog_class"),
                       x.get("start"), x.get("end"))
        wall = None
        wf = os.path.join(os.path.dirname(f), "wall_seconds")
        if os.path.exists(wf):
            try:
                wall = float(open(wf).read())
            except ValueError:
                pass
        out[g] = (cs, meta, wall)
    return out


def wall_total(runs, genomes):
    return sum(runs[g][2] or 0 for g in genomes)


def stress():
    lines = []
    for s, sub in SUBS.items():
        ref = load(os.path.join(S, "off_plain", sub))
        for cap in ("c3", "c2"):
            plain = load(os.path.join(S, f"{cap}_plain", sub))
            common = set(ref) & set(plain)
            lost = [(g, k) for g in common for k in ref[g][0] - plain[g][0]]
            extra_plain = [(g, k) for g in common for k in plain[g][0] - ref[g][0]]
            lines.append(f"\n== {s} {cap}: genomes {len(common)}; cap-off calls "
                         f"{sum(len(ref[g][0]) for g in common)}; lost by plain {len(lost)}; "
                         f"plain calls not in cap-off {len(extra_plain)}; "
                         f"wall off {wall_total(ref, common):.0f}s plain {wall_total(plain, common):.0f}s")
            for v in ("V1", "V2", "V3"):
                var = load(os.path.join(S, f"{cap}_{v}", sub))
                c2 = common & set(var)
                rec = [(g, k) for g, k in lost if g in c2 and k in var[g][0]]
                new = [(g, k) for g in c2 for k in var[g][0] - plain[g][0]]
                false = []
                for g, k in new:
                    para = var[g][1][k][2]
                    if k not in ref[g][0] or para:
                        false.append((g, k, var[g][1][k]))
                still_lost = [(g, k) for g, k in lost if g in c2 and k not in var[g][0]]
                dropped = [(g, k) for g in c2 for k in plain[g][0] - var[g][0]]
                lines.append(f"   {v}: recovered {len(rec)}/{len(lost)}; false promotions {len(false)}; "
                             f"plain calls dropped {len(dropped)}; wall {wall_total(var, c2):.0f}s")
                for g, k, m in false[:12]:
                    lines.append(f"      FALSE {g} {k} {m}")
                for g, k in dropped[:6]:
                    lines.append(f"      DROPPED {g} {k}")
                if v == "V3":
                    for g, k in still_lost[:8]:
                        lines.append(f"      STILL_LOST {g} {k}")
    for cfg in os.listdir(S):
        p = os.path.join(S, cfg, "zygo23", "score.txt")
        if os.path.exists(p):
            lines.append(f"\n== zygo {cfg}: " + open(p).read().strip().replace("\n", " | "))
    return lines


def early_diverging():
    lines = []
    E = os.path.join(HERE, "ed")
    for cap in ("c6", "c3"):
        for grp in ("mort_00", "kickx_00", "kickx_01"):
            p = load(os.path.join(E, f"{cap}_plain", grp, "runs"))
            v = load(os.path.join(E, f"{cap}_V3", grp, "runs"))
            c = set(p) & set(v)
            gained = [(g, k, v[g][1][k]) for g in c for k in v[g][0] - p[g][0]]
            lostc = [(g, k) for g in c for k in p[g][0] - v[g][0]]
            lines.append(f"\n== ED {cap} {grp}: genomes {len(c)}; called plain "
                         f"{sum(1 for g in c if p[g][0])} V3 {sum(1 for g in c if v[g][0])}; "
                         f"gained {len(gained)} lost {len(lostc)}; wall {wall_total(p, c):.0f}s -> "
                         f"{wall_total(v, c):.0f}s")
            for g, k, m in gained[:10]:
                lines.append(f"      GAINED {g} {k} {m}")
            for g, k in lostc[:6]:
                lines.append(f"      LOST {g} {k}")
    return lines


if __name__ == "__main__":
    out = stress() + early_diverging()
    print("\n".join(out))
