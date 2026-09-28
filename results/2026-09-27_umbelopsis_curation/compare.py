"""Compare Mucoromycota scans before (4457c3a) and after (3adfd90) the Umbelopsis records.

Usage: compare.py BEFORE_RUNS AFTER_RUNS  -> prints per-Umbelopsidaceae table and whole-scan changes.
"""
import csv
import glob
import os
import sys

import yaml

samp = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}


def load(d):
    out = {}
    for f in glob.glob(f"{d}/*/detection_report.yaml"):
        a = os.path.basename(os.path.dirname(f))
        r = yaml.safe_load(open(f))
        calls = []
        for x in r.get("detected") or []:
            ev = {e["gene"]: e for e in x.get("gene_evidence") or []}
            core = [g for g in ("sexM", "sexP") if g in ev]
            clf = x.get("idiomorph_classifier") or {}
            calls.append(dict(
                contig=x["contig"], start=x["start"], idiom=x["idiomorph"], conf=x["confidence"],
                cls=x["locus_class"],
                core=";".join(f"{g}:{ev[g]['status']}:{ev[g]['identity']}" for g in core),
                clf=clf.get("method", "-"), clf_input=clf.get("input", clf.get("classifier_input", "model" if clf else "-")),
                margin=clf.get("margin"), routing=r.get("routing_mode")))
        out[a] = (r.get("routing_mode"), calls)
    return out


before, after = load(sys.argv[1]), load(sys.argv[2])
umb = sorted(a for a in set(before) | set(after) if samp.get(a, {}).get("FAMILY") == "Umbelopsidaceae")
print("## Umbelopsidaceae (before 4457c3a -> after 3adfd90)")
for a in umb:
    s = samp[a]
    b = before.get(a, (None, []))
    f = after.get(a, (None, []))
    fmt = lambda c: "; ".join(f"{x['idiom']}/{x['conf']}/{x['cls']} [{x['core']}] clf={x['clf_input']} m={x['margin']}" for x in c) or "uncalled"
    print(f"{s.get('SPECIES')} {s.get('STRAIN')} {a}\n  before ({b[0]}): {fmt(b[1])}\n  after  ({f[0]}): {fmt(f[1])}")

print("\n## Whole scan")
common = set(before) & set(after)
chg = {"unchanged": 0, "gained": 0, "lost": 0, "label": 0, "conf": 0}
for a in sorted(common):
    b = [(x["contig"], x["idiom"], x["conf"]) for x in before[a][1]]
    f = [(x["contig"], x["idiom"], x["conf"]) for x in after[a][1]]
    if b == f:
        chg["unchanged"] += 1
        continue
    bl, fl = {x[1] for x in b}, {x[1] for x in f}
    kind = "gained" if not b and f else "lost" if b and not f else "label" if bl != fl else "conf"
    chg[kind] += 1
    fam = samp.get(a, {}).get("FAMILY")
    print(f"  {kind:7s} {fam} {samp.get(a, {}).get('SPECIES')} {a}: {b} -> {f}")
print("genomes compared", len(common), chg)
print("called before", sum(1 for a in common if before[a][1]), "after", sum(1 for a in common if after[a][1]))
