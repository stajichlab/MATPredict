"""Compare the Dothideomycetes panel before (f81dad1) and after (dca6ccc) the first class records."""
import collections
import csv
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).parent
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"
order = {r["ASMID"]: r["ORDER"] or "?" for r in csv.DictReader(open(SAMPLES))}


def load(tag):
    out = {}
    for d in sorted((HERE / tag / "runs").iterdir()):
        p = d / "detection_report.yaml"
        if not p.exists():
            out[d.name] = None
            continue
        r = yaml.safe_load(p.read_text())
        calls = [x for x in (r.get("detected") or []) if x["family"] == "Ascomycota:MAT"]
        out[d.name] = calls
    return out


def summary(tag, runs):
    called = sum(1 for c in runs.values() if c)
    conf = collections.Counter(x["confidence"] for c in runs.values() if c for x in c)
    idio = collections.Counter(x["idiomorph"] for c in runs.values() if c for x in c)
    both = sum(1 for c in runs.values() if c and {x["idiomorph"] for x in c} >= {"MAT1-1", "MAT1-2"})
    print(f"{tag}: genomes {len(runs)} no_report {sum(1 for c in runs.values() if c is None)} "
          f"called {called} loci {sum(len(c) for c in runs.values() if c)} conf {dict(conf)} "
          f"idiomorph {dict(idio)} genomes_with_both {both}")


def genotype(calls):
    if not calls:
        return "none"
    return "+".join(sorted({f"{x['idiomorph']}/{x['confidence']}" for x in calls}))


a, b = load("f81dad1"), load("dca6ccc")
summary("before f81dad1", a)
summary("after  dca6ccc", b)

by_order = collections.defaultdict(lambda: [0, 0, 0, 0, 0])  # n, called_before, called_after, high_before, high_after
changes = []
for g in sorted(set(a) | set(b)):
    o = order.get(g, "?")
    ca, cb = a.get(g), b.get(g)
    row = by_order[o]
    row[0] += 1
    row[1] += bool(ca)
    row[2] += bool(cb)
    row[3] += bool(ca) and any(x["confidence"] == "high" for x in ca)
    row[4] += bool(cb) and any(x["confidence"] == "high" for x in cb)
    if genotype(ca) != genotype(cb):
        changes.append((o, g, genotype(ca), genotype(cb)))

print("\norder                 n  called_before  called_after  high_before  high_after")
for o, (n, c1, c2, h1, h2) in sorted(by_order.items(), key=lambda kv: -kv[1][0]):
    print(f"{o:20s} {n:3d} {c1:14d} {c2:13d} {h1:12d} {h2:11d}")

print(f"\ngenomes whose idiomorph/confidence set changed: {len(changes)}")
for c in changes:
    print("  ", *c)
