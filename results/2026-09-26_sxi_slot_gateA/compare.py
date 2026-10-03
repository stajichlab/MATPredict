"""Gate A: compare Cryptococcus calls, slot (9b8d945) vs baseline (088c055)."""
import os, sys, yaml, collections
R = os.path.dirname(os.path.abspath(__file__))
NEW, OLD = os.path.join(R, os.environ.get("NEWRUN", "9b8d945"), "runs"), os.path.join(R, "088c055", "runs")
FOCUS = {"GCA_002221985.1": "A1-35-8", "GCA_002217545.1": "cng10", "GCA_002220035.1": "MW_RSA852",
         "GCA_002222025.1": "cng4", "GCA_000836315.1": "LA55", "GCA_014964225.1": "WM1802"}

def calls(d, g):
    p = os.path.join(d, g, "detection_report.yaml")
    if not os.path.exists(p):
        return None
    r = yaml.safe_load(open(p)) or {}
    return sorted((x["idiomorph"], x["confidence"], x["locus_class"]) for x in (r.get("detected") or []))

def genotype(c):
    return "none" if not c else "+".join(sorted({x[0] for x in c}))

genomes = sorted(set(os.listdir(NEW)) & set(os.listdir(OLD)))
tally = collections.Counter(); lost = []; changed = []; gained = []; other = []
for g in genomes:
    a, b = calls(OLD, g), calls(NEW, g)
    if a is None or b is None:
        tally["missing report"] += 1; continue
    ga, gb = genotype(a), genotype(b)
    if a == b: tally["identical"] += 1
    elif ga == gb: tally["same genotype, calls differ"] += 1; other.append((g, a, b))
    elif ga == "none": tally["gained"] += 1; gained.append((g, gb, b))
    elif gb == "none": tally["lost"] += 1; lost.append((g, ga))
    else: tally["genotype changed"] += 1; changed.append((g, ga, gb))
print(f"genomes compared: {len(genomes)}  (new {len(os.listdir(NEW))}, old {len(os.listdir(OLD))})")
for k, v in tally.most_common(): print(f"  {k}: {v}")
print("\nfocus genomes (baseline -> slot):")
for g in genomes:
    acc = "_".join(g.split("_")[:2])
    if acc in FOCUS:
        print(f"  {FOCUS[acc]:10s} {g}: {calls(OLD, g)} -> {calls(NEW, g)}")
for name, rows in (("LOST", lost), ("GENOTYPE CHANGED", changed), ("GAINED", gained), ("SAME GENOTYPE, CALLS DIFFER", other)):
    print(f"\n{name} ({len(rows)}):")
    for r in rows: print("  ", *r)
