"""Tree-independent check: idiomorph label vs the core genes the call itself found.

For every called locus in the early-diverging reports: if genes_found holds sexM but
not sexP, the expected label is Minus; sexP but not sexM, Plus. Counts agreement.
"""
import collections, glob, os, yaml
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_early_diverging"
c = collections.Counter(); ex = collections.defaultdict(collections.Counter)
for group in ("Mucoromycota", "Mortierellomycota", "Kickxellomycota"):
    for rep in glob.glob(f"{R}/{group}/runs/*/detection_report.yaml"):
        r = yaml.safe_load(open(rep)) or {}
        for x in r.get("detected") or []:
            g = set(x.get("genes_found") or [])
            exp = "Minus" if "sexM" in g and "sexP" not in g else "Plus" if "sexP" in g and "sexM" not in g else "both/none"
            lab = x.get("idiomorph")
            k = (group, exp, lab); c[k] += 1
print("group\tcore genes imply\tlabel\tcalls")
for k, v in sorted(c.items()): print("\t".join(map(str, k)) + f"\t{v}")
