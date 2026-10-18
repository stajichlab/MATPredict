#!/usr/bin/env python
"""Coverage x identity grid on the truth-set `high` Ascomycota:MAT/MATtub loci (non-Botrytis). python truth_grid.py loci_features.tsv genomes.tsv agreement.tsv"""
import csv, statistics, sys
sys.argv = [sys.argv[0]] + sys.argv[1:4]
exec(open("replay_analysis.py").read().split("if __name__")[0].replace("print(f\"{len(loc)}", "pass #"))
truth = {r["genome"]: r for r in csv.DictReader(open(sys.argv[3]), delimiter="\t")}
T = []
for r in high:
    t = truth.get(r["genome"])
    if t and t["stated"] in ("MAT1-1", "MAT1-2") and r["family"] in ("Ascomycota:MAT", "Ascomycota:MATtub") and "Botrytis" not in t["organism"]:
        r["ok"] = r["idiomorph"] == t["stated"]; T.append(r)
print(f"non-Botrytis truth high loci: {len(T)}, right {sum(r['ok'] for r in T)}, wrong {sum(not r['ok'] for r in T)}")
print("median best-core coverage, all high:", statistics.median(r["best_cov"] for r in high), "; high mat_locus:", statistics.median(r["best_cov"] for r in high if r["locus_class"] == "mat_locus"))
print("\n| coverage below | identity below | demoted right | demoted wrong | all high demoted |\n|---|---|---|---|---|")
for c in (30, 50, 80):
    for i in (45, 50, 55, 60, 70):
        d = [r for r in T if r["best_cov"] < c and r["best_ident"] < i]
        print(f"| {c} | {i} | {sum(r['ok'] for r in d)} | {sum(not r['ok'] for r in d)} | {sum(1 for r in high if r['best_cov'] < c and r['best_ident'] < i)} |")
