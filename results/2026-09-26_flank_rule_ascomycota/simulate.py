"""Simulate alternative flank-carried rules on the 108 changed calls.
Per call, per core gene: detection's coordinates (changed.tsv core_detail) and the best
in-locus tblastn e-value under detect's settings (classified3.tsv).
Rule A = current (ALL core hits inside flank span +-3 kb -> low, else withheld).
Rule B(W, E) = keep at low if the STRONGEST core hit (lowest e-value) has e <= E and lies
within W bp of the flank span on the same contig; else withhold."""
import csv, collections
rows = list(csv.DictReader(open("changed.tsv"), delimiter="\t"))
ev = collections.defaultdict(dict)
for r in csv.DictReader(open("classified3.tsv"), delimiter="\t"):
    if r["here_evalue"]:
        ev[(r["panel"], r["genome"], r["mat_family"], r["call"])][r["gene"]] = float(r["here_evalue"])
def best_core(r):
    key = (r["panel"], r["genome"], r["mat_family"], f'{r["contig"]}:{r["start"]}-{r["end"]}')
    lo, hi = int(r["flank_lo"]), int(r["flank_hi"])
    out = []
    for d in r["core_detail"].split(";"):
        g, pos = d.split(":")[0], d.split(":")[1]
        s, e = map(int, pos.split("-"))
        dist = max(lo - e, s - hi, 0)
        out.append((ev[key].get(g, 99.0), dist, g))
    return min(out)
res = collections.Counter()
out = csv.writer(open("simulated.tsv", "w"), delimiter="\t")
out.writerow(["panel", "genome", "order", "mat_family", "call", "ruleA", "best_core_gene", "best_core_evalue", "best_core_dist_bp",
              "B_3kb_noE", "B_20kb_noE", "B_20kb_E1e-5", "B_3kb_E1e-5"])
for r in rows:
    e, dist, g = best_core(r)
    b = lambda W, E: "low" if dist <= W and e <= E else "withheld"
    rec = [r["outcome"], b(3000, 99), b(20000, 99), b(20000, 1e-5), b(3000, 1e-5)]
    out.writerow([r["panel"], r["genome"], r["order"], r["mat_family"], f'{r["contig"]}:{r["start"]}-{r["end"]}', rec[0], g, e, dist] + rec[1:])
    grp = "Asco" if r["panel"] == "cap6" else "MucoroG"
    for name, v in zip(["A_current", "B_3kb_noE", "B_20kb_noE", "B_20kb_E1e-5", "B_3kb_E1e-5"], rec):
        res[(grp, name, v)] += 1
for grp in ["Asco", "MucoroG"]:
    print(grp, {n: (res[(grp, n, "low")], res[(grp, n, "withheld")]) for n in ["A_current", "B_3kb_noE", "B_20kb_noE", "B_20kb_E1e-5", "B_3kb_E1e-5"]})
