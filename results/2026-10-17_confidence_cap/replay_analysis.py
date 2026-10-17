#!/usr/bin/env python
"""Replay candidate confidence-cap rules on loci_features.tsv (v0.6.0 reports). Plain python.
    replay_analysis.py loci_features.tsv ../2026-10-03_ascomycota_v060/genomes.tsv
Prints markdown tables. A rule is a predicate on a `high` locus; it would drop that call to `medium`."""
import csv, sys
from collections import Counter, defaultdict

loc = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
order = {r["genome"]: (r["class_"], r["order"]) for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
f = lambda x: float(x) if x not in ("", None) else None
for r in loc:
    for k in ("best_cov", "best_ident", "max_cov", "max_ident", "edge_dist", "best_gene_edge"):
        r[k] = f(r[k])
    r["n_flank"] = int(r["n_flank"]); r["n_modelled_core"] = int(r["n_modelled_core"])
    r["order"] = order.get(r["genome"], ("?", "?"))[1] or "(none)"
    r["class_"] = order.get(r["genome"], ("?", "?"))[0]
high = [r for r in loc if r["confidence"] == "high"]
print(f"{len(loc)} loci, {len(high)} high\n")

def lt(v, t): return v is not None and v < t
RULES = {}
for t in (50, 60, 70, 80, 90):
    RULES[f"P{t}: best-core coverage <{t}"] = lambda r, t=t: lt(r["best_cov"], t)
for t in (60, 80):
    RULES[f"P{t}+locus-edge: best-core coverage <{t} and LOCUS segment within 500 bp of a contig end"] = lambda r, t=t: lt(r["best_cov"], t) and lt(r["edge_dist"], 500)
for t in (40, 45, 50):
    RULES[f"W{t}: best-core identity <{t}, no flank gene"] = lambda r, t=t: lt(r["best_ident"], t) and r["n_flank"] == 0
    RULES[f"W{t}max: highest core identity <{t}, no flank gene"] = lambda r, t=t: lt(r["max_ident"], t) and r["n_flank"] == 0
RULES["W45+edge: best-core identity <45, no flank, within 500 bp of contig end"] = lambda r: lt(r["best_ident"], 45) and r["n_flank"] == 0 and lt(r["edge_dist"], 500)
RULES["W45+gene_only: best-core identity <45, no flank, locus_class idiomorph_gene_only"] = lambda r: lt(r["best_ident"], 45) and r["n_flank"] == 0 and r["locus_class"] == "idiomorph_gene_only"
#RULES["PW(80,45): P80 OR W45"] = lambda r: lt(r["best_cov"], 80) and lt(r["best_ident"], 45) or (lt(r["best_ident"], 45) and r["n_flank"] == 0)
#RULES["PW2: (cov<80 AND ident<45) AND no flank"] = lambda r: lt(r["best_cov"], 80) and lt(r["best_ident"], 45) and r["n_flank"] == 0
#RULES["P80&edge OR W45"] = lambda r: (lt(r["best_cov"], 80) and lt(r["edge_dist"], 500)) or (lt(r["best_ident"], 45) and r["n_flank"] == 0)


# ---- the implemented criteria (confidence_cap.py), replayed ----
RULES["PE: partial_edge = best-core cov <80, no core gene >=90, best gene within 500 bp of a contig end"] = lambda r: lt(r["best_cov"], 80) and lt(r["max_cov"], 90) and lt(r["best_gene_edge"], 500)
RULES["WK: weak_identity = best-core identity <45, no flank gene"] = lambda r: lt(r["best_ident"], 45) and r["n_flank"] == 0
IMPL = lambda r: RULES["PE: partial_edge = best-core cov <80, no core gene >=90, best gene within 500 bp of a contig end"](r) or RULES["WK: weak_identity = best-core identity <45, no flank gene"](r)

def by(rows, key, top=12):
    c = Counter(key(r) for r in rows)
    return ", ".join(f"{k} {n}" for k, n in c.most_common(top))

def main():
    truth = {r["genome"]: r for r in csv.DictReader(open(sys.argv[3]), delimiter="\t")} if len(sys.argv) > 3 else {}
    print("## Candidate rules on the v0.6.0 `high` calls\n")
    print("| rule | drops | % of high | mat_locus | gene_only | homothallic |\n|---|---|---|---|---|---|")
    for n, p in RULES.items():
        d = [r for r in high if p(r)]
        c = Counter(r["locus_class"] for r in d)
        print(f"| {n} | {len(d)} | {100*len(d)/len(high):.1f} | {c['mat_locus']} | {c['idiomorph_gene_only']} | {c['homothallic_candidate']} |")
    print("\n## Truth set (BioSample MAT1-1/MAT1-2 stated, Ascomycota:MAT/MATtub high loci)\n")
    T = []
    for r in high:
        t = truth.get(r["genome"])
        if t and t["stated"] in ("MAT1-1", "MAT1-2") and r["family"] in ("Ascomycota:MAT", "Ascomycota:MATtub"):
            r["ok"] = r["idiomorph"] == t["stated"]; r["bot"] = "Botrytis" in t["organism"]; T.append(r)
    print(f"{len(T)} high loci in truth genomes: right {sum(r['ok'] for r in T)} (Botrytis {sum(r['ok'] and r['bot'] for r in T)}), wrong {sum(not r['ok'] for r in T)} (Botrytis {sum((not r['ok']) and r['bot'] for r in T)})\n")
    print("| rule | demoted right | demoted wrong (non-Botrytis) | demoted wrong (Botrytis, fixed by PR #95) |\n|---|---|---|---|")
    for n, p in RULES.items():
        d = [r for r in T if p(r)]
        print(f"| {n} | {sum(r['ok'] for r in d)} | {sum((not r['ok']) and not r['bot'] for r in d)} | {sum((not r['ok']) and r['bot'] for r in d)} |")
    print("\n## Sensitivity of the weak-identity threshold (high, no flank gene)\n")
    print("| identity below | loci |\n|---|---|")
    for t in (35, 40, 45, 50, 55, 60, 65, 70, 80):
        print(f"| {t} | {sum(1 for r in high if r['n_flank']==0 and lt(r['best_ident'], t))} |")
    print("\nHistogram of best-core identity, high gene-only loci (5-point bins):\n")
    h = Counter(int(r["best_ident"] // 5) * 5 for r in high if r["n_flank"] == 0 and r["best_ident"] is not None)
    print(", ".join(f"{k}:{h[k]}" for k in sorted(h)))
    print("\nSensitivity of the edge distance (partial_edge: best cov <80, no full core gene, gene-level distance):\n")
    print("| edge below (bp) | loci |\n|---|---|")
    for e in (100, 500, 1000, 2500, 5000, 10000):
        print(f"| {e} | {sum(1 for r in high if lt(r['best_cov'],80) and lt(r['max_cov'],90) and lt(r['best_gene_edge'], e))} |")
    d = [r for r in high if IMPL(r)]
    pe = [r for r in high if RULES["PE: partial_edge = best-core cov <80, no core gene >=90, best gene within 500 bp of a contig end"](r)]
    wk = [r for r in high if RULES["WK: weak_identity = best-core identity <45, no flank gene"](r)]
    print(f"\n## Implemented rule (partial_weak = PE or WK): {len(d)} of {len(high)} high ({100*len(d)/len(high):.2f}%); PE {len(pe)}, WK {len(wk)}, both {len(set(map(id,pe))&set(map(id,wk)))}\n")
    for lab, rows in (("all", d), ("partial_edge", pe), ("weak_identity", wk)):
        print(f"**{lab}** ({len(rows)})\n")
        print("- family: " + by(rows, lambda r: r["family"].replace("Ascomycota:", "")))
        print("- locus_class: " + by(rows, lambda r: r["locus_class"]))
        print("- idiomorph: " + by(rows, lambda r: r["idiomorph"]))
        print("- class: " + by(rows, lambda r: r["class_"]))
        print("- order: " + by(rows, lambda r: r["order"], 15))
        print("- genome with another call of the same family at `high`: " + str(sum(1 for r in rows if any(o is not r and o["genome"] == r["genome"] and o["family"] == r["family"] and o["confidence"] == "high" for o in loc))))
        print()
    print("\n## High calls per family (denominator) and share dropped\n")
    fh = Counter(r["family"] for r in high); fd = Counter(r["family"] for r in d)
    print("| family | high | dropped | % |\n|---|---|---|---|")
    for k, n in fh.most_common():
        print(f"| {k} | {n} | {fd[k]} | {100*fd[k]/n:.1f} |")
    print("\n## High calls per order (top, by drops)\n")
    oh = Counter(r["order"] for r in high); od = Counter(r["order"] for r in d)
    print("| order | high | dropped | % |\n|---|---|---|---|")
    for k, n in od.most_common(15):
        print(f"| {k} | {oh[k]} | {n} | {100*n/oh[k]:.1f} |")
    print("\n## locus_class of dropped\n")
    lh = Counter(r["locus_class"] for r in high); ld = Counter(r["locus_class"] for r in d)
    print("| locus_class | high | dropped |\n|---|---|---|")
    for k, n in lh.most_common():
        print(f"| {k} | {n} | {ld[k]} |")
    w = [r for r in high if r["locus_class"]=="idiomorph_gene_only"]
    print(f"\nGene-only `high` loci: {len(w)}; dropped {sum(1 for r in w if IMPL(r))} ({100*sum(1 for r in w if IMPL(r))/len(w):.1f}%); kept {sum(1 for r in w if not IMPL(r))}.")
    open("replay_dropped_loci.tsv", "w").write("genome\tfamily\torder\tidiomorph\tlocus_class\tcriteria\tbest_gene\tbest_ident\tbest_cov\tbest_aa\tbest_gene_edge\ttruth\n" + "".join(
        f"{r['genome']}\t{r['family']}\t{r['order']}\t{r['idiomorph']}\t{r['locus_class']}\t{'+'.join(x for x,pp in (('partial_edge', r in pe), ('weak_identity', r in wk)) if pp)}\t{r['best_gene']}\t{r['best_ident']}\t{r['best_cov']}\t{r['best_aa']}\t{r['best_gene_edge']}\t{truth.get(r['genome'],{}).get('verdict','')}\n" for r in d))

if __name__ == "__main__":
    main()
