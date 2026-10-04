#!/usr/bin/env python3
"""Patterns in the rooted sexP/sexM gene trees (gene_tree_{full,hmg}.rooted.nwk).

For each tree and idiomorph:
- main clade = the largest clade with no tip of the other idiomorph (outgroup
  tips allowed); its size, UFBoot and the outgroup tips inside it;
- genus monophyly among that idiomorph's tips (genera with >= 2 tips; a genus is
  monophyletic when the MRCA of its tips holds no same-idiomorph tip of another
  genus; outgroup tips ignored);
- median patristic distance between 2,000 random tip pairs of the main clade;
- Umbelopsis: the MRCA of its tips, and the sister group of that MRCA.
Also names the source genome of the outgroup tips nested in the main clades.
Prints to stdout (saved as patterns.txt).
"""
import csv, re, collections, statistics, random
from Bio import Phylo

T = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/tree"
tips = {r["tip"]: r for r in csv.DictReader(open(f"{T}/tips.tsv"), delimiter="\t")}
samp = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}


def idio(x):
    return tips[x.name]["idiomorph"]


def genus(x):
    if idio(x) == "outgroup":
        return "OUTGROUP"
    n = tips[x.name]["name"] or ""
    return n.split()[0] if n and not re.match(r"^\d", n) else "record"


def sup(c):
    n = str(c.name or "")
    return float(n.split("/")[-1]) if "/" in n else None


def species_of(tip):
    asm = "_".join(tip.split("_")[1:3])
    for k, v in samp.items():
        if k.startswith(asm):
            return f"{v['SPECIES']} ({v['ORDER']})"
    return "?"


nested = set()
for aln in ("full", "hmg"):
    t = Phylo.read(f"{T}/gene_tree_{aln}.rooted.nwk", "newick")
    L = t.get_terminals()
    print(f"\n===== {aln}")
    for want, other in (("Plus", "Minus"), ("Minus", "Plus")):
        allw = [x for x in L if idio(x) == want]
        best = max((c for c in t.find_clades() if not any(idio(x) == other for x in c.get_terminals())),
                   key=lambda c: sum(idio(x) == want for x in c.get_terminals()))
        main = [x for x in best.get_terminals() if idio(x) == want]
        outg = [x.name for x in best.get_terminals() if idio(x) == "outgroup"]
        if len(outg) <= 10:
            nested.update(outg)
        by = collections.defaultdict(list)
        for x in allw:
            by[genus(x)].append(x)
        mono, tested, nonmono = 0, 0, []
        for g, xs in by.items():
            if len(xs) < 2 or g == "record":
                continue
            tested += 1
            gs = [genus(y) for y in t.common_ancestor(xs).get_terminals() if idio(y) == want]
            if all(gg == g for gg in gs):
                mono += 1
            else:
                nonmono.append((g, len(xs)))
        random.seed(1)
        d = [t.distance(*random.sample(main, 2)) for _ in range(2000)]
        q = statistics.quantiles(d, n=4)
        print(f"{want}: {len(allw)} tips; main clade {len(main)} (UFBoot {sup(best)}), outgroup tips in it {len(outg)}")
        print(f"   genera monophyletic {mono}/{tested}; not: {sorted(nonmono, key=lambda x: -x[1])}")
        print(f"   patristic distance in main clade: median {statistics.median(d):.2f} (IQR {q[0]:.2f}-{q[2]:.2f})")
        um = [x for x in allw if genus(x) == "Umbelopsis"]
        if um:
            m = t.common_ancestor(um)
            path = t.get_path(m)
            parent = path[-2] if len(path) > 1 else t.root
            inside = collections.Counter(genus(x) for x in m.get_terminals())
            sis = collections.Counter(genus(x) for x in parent.get_terminals() if x not in set(m.get_terminals()))
            print(f"   Umbelopsis MRCA: {dict(inside)}; sister group: {sum(sis.values())} tips {sis.most_common(5)}")
            for x in m.get_terminals():
                if genus(x) != "Umbelopsis":
                    r = tips[x.name]
                    print(f"     non-Umbelopsis tip in Umbelopsis clade: {r['label']} ({r['genome']})")

print("\n== outgroup tips nested in a main clade (<= 10 per clade)")
for n in sorted(nested):
    print(" ", n, "->", species_of(n))
