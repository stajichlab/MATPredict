import csv, sys, collections
D, REF, V = sys.argv[1:4]
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
work = {r[0]: r[2] for r in csv.reader(open(f"{D}/work.tsv"), delimiter="\t")}
gen = {r["genome"]: r for r in rd(f"{V}/genomes.tsv")}
def load(pol):
    loci = {}; genomes = set()
    for r in rd(f"{D}/loci_{pol}.tsv"):
        if r["family"] == "__genome__": genomes.add(r["genome"]); continue
        if r["family"] != "Ascomycota:MAT": continue
        loci[(r["genome"], r["contig"], r["start"])] = r
    return loci, genomes
L, gl = load("legacy"); C, gc = load("core_first")
both = gl & gc
print("genomes with both reports:", len(both), "fix:", sum(work[g]=="fix" for g in both), "ctl:", sum(work[g]=="ctl" for g in both))
def lab(r): return r["idiomorph"] if r else "(absent)"
def summarize(kind):
    gs = {g for g in both if work[g]==kind}
    keys = [k for k in set(L)|set(C) if k[0] in gs]
    trans = collections.Counter(); confch = collections.Counter(); chg_g = set(); lab_g = set(); byorder = collections.Counter(); bygenus = collections.Counter(); nloc = collections.Counter()
    gbyorder = collections.defaultdict(set); gbygenus = collections.defaultdict(set)
    for k in keys:
        a, b = L.get(k), C.get(k)
        nloc["legacy"] += a is not None; nloc["core_first"] += b is not None
        g = gen.get(k[0], {}); order = g.get("order") or g.get("class_") or "NA"; genus = (g.get("species") or "?").split()[0]
        changed = False
        if a and b:
            if a["idiomorph"] != b["idiomorph"]:
                trans[(a["idiomorph"], b["idiomorph"])] += 1; changed = True
            if a["confidence"] != b["confidence"]:
                confch[(a["confidence"], b["confidence"])] += 1; chg_g.add(k[0])
            elif a["locus_class"] != b["locus_class"]:
                confch[("class:"+a["locus_class"], "class:"+b["locus_class"])] += 1; chg_g.add(k[0])
        else:
            trans[(lab(a), lab(b))] += 1; changed = True
        if changed:
            chg_g.add(k[0]); lab_g.add(k[0]); byorder[order] += 1; bygenus[genus] += 1; gbyorder[order].add(k[0]); gbygenus[genus].add(k[0])
    print(f"\n== {kind}: {len(gs)} genomes; loci legacy {nloc['legacy']} core_first {nloc['core_first']}; genomes with a label/loci change {len(lab_g)}; with any label/conf/class change {len(chg_g)}")
    for t, n in trans.most_common(): print("  label", t, n)
    for t, n in confch.most_common(): print("  conf", t, n)
    print("  genomes by order:", sorted(((o, len(s)) for o, s in gbyorder.items()), key=lambda x: -x[1])[:15])
    print("  genomes by genus:", sorted(((o, len(s)) for o, s in gbygenus.items()), key=lambda x: -x[1])[:20])
    return lab_g
summarize("fix"); summarize("ctl")
ctlg = {g for g in both if work[g]=="ctl"}
nd = lambda r: {x: y for x, y in (r or {}).items() if x != "policy"}
print("control loci differing in any column except policy:", sum(1 for k in set(L)|set(C) if k[0] in ctlg and nd(L.get(k)) != nd(C.get(k))))
def gcall(loci, g):
    labs = {r["idiomorph"] for k, r in loci.items() if k[0]==g}
    if not labs: return "none"
    return next(iter(labs)) if len(labs)==1 else "mixed"
v6 = collections.defaultdict(set)
for r in rd(f"{V}/loci.tsv"):
    if r["family_called"]=="Ascomycota:MAT": v6[r["genome"]].add(r["idiomorph"])
def v6call(g):
    labs = v6.get(g, set())
    return "none" if not labs else (next(iter(labs)) if len(labs)==1 else "mixed")
T = [r for r in rd(f"{REF}/truth_set.tsv") if r["stated"] in ("MAT1-1","MAT1-2")]
print("\ntruth rows with stated MAT1-1/MAT1-2:", len(T))
tab = {"before": collections.Counter(), "after": collections.Counter()}; ranset = 0; flips = []
for r in T:
    g = r["genome"]
    if g in both: ranset += 1; b = gcall(L, g); a = gcall(C, g)
    else: b = a = v6call(g)
    for nm, c in (("before", b), ("after", a)):
        tab[nm]["agree" if c == r["stated"] else (c if c in ("undetermined","none","mixed") else "wrong")] += 1
    if a != b: flips.append((g, r["organism"], r["stated"], b, a))
print(" in run set:", ranset, "(others: v0.6.0 label, unchanged by construction)")
for nm in ("before","after"): print(" ", nm, dict(tab[nm]))
for f in flips: print("  flip", f)
P = [g for g, r in gen.items() if r["species"].startswith("Parastagonospora")]
print("\nParastagonospora genomes:", len(P), "in run set:", sum(g in both for g in P))
for nm, loci in (("legacy", L), ("core_first", C)):
    cnt = collections.Counter(); rec = 0
    for g in P:
        if g not in both: continue
        ls = [r for k, r in loci.items() if k[0]==g]
        if not ls: cnt["no_record"] += 1; continue
        rec += 1; r = ls[0]; genes = set(r["genes_found"].split("|"))
        cnt[(r["idiomorph"], "both_core_genes" if ("MAT1-2-1" in genes and "MAT1-1-1" in genes) else "one")] += 1
    print(" ", nm, "with records:", rec, dict(cnt))
S = {(r["genome"], r["contig"], r["start"]): r for r in rd(f"{REF}/scored.tsv")}
for nm, loci in (("legacy", L), ("core_first", C)):
    c = collections.Counter()
    for k, r in S.items():
        if k[0] in both and k in loci:
            c["scored"] += 1; c["agree"] += loci[k]["idiomorph"] == r["profile_call"]
    print("profile", nm, dict(c))
c = collections.Counter()
for r in rd(f"{V}/loci.tsv"):
    if r["family_called"]=="Ascomycota:MAT" and r["genome"] in both:
        k=(r["genome"],r["contig"],r["start"]); c["v6"]+=1
        if k in L: c["in_legacy"]+=1; c["same_label"] += L[k]["idiomorph"]==r["idiomorph"]
print("v0.6.0 vs legacy-on-main", dict(c))
