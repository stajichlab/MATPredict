import csv, sys, collections
D, REF, V = sys.argv[1:4]
B="resolution_only"
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
sets = collections.defaultdict(set)
for r in csv.reader(open(f"{D}/work2.tsv"), delimiter="\t"):
    for s in r[2].split(","): sets[s].add(r[0])
gen = {r["genome"]: r for r in rd(f"{V}/genomes.tsv")}
def load(pol):
    loci = {}; genomes = set()
    for r in rd(f"{D}/loci2_{pol}.tsv"):
        if r["family"] == "__genome__": genomes.add(r["genome"]); continue
        if r["family"] != "Ascomycota:MAT": continue
        loci[(r["genome"], r["contig"], r["start"])] = r
    return loci, genomes
L, gl = load("legacy"); R, gr = load(B)
both = gl & gr
print("genomes with both reports:", len(both))
names={"set1":"473 changed under core_first","set2":"497 controls","set3":"203 truth","set4":"~300 no-locus-in-v0.6.0"}
for s in ("set1","set2","set3","set4"):
    gs = sets[s] & both
    trans = collections.Counter(); conf = collections.Counter(); chg=set(); n=collections.Counter()
    bygen=collections.Counter(); ident=0
    for k in set(L)|set(R):
        if k[0] not in gs: continue
        a, b = L.get(k), R.get(k); n["legacy"]+=a is not None; n["res"]+=b is not None
        la = a["idiomorph"] if a else "absent"; lb = b["idiomorph"] if b else "absent"
        if la != lb: trans[(la, lb)] += 1; chg.add(k[0]); bygen[(gen.get(k[0],{}).get("species") or "?").split()[0]]+=1
        if a and b:
            if a["confidence"]!=b["confidence"]: conf[(a["confidence"],b["confidence"])]+=1; chg.add(k[0])
            if a["locus_class"]!=b["locus_class"]: conf[("class:"+a["locus_class"],"class:"+b["locus_class"])]+=1; chg.add(k[0])
    print(f"\n== {s} {names[s]}: {len(gs)} run; loci legacy {n['legacy']} res {n['res']}; genomes with label/loci change {sum(1 for g in chg if any(((L.get(k) or {}).get('idiomorph'))!=((R.get(k) or {}).get('idiomorph')) for k in set(L)|set(R) if k[0]==g))}; any change {len(chg)}")
    for t, v in trans.most_common(): print("  label", t, v)
    for t, v in conf.most_common(): print("  conf/class", t, v)
    print("  genera:", bygen.most_common(12))
def gcall(loci, g):
    labs = {r["idiomorph"] for k, r in loci.items() if k[0]==g}
    return "none" if not labs else (next(iter(labs)) if len(labs)==1 else "mixed")
T = [r for r in rd(f"{REF}/truth_set.tsv") if r["stated"] in ("MAT1-1","MAT1-2")]
tab = {"legacy": collections.Counter(), B: collections.Counter()}; flips=[]
for r in T:
    g = r["genome"]
    if g not in both: tab["legacy"]["not_run"]+=1; tab[B]["not_run"]+=1; continue
    a, b = gcall(L, g), gcall(R, g)
    for nm, c in (("legacy", a), (B, b)): tab[nm]["agree" if c==r["stated"] else (c if c in("undetermined","none","mixed") else "wrong")]+=1
    if a!=b: flips.append((g, r["organism"], r["stated"], a, b))
print("\ntruth:", {k: dict(v) for k,v in tab.items()})
for f in flips: print("  flip", f)
P=[g for g,r in gen.items() if r["species"].startswith("Parastagonospora") and g in both]
print("\nParastagonospora run:", len(P))
for nm, loci in (("legacy", L), (B, R)):
    c=collections.Counter()
    for g in P:
        ls=[r for k,r in loci.items() if k[0]==g]
        if not ls: c["no_record"]+=1; continue
        genes=set(ls[0]["genes_found"].split("|"))
        c[(ls[0]["idiomorph"], "both" if {"MAT1-1-1","MAT1-2-1"}<=genes else "one")]+=1
    print(" ", nm, dict(c))
S={(r["genome"],r["contig"],r["start"]):r for r in rd(f"{REF}/scored.tsv")}
for nm, loci in (("legacy", L), (B, R)):
    c=collections.Counter()
    for k,r in S.items():
        if k[0] in both and k in loci: c["scored"]+=1; c["agree"]+=loci[k]["idiomorph"]==r["profile_call"]
    print("profile", nm, dict(c))
# Pyricularia-type: gained loci with MAT1-2-1+MAT1-1-3 genes
gain=[(k,R[k]) for k in R if k not in L and k[0] in both]
print("\nloci gained overall:", len(gain), "with MAT1-2-1+MAT1-1-3 genes:", sum(1 for k,r in gain if {"MAT1-2-1","MAT1-1-3"}<=set(r["genes_found"].split("|"))))
print("gained by genus:", collections.Counter((gen.get(k[0],{}).get("species") or "?").split()[0] for k,r in gain).most_common(10))
lost=[k for k in L if k not in R and k[0] in both]; print("loci lost:", len(lost))
