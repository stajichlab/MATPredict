"""Do curated mating receptors form clades that non-mating STE3 copies stay out of?"""
import csv, collections
from Bio import Phylo
hits = {r["seq_id"]: r for r in csv.DictReader(open("ste3_hits.tsv"), delimiter="\t")}
t = Phylo.read("tree_ft.treefile", "newick")
out = [x for x in t.get_terminals() if x.name.startswith("OUT|")]
t.root_with_outgroup(out[0])
leaves = {x.name: x for x in t.get_terminals()}
def grp(n):
    if not n.startswith("REF|"): return None
    g = n.split("|")[2]
    if g in ("bar3", "bbr2", "pheromone_receptor"): return "Agaricales_PR"
    if g == "pra1": return "Ustilaginales_pra"
    if g == "STE3" and "_MAT_" in n: return "Tremellales_STE3"
    if g == "STE3a1": return "Sporidiobolales_A1"
    if g == "STE3a2": return "Sporidiobolales_A2"
    if g in ("STE3", "STE3v2") and "wallMAT" in n: return "Wallemia_STE3"
known = collections.defaultdict(list)
for n in leaves:
    g = grp(n)
    if g: known[g].append(n)
rep = open("mating_clades.tsv", "w")
rep.write("group\tn_known\tmonophyletic\tclade_size\tsupport\tqueries_in_clade\tquery_orders\n")
inclade = {}
for g, ks in sorted(known.items()):
    ca = t.common_ancestor([leaves[k] for k in ks])
    tips = [x.name for x in ca.get_terminals()]
    q = [n for n in tips if not n.startswith(("REF|", "OUT|"))]
    other_ref = [n for n in tips if n.startswith("REF|") and grp(n) != g]
    mono = "yes" if not other_ref else f"no (+{len(other_ref)} other-group refs)"
    orders = collections.Counter(hits[n]["order"] for n in q)
    sup = ca.confidence
    rep.write(f"{g}\t{len(ks)}\t{mono}\t{len(tips)}\t{sup}\t{len(q)}\t{dict(orders.most_common(6))}\n")
    for n in q: inclade.setdefault(n, set()).add(g)
rep.close()
print(open("mating_clades.tsv").read())
# copies per genome and in-clade per genome, by order
per = collections.defaultdict(lambda: [0, 0])
for n, r in hits.items():
    if r["kind"] != "query" or n not in leaves: continue
    per[(r["order"], r["ASMID"])][0] += 1
    if n in inclade: per[(r["order"], r["ASMID"])][1] += 1
agg = collections.defaultdict(list)
for (o, a), (tot, inc) in per.items(): agg[o].append((tot, inc))
print("order\tgenomes\tmedian_copies\tmedian_in_mating_clades\tgenomes_with_>=1_in_clade")
for o, v in sorted(agg.items(), key=lambda kv: -len(kv[1])):
    tots = sorted(x[0] for x in v); incs = sorted(x[1] for x in v)
    print(f"{o}\t{len(v)}\t{tots[len(v)//2]}\t{incs[len(v)//2]}\t{sum(1 for x in incs if x>0)}")
with open("query_clade.tsv", "w") as fo:
    fo.write("seq_id\torder\tspecies\tstrain\tmating_clade\n")
    for n, r in hits.items():
        if r["kind"] == "query" and n in leaves:
            fo.write(f"{n}\t{r['order']}\t{r['species']}\t{r['strain']}\t{','.join(sorted(inclade.get(n, [])))}\n")
