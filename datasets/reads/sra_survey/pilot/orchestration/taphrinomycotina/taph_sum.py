import csv,collections
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates"
r=list(csv.DictReader(open(O+"/candidates.tsv"),delimiter="\t"))
print(len(r),"candidates in",len({x["genome"] for x in r}),"genomes")
print(collections.Counter(x["ref_class"] for x in r))
print("by gene:",collections.Counter(x["ref_gene"] for x in r).most_common(14))
print("in report loci:",collections.Counter(x["report_locus"].split(":")[0] for x in r))
sp=collections.defaultdict(list)
for x in r: sp[x["species"]].append(x)
for s,v in sorted(sp.items(),key=lambda kv:-len(kv[1]))[:14]: print(" ",s,len(v),dict(collections.Counter(x["ref_class"] for x in v)),"median id",sorted(float(x["identity"]) for x in v)[len(v)//2])
print("Neolecta:")
for x in r:
    if x["species"].startswith("Neolecta"): print("  ",{k:x[k] for k in ("contig","start","end","best_ref","ref_class","identity","query_coverage","protein_aa","has_stop","report_locus")})
ng=[g for g in {x["genome"] for x in csv.DictReader(open(O+"/genomes.tsv"),delimiter="\t") for x in [x]} ]
have={x["genome"] for x in r}
none=[x for x in csv.DictReader(open(O+"/genomes.tsv"),delimiter="\t") if x["genome"] not in have]
print("genomes with no candidate:",len(none),collections.Counter(x["species"] for x in none).most_common(8))
