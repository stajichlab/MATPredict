import csv,collections,os
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates"
def fa(p):
    d={};n=None
    for l in open(p):
        if l.startswith(">"): n=l[1:].split()[0]; d[n]=[]
        elif n: d[n].append(l.strip())
    return {k:"".join(v) for k,v in d.items()}
for tag,dr in (("strict (id>=25%, cov>=40%)",O),("loose (id>=20%, cov>=20%)",O+"/loose")):
    if not os.path.exists(dr+"/candidates.tsv") or not os.path.exists(dr+"/run.log") or " candidates" not in open(dr+"/run.log").read().splitlines()[-1]: print(tag,"not finished"); continue
    r=list(csv.DictReader(open(dr+"/candidates.tsv"),delimiter="\t")); seq=fa(dr+"/candidates.faa")
    uniq=set(seq.values()); sp=collections.defaultdict(set)
    for x in r: sp[(x["species"],x["ref_gene"])].add(seq[x["candidate"]])
    species={x["species"] for x in r}
    one_per_sp_gene=sum(1 for v in sp.values())
    nonpombe=[x for x in r if not x["species"].startswith("Schizosaccharomyces pombe")]
    print(tag,":",len(r),"candidate proteins in",len({x["genome"] for x in r}),"genomes,",len(species),"species;",len(uniq),"distinct sequences;",one_per_sp_gene,"species x gene groups;",len(nonpombe),"outside S. pombe")
    print("   classes",dict(collections.Counter(x["ref_class"] for x in r)),"| with stop codon",sum(int(x["has_stop"]) for x in r),"| length<80 aa",sum(1 for x in r if int(x["protein_aa"])<80))
    print("   genes",dict(collections.Counter(x["ref_gene"] for x in r).most_common(12)))
    print("   report locus:",dict(collections.Counter(x["report_locus"].split(":")[0] for x in r)))
