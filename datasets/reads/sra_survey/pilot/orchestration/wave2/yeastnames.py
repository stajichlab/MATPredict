import csv,collections
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/wave_inputs"
G={r["genome"]:r for r in csv.DictReader(open(W+"/asco_genomes.tsv"),delimiter="\t")}
sp=[l.split("\t")[0] for l in open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2/species.tsv")][1:]
held=set()
import os
P="/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/waves/wave2/plan"
for s in os.listdir(P):
    if open(P+"/"+s+"/status").read().strip()=="held:panel_build_failed": held.add(s)
c=collections.Counter(); fam=collections.Counter(); genes=collections.Counter(); per=collections.defaultdict(collections.Counter)
for r in csv.DictReader(open(W+"/asco_loci.tsv"),delimiter="\t"):
    g=G.get(r["genome"]); 
    if not g: continue
    slug=g["species"].lower().replace(" ","_")
    if slug in held:
        c[r["idiomorph"]]+=1; fam[r["family_called"]]+=1; per[slug][r["idiomorph"]]+=1; genes[r["genes_found"][:60]]+=1
print("idiomorph names in held species:",c.most_common(12)); print("families:",fam.most_common(8)); print("genes:",genes.most_common(8))
print({k:dict(v) for k,v in list(per.items())[:8]})
print("held species with any loci:",len(per),"of",len(held))
