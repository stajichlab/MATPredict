#!/bin/bash -l
P=/bigdata/stajichlab/jstajich/projects/MATPredict; R=$P/results/sra_reads_typing; W=$R/waves/wave2
export PATH=$P/.pixi/envs/default/bin:$PATH
echo "== read results for B05.10 runs"
for r in SRR29713114 SRR29713115 SRR29713116; do head -2 $W/out/$r.botrytis.tsv | tr -d '\r' | cut -f2-11; done | awk 'NR%2==1&&NR>1{next}1'
echo "== detect loci for B05.10 assemblies"
python3 - <<PY
import csv
G=[r for r in csv.DictReader(open("$R/wave_inputs/asco_genomes.tsv"),delimiter="\t") if r["species"].startswith("Botrytis cinerea")]
print(len(G),"Botrytis cinerea genomes; B05.10-like:",[g["genome"] for g in G if "05" in g["genome"] or "B05" in g["genome"]][:5])
names=[g["genome"] for g in G]
L=[r for r in csv.DictReader(open("$R/wave_inputs/asco_loci.tsv"),delimiter="\t") if r["genome"] in names]
import collections
print(collections.Counter(r["idiomorph"] for r in L))
for g in G:
    if "GCF_000143535" in g["genome"]: print(g["genome"],g["status"],[ (r["contig"],r["start"],r["end"],r["idiomorph"],r["locus_class"],r["genes_found"]) for r in L if r["genome"]==g["genome"]])
PY
echo "== panel"
ls $W/plan/botrytis_cinerea/panel; cat $W/plan/botrytis_cinerea/panel/botrytis_lineages.tsv | cut -f1-7
