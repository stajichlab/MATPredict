#!/bin/bash
#SBATCH -p exfab -A exfab -c 4 --mem 8G -t 01:00:00 -J pf-miniprot
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/pf_miniprot_%j.log
# P. fabianii: curated Saccharomycotina MAT proteins + ONH67106 aligned with miniprot to four assemblies (no RNA-seq exists for the species).
set -euo pipefail
W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-30_petasospora/work; cd $W
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
for f in GCA_*.fna.gz; do n=${f%.fna.gz}; zcat $f > $n.fa; miniprot -t4 -d $n.mpi $n.fa 2>/dev/null; miniprot -t4 --gff --outs=0.5 -P PF $n.mpi refs.faa > $n.miniprot.gff 2> $n.miniprot.log; done
touch ../DONE
