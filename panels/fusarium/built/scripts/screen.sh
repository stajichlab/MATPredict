#!/bin/bash
#SBATCH -p short -A stajichlab -c 8 --mem 16G -t 2:00:00 -J mrx_screen
# usage: sbatch screen.sh RUN   -> download, every-Nth subsample (~2M pairs), diamond blastx vs MAT proteins, counts per gene
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/common.sh
run=$1; mkdir -p $R/res/screen
fetch $run || exit 1
subsample $run 2000000
cd $S/$run
for m in 1 2; do
  diamond blastx -d $R/res/ref/matprot -q sub_$m.fastq.gz --threads 8 -e 1e-5 -k 1 --outfmt 6 qseqid sseqid pident length evalue bitscore -o hits_$m.tsv 2>&1 | tail -1
done
cat hits_1.tsv hits_2.tsv | awk -F'\t' '{split($2,a,"name="); split(a[2],b,"|"); split($2,c,"_MAT_"); print $1"\t"b[1]"\t"$3"\t"$4}' | sort -u -k1,2 > $R/res/screen/$run.hits.tsv
awk -F'\t' '{n[$2]++; p[$2]+=$3} END{for(g in n) printf "%s\t%d\t%.1f\n", g, n[g], p[g]/n[g]}' $R/res/screen/$run.hits.tsv | sort > $R/res/screen/$run.genes.tsv
cat $R/res/screen/$run.genes.tsv
