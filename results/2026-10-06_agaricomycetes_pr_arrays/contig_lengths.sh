#!/usr/bin/bash -l
#SBATCH -N 1 -c 1 --mem 4gb --time=0-01:00:00 -p short
#SBATCH --array=1-62
#SBATCH --out /bigdata/stajichlab/jstajich/agari_work/slurm_clen_%A_%a.log
# contig name + length for every scanned Agaricomycetes genome (30 per task); needed for contig-matched null models.
ROOT=/bigdata/stajichlab/jstajich/agari_work
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
mkdir -p $ROOT/clen
PER=30
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 1 )); LAST=$(( FIRST + PER - 1 ))
cd $ROOT
awk -F'\t' 'NR>1 && $11=="True"{print $1}' agari_genomes.tsv | sed -n "${FIRST},${LAST}p" | while read asm; do
    [ -s clen/$asm.tsv ] && continue
    zcat $LIB/$asm.fa.gz | awk '/^>/{if(n)print name"\t"len; name=substr($1,2); len=0; n=1; next}{len+=length($0)}END{if(n)print name"\t"len}' > clen/$asm.tsv.tmp && mv clen/$asm.tsv.tmp clen/$asm.tsv
done
