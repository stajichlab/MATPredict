#!/bin/bash
#SBATCH -p short -A stajichlab -c 4 --mem 8G -t 1:00:00 -J mrx_cov
# usage: sbatch covcheck.sh RUN CONTIG.fa   -> map every-Nth-pair subsample to the full recruited contig (AAFTF image bwa+samtools), depth/pairing stats
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/common.sh
run=$1; ctg=$2; mkdir -p $R/res/cov
fetch $run || exit 1
subsample $run 1600000
cd $S/$run; cp $ctg ref.fa
$AP bwa index ref.fa 2>/dev/null
$AP bwa mem -t 4 ref.fa sub_1.fastq.gz sub_2.fastq.gz 2>/dev/null | $AP samtools sort -o m.bam - 
$AP samtools index m.bam
$AP samtools flagstat m.bam > $R/res/cov/$run.flagstat.txt
$AP samtools depth -a m.bam | awk '{n++; s+=$3; if($3==0)z++; if(min==""||$3<min)min=$3} END{printf "bases %d mean_depth %.1f min %d zero_cov_bases %d\n", n, s/n, min, z}' > $R/res/cov/$run.depth.txt
$AP samtools depth -a m.bam | awk '{print $1"\t"$2"\t"$3}' | awk 'NR%50==1' > $R/res/cov/$run.depth_every50.tsv
cat $R/res/cov/$run.depth.txt; grep -E "mapped|properly" $R/res/cov/$run.flagstat.txt
