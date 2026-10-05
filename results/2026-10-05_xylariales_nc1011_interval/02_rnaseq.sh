#!/bin/bash
#SBATCH -p batch -c 8 --mem 32G -t 12:00:00 -J xyl-rnaseq
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/xyl_rnaseq_%j.log
# Map NC1011 RNA-seq (SRR8861595) to the NC1011 assembly; depth over the SLA2-APN2 block.
set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval/results/2026-10-05_xylariales_nc1011_interval/common.sh
set +eu; source /etc/profile 2>/dev/null; module load sratoolkit/3.2.0 hisat2/2.2.1 samtools/1.19.2; set -eu
RUN=${RUN:-SRR8861595}
cd "$D"; mkdir -p work/rna; W=$D/work/rna
cd "$SCRATCH"
zcat "$GENOME_GZ" > genome.fa
hisat2-build -p 8 -q genome.fa idx > /dev/null
prefetch --max-size 100G -O sra "$RUN"
fasterq-dump -e 8 -t tmp --split-files -O fq "sra/$RUN/$RUN.sra" 2> /dev/null || fasterq-dump -e 8 -t tmp --split-files -O fq "$RUN"
ls -l fq
if [[ -f fq/${RUN}_2.fastq ]]; then IN="-1 fq/${RUN}_1.fastq -2 fq/${RUN}_2.fastq"; else IN="-U fq/${RUN}.fastq"; fi
echo "layout: $IN"
hisat2 --dta -p 8 --max-intronlen 5000 --summary-file "$W/hisat2_summary.txt" -x idx $IN \
  | samtools sort -@ 4 -o "$W/$RUN.bam" -
samtools index "$W/$RUN.bam"
REG="$CONTIG:$REG_START-$REG_END"
samtools depth -a -r "$REG" "$W/$RUN.bam" > "$W/${RUN}_region_depth.tsv"
# spliced reads (CIGAR N) in the region
samtools view "$W/$RUN.bam" "$CONTIG:247000-250500" | awk '$6 ~ /N/' | cut -f1-9 > "$W/${RUN}_spliced_reads_247000-250500.tsv"
samtools flagstat "$W/$RUN.bam" > "$W/${RUN}_flagstat.txt"
cat "$W/hisat2_summary.txt"; wc -l "$W/${RUN}_spliced_reads_247000-250500.tsv"
