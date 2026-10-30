#!/bin/bash
#SBATCH -p exfab -A exfab -c 8 --mem 32G -t 08:00:00 -J ss-models
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/ss_models_%j.log
# S. schoenii (GCA_054953955.1, strain GDD12): public RNA-seq (PRJNA1364152) mapped to the assembly, and the curated
# Saccharomycotina MAT proteins aligned with miniprot, to build gene models for the MAT loci (the paper's GTF is not public).
set -euo pipefail
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-29_ss_records
W=$D/work; mkdir -p $W; cd $W
set +eu; source /etc/profile 2>/dev/null; module load hisat2/2.2.1 samtools/1.19.2; set -eu
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
ASM=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/054/953/955/GCA_054953955.1_ASM5495395v1/GCA_054953955.1_ASM5495395v1_genomic.fna.gz
[ -s genome.fa ] || curl -sL $ASM | zcat > genome.fa
samtools faidx genome.fa
miniprot -t8 -d genome.mpi genome.fa
miniprot -t8 --gff -P SS genome.mpi $D/refs.faa > miniprot_all.gff 2> miniprot.log
[ -s idx.1.ht2 ] || hisat2-build -p 8 -q genome.fa idx > /dev/null
for RUN in SRR37453601 SRR37453602 SRR37453603; do
  P=${RUN:0:6}; L=00${RUN: -1}
  for m in 1 2; do curl -sL --retry 5 --retry-all-errors -C - ftp.sra.ebi.ac.uk/vol1/fastq/$P/$L/$RUN/${RUN}_$m.fastq.gz -o ${RUN}_$m.fastq.gz; done
  [ -s $RUN.bam ] || hisat2 --dta -p 8 --max-intronlen 5000 --summary-file $RUN.hisat2.txt -x idx -1 ${RUN}_1.fastq.gz -2 ${RUN}_2.fastq.gz | samtools sort -@4 -o $RUN.bam -
  samtools index $RUN.bam
done
samtools merge -f -@4 all.bam SRR3745360?.bam; samtools index all.bam
# per-base depth and spliced-read junctions over every contig carrying an MAT call
for R in CM143126.1 CM143127.1 CM143128.1 CM143129.1 CM143130.1 CM143131.1; do
  samtools depth -a -r $R all.bam | awk -v r=$R '$3>0' > depth_$R.tsv
  samtools view all.bam $R | awk '$6 ~ /N/{print $3,$4,$6}' > spliced_$R.txt
done
touch $D/DONE
