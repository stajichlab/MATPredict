#!/bin/bash
#SBATCH -p exfab -A exfab -c 8 --mem 24G -t 04:00:00 -J ss-stringtie
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/ss_stringtie_%j.log
# Transcript models (StringTie) and ORFs (TransDecoder) from the merged S. schoenii RNA-seq alignment (non-stranded library).
set -euo pipefail
W=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-29_ss_records/work; cd $W
set +eu; source /etc/profile 2>/dev/null; module load stringtie/3.0.3 transdecoder/5.5.0 gffread/0.12.7; set -eu
stringtie all.bam -p 8 -o stringtie.gtf -l SS
gffread -g genome.fa -w transcripts.fa stringtie.gtf
gffread stringtie.gtf -T -o /dev/null 2>/dev/null || true
rm -rf td; mkdir td; cd td
TransDecoder.LongOrfs -t ../transcripts.fa -m 60 > /dev/null
TransDecoder.Predict -t ../transcripts.fa --single_best_only > /dev/null
cp transcripts.fa.transdecoder.pep transcripts.fa.transdecoder.gff3 .. 2>/dev/null || cp *.transdecoder.* ..
touch $W/../DONE2
