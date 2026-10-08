#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 8gb -t 4:00:00 -J fola-blastx --array=0-14 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/hits/slurm_%A_%a.log
# DIAMOND blastx of the first 4M reads (R1) of each strain against tiered protein references (tiers_v2).
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
N=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
for s in $(tail -n +2 /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/chunks/c_$N.tsv | cut -f1); do
  zcat /rhome/nicolel/bigdata/Fola/05_Pre-Processing/${s}_R1_trimmed.fastq.gz | head -n 16000000 > ${SCRATCH:?}/$s.fq
  for t in T1 T2 T4; do
    diamond blastx -d /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/tiers_v2/$t -q $SCRATCH/$s.fq -o /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/hits/${s}_$t.tsv -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore -k 3 -e 1e-5 -p 2 --quiet
  done
  rm -f $SCRATCH/$s.fq
done
echo done
