#!/bin/bash
#SBATCH -p exfab -A exfab -c 8 --mem 32gb -t 8:00:00 -J fola-recruit --array=0-3 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/recruit/slurm_%A_%a.log
# Recruit read pairs with a blastx hit to MAT or flank proteins (any tier), assemble with SPAdes, run detect on the contigs.
set +eu; source /etc/profile; module load spades/4.2.0; set -eu
S=(VSP-0931 50a VSP-0947 VSP-0980); s=${S[$SLURM_ARRAY_TASK_ID]}
PRE=/rhome/nicolel/bigdata/Fola/05_Pre-Processing; O=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/recruit/$s; mkdir -p $O
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
for m in 1 2; do
  diamond blastx -d /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/recruit/recruit_db -q $PRE/${s}_R${m}_trimmed.fastq.gz -o $SCRATCH/h$m.tsv -f 6 qseqid -k 1 -e 1e-5 -p 8 --quiet
done
cat $SCRATCH/h1.tsv $SCRATCH/h2.tsv | sed -E 's#/[12]$##' | sort -u > $O/hit_ids.txt
/usr/bin/python3.12 /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/../../.claude/worktrees/reads-type/scripts/recruit_pairs.py $PRE/${s}_R1_trimmed.fastq.gz $PRE/${s}_R2_trimmed.fastq.gz $O/hit_ids.txt $O/rec_R1.fastq.gz $O/rec_R2.fastq.gz 2> $O/recruit.log
spades.py -1 $O/rec_R1.fastq.gz -2 $O/rec_R2.fastq.gz --careful -t 8 -m 30 -o $SCRATCH/sp > $O/spades.stdout 2>&1
cp $SCRATCH/sp/{contigs.fasta,scaffolds.fasta,spades.log} $O/
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/../../.claude/worktrees/run-94c1a3b/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/../../.claude/worktrees/run-94c1a3b/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/../../.matpredict_cache
python -m MATPredict detect --genome $O/scaffolds.fasta --taxid 5506 --out-dir $O/detect > $O/detect.log 2>&1 || echo "detect failed"
echo done
