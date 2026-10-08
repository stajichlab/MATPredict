#!/bin/bash
#SBATCH -p short -A stajichlab -c 4 --mem 16gb -t 2:00:00 -J afum-recruit --array=0-7 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/recruit/slurm_%A_%a.log
# Recruit read pairs onto the two curated A. fumigatus locus segments (minimap2 -ax sr), extend 3 more rounds against the
# growing contigs, assemble each round (SPAdes --careful), then run detect on the final scaffolds.
set +eu; source /etc/profile; module load spades/4.2.0 minimap2/2.30 samtools; set -eu
export PATH=$PATH:/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-2d1860a/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-2d1860a/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
REC=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/scripts/recruit_pairs.py
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/recruit/strains.tsv); s=$(echo "$line" | cut -f1); R1=$(echo "$line" | cut -f2); R2=$(echo "$line" | cut -f3)
O=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/recruit/$s; mkdir -p $O
cat /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/panels/aspergillus_fumigatus/MAT1-1.fasta /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/panels/aspergillus_fumigatus/MAT1-2.fasta > $O/bait.fasta
ref=$O/bait.fasta
echo -e "round\tpairs\tscaffolds\ttotal_bp\tlongest" > $O/rounds.tsv
stats(){ awk '/^>/{if(n)print l; n=1; l=0; next}{l+=length}END{print l}' $1 | sort -nr | awk '{a[NR]=$1; t+=$1}END{printf "%d\t%d\t%d", NR, t, a[1]}'; }
for r in 1 2 3 4; do
  q=10; [ $r -eq 1 ] && q=0
  minimap2 -ax sr -t 4 $ref $R1 $R2 2>/dev/null | samtools view -q $q -F 4 | cut -f1 | sed -E 's#/[12]$##' | sort -u > $O/ids_r$r.txt
  n=$(wc -l < $O/ids_r$r.txt)
  /usr/bin/python3.12 $REC $R1 $R2 $O/ids_r$r.txt $O/r${r}_R1.fastq.gz $O/r${r}_R2.fastq.gz 2> $O/recruit_r$r.log
  spades.py -1 $O/r${r}_R1.fastq.gz -2 $O/r${r}_R2.fastq.gz --careful -t 4 -m 14 -o $SCRATCH/sp$r > $O/spades_r$r.stdout 2>&1
  cp $SCRATCH/sp$r/scaffolds.fasta $O/round${r}_scaffolds.fasta
  echo -e "$r\t$n\t$(stats $O/round${r}_scaffolds.fasta)" >> $O/rounds.tsv
  ref=$O/round${r}_scaffolds.fasta
done
python -m MATPredict detect --genome $O/round4_scaffolds.fasta --taxid 746128 --out-dir $O/detect > $O/detect.log 2>&1 || true
echo done
