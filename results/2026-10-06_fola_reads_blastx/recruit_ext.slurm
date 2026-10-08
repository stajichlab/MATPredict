#!/bin/bash
#SBATCH -p exfab -A exfab -c 8 --mem 32gb -t 12:00:00 -J fola-recruit-ext --array=0-3 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/recruit_ext/slurm_%A_%a.log
# Iterative extension of the round-1 MAT/flank scaffolds: map all pairs (minimap2, MAPQ>=10), keep pairs with a mapped mate,
# reassemble (SPAdes --careful), repeat. Stops after 4 rounds, if recruited pairs grow < 5%, or if > 1,000,000 pairs (guard).
set +eu; source /etc/profile; module load spades/4.2.0 minimap2/2.30 samtools; set -eu
S=(VSP-0931 50a VSP-0947 VSP-0980); s=${S[$SLURM_ARRAY_TASK_ID]}
PRE=/rhome/nicolel/bigdata/Fola/05_Pre-Processing; O=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/recruit_ext/$s; mkdir -p $O
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
REC=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/reads-type/scripts/recruit_pairs.py
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-94c1a3b/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-94c1a3b/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
cp /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/recruit/$s/scaffolds.fasta $O/round1_scaffolds.fasta
echo -e "round\tpairs\tscaffolds\ttotal_bp\tlongest\tdetect" > $O/rounds.tsv
stats(){ awk '/^>/{if(n)print l; n=1; l=0; next}{l+=length}END{print l}' $1 | sort -nr | awk '{a[NR]=$1; t+=$1}END{printf "%d\t%d\t%d", NR, t, a[1]}'; }
detect_line(){ /usr/bin/python3.12 - $1 <<'PY'
import sys,yaml,os
f=sys.argv[1]+'/detection_report.yaml'
if not os.path.exists(f): print('no_report'); sys.exit()
r=yaml.safe_load(open(f)); d=r.get('detected') or []
print(';'.join(f"{x['idiomorph']}:{x['confidence']}:{x['contig'].split('_length')[0]}:{x['end']-x['start']+1}" for x in d) or 'none')
PY
}
prev=$(wc -l < /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/recruit/$s/hit_ids.txt)
echo -e "1\t$prev\t$(stats $O/round1_scaffolds.fasta)\t$(detect_line /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_blastx/recruit/$s/detect)" >> $O/rounds.tsv
ref=$O/round1_scaffolds.fasta
for r in 2 3 4 5; do
  minimap2 -ax sr -t 8 $ref $PRE/${s}_R1_trimmed.fastq.gz $PRE/${s}_R2_trimmed.fastq.gz 2>/dev/null | samtools view -q 10 -F 4 | cut -f1 | sed -E 's#/[12]$##' | sort -u > $O/ids_r$r.txt
  n=$(wc -l < $O/ids_r$r.txt)
  if [ "$n" -gt 1000000 ]; then echo -e "$r\t$n\tGUARD_stop_over_1M_pairs" >> $O/rounds.tsv; break; fi
  /usr/bin/python3.12 $REC $PRE/${s}_R1_trimmed.fastq.gz $PRE/${s}_R2_trimmed.fastq.gz $O/ids_r$r.txt $O/r${r}_R1.fastq.gz $O/r${r}_R2.fastq.gz 2> $O/recruit_r$r.log
  spades.py -1 $O/r${r}_R1.fastq.gz -2 $O/r${r}_R2.fastq.gz --careful -t 8 -m 30 -o $SCRATCH/sp$r > $O/spades_r$r.stdout 2>&1
  cp $SCRATCH/sp$r/scaffolds.fasta $O/round${r}_scaffolds.fasta
  python -m MATPredict detect --genome $O/round${r}_scaffolds.fasta --taxid 5506 --out-dir $O/detect_r$r > $O/detect_r$r.log 2>&1 || true
  echo -e "$r\t$n\t$(stats $O/round${r}_scaffolds.fasta)\t$(detect_line $O/detect_r$r)" >> $O/rounds.tsv
  growth=$(awk -v a=$n -v b=$prev 'BEGIN{printf "%.3f", (a-b)/b}')
  ref=$O/round${r}_scaffolds.fasta
  if awk -v g=$growth 'BEGIN{exit !(g<0.05)}'; then echo "converged at round $r (growth $growth)" >> $O/rounds.tsv; break; fi
  prev=$n
done
echo done
