#!/bin/bash -l
#SBATCH -A preempt -p preempt -c 2 --mem 16G -t 06:00:00 --requeue
# usage: JOBS=jobs_sample.tsv sbatch --array=1-N%15 -J conf-sample conf_run.slurm
# jobs.tsv columns: id, taxid, stated, organism, fasta. Idempotent per variant; output goes to a temp dir and is
# renamed only when detection_report.yaml exists, so a requeue is harmless. Same code, only the flag differs.
W=/bigdata/stajichlab/jstajich/projects/conf
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/default
export PATH=$ENV/bin:$PATH
IFS=$'\t' read ID TAXID STATED ORG FA < <(sed -n ${SLURM_ARRAY_TASK_ID}p $JOBS)
T=${TMPDIR:-/scratch/$USER}/conf_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p $T
for V in off partial_weak; do
  OUT=$W/runs/$V/$ID
  [ -s $OUT/detection_report.yaml ] && continue
  [ -s $T/genome.fa ] || zcat -f $FA > $T/genome.fa
  rm -rf $OUT.tmp; mkdir -p $OUT.tmp
  MATPREDICT_HTML=0 MATPREDICT_DB_ROOT=$W/db PYTHONPATH=$W/src python -m MATPredict detect --taxid $TAXID \
     --genome $T/genome.fa --out-dir $OUT.tmp --confidence-cap-policy $V > $OUT.tmp/stdout.txt 2>&1
  if [ -s $OUT.tmp/detection_report.yaml ]; then rm -rf $OUT; mv $OUT.tmp $OUT; else echo "FAILED $ID $V" >&2; fi
done
rm -rf $T
