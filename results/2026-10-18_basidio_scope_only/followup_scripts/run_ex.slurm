#!/usr/bin/bash
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J fam18-eq -p short
# usage: EXCLUDE=id1,id2 sbatch --export=ALL --array=0-9%10 run_ex.slurm LISTPREFIX OUTDIR
set -uo pipefail
W=/bigdata/stajichlab/jstajich/projects/mat-fam18
cd $W
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH
export PYTHONPATH=$W/src
export MATPREDICT_HTML=0
SC=${SCRATCH:-/tmp}/fam18eq_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $SC
L=$(printf "%s.%02d" $1 $SLURM_ARRAY_TASK_ID)
OUT=$2
mkdir -p $OUT
while IFS=$'\t' read -r acc taxid fa; do
  d=$OUT/${taxid}_${acc}
  [ -s $d/detection_report.yaml ] && continue
  mkdir -p $d
  zcat $fa > $SC/g.fa
  python -m MATPredict detect --genome $SC/g.fa --taxid $taxid --out-dir $d --exclude-records "$EXCLUDE" > $d/stdout.log 2> $d/stderr.log || echo "FAIL $acc" >> $OUT/failed.$SLURM_ARRAY_TASK_ID
done < $L
rm -rf $SC
