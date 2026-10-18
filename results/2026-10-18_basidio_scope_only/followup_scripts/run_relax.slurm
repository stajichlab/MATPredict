#!/usr/bin/bash
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J fam18-relax -p short
# usage: sbatch --array=0-9%2 run_relax.slurm BAR   (genomes of lists/rc.NN, reports in runs_relax/bar<BAR>/)
set -uo pipefail
W=/bigdata/stajichlab/jstajich/projects/mat-fam18
cd $W
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH
export PYTHONPATH=$W/src
export MATPREDICT_HTML=0
BAR=$1
SC=${SCRATCH:-/tmp}/fam18rl_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $SC
L=$(printf "lists/rc.%02d" $SLURM_ARRAY_TASK_ID)
OUT=runs_relax/bar$BAR
mkdir -p $OUT
while IFS=$'\t' read -r acc taxid fa; do
  d=$OUT/${taxid}_${acc}
  [ -s $d/detection_report.yaml ] && continue
  mkdir -p $d
  zcat $fa > $SC/g.fa
  python relax_driver.py $BAR $SC/g.fa $taxid $d > $d/stdout.log 2> $d/stderr.log || echo "FAIL $acc" >> $OUT/failed.$SLURM_ARRAY_TASK_ID
done < $L
rm -rf $SC
