#!/usr/bin/bash
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J fam18-rc -p short
# usage: sbatch --array=0-9%2 run_rc.slurm
set -uo pipefail
W=/bigdata/stajichlab/jstajich/projects/mat-fam18
cd $W
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test/bin:$PATH
export SCRATCH=${SCRATCH:-/tmp}
L=$(printf "lists/rc.%02d" $SLURM_ARRAY_TASK_ID)
mkdir -p rc/out
while IFS=$'\t' read -r acc taxid fa; do
  [ -s rc/out/$acc.json ] && continue
  python rust_check.py $acc $taxid $fa rc/out/$acc.json 2> rc/out/$acc.err || echo FAIL $acc >> rc/failed.$SLURM_ARRAY_TASK_ID
done < $L
