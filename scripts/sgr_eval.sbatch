#!/bin/bash
#SBATCH -p batch -c 2 --mem 16G -t 12:00:00
# usage: sbatch -J sgr-X [--array=0-19%20] sgr_eval.sbatch LISTPREFIX OUT   (chunk files LISTPREFIX.NN)
cd /bigdata/stajichlab/jstajich/projects/sgr
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/sgr/src
export SCRATCH=${SCRATCH:-/tmp}
L=$1
[ -n "$SLURM_ARRAY_TASK_ID" ] && L=$(printf "%s.%02d" $1 $SLURM_ARRAY_TASK_ID)
cd /bigdata/stajichlab/jstajich/projects/sgr
/rhome/jstajich/.pixi/bin/pixi run --manifest-path /bigdata/stajichlab/jstajich/projects/mat-prof/pixi.toml -e test python scripts/sgr_eval_driver.py $L $2
