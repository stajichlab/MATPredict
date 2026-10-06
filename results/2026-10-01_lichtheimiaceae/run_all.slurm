#!/bin/bash
#SBATCH -p short -N 1 -c 8 --mem=16G --time=2:00:00 -J licht_all -o logs/licht_%A_%a.log --array=0-7
set -euo pipefail
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_lichtheimiaceae
: "${SCRATCH:?}"
cd $D
/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python $D/worker.py $D/all.tsv $D/out $SLURM_ARRAY_TASK_ID 8 8
