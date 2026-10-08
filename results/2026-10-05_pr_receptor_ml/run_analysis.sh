#!/usr/bin/bash -l
#SBATCH -N 1 -c 4 --mem 12gb --time=0-02:00:00 -p short
#SBATCH --array=0-3
#SBATCH --out /bigdata/stajichlab/jstajich/prml_work/slurm_analysis_%A_%a.log
cd /bigdata/stajichlab/jstajich/prml_work
source /bigdata/stajichlab/jstajich/envs/prml/bin/activate
export OMP_NUM_THREADS=2 NBOOT=1000
SETS=(H Hq Hnocc ALL)
ONLY=${SETS[$SLURM_ARRAY_TASK_ID]} python analysis.py
