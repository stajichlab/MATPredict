#!/usr/bin/bash -l
#SBATCH -N 1 -c 4 --mem 12gb --time=0-02:00:00 -p short
#SBATCH --out /bigdata/stajichlab/jstajich/prml_work/slurm_extra_%j.log
# Ablation of the MAT-gene neighbourhood features (headline set only).
cd /bigdata/stajichlab/jstajich/prml_work
source /bigdata/stajichlab/jstajich/envs/prml/bin/activate
export OMP_NUM_THREADS=2 NBOOT=1000 EXTRA=_extra ONLY=H
python analysis.py
