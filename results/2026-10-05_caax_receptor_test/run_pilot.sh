#!/usr/bin/bash -l
#SBATCH -N 1 -c 8 --mem 16gb --time=0-01:00:00 -p short
#SBATCH --out /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/caax-receptor-test/results/2026-10-05_caax_receptor_test/slurm_pilot_%j.log
# Receptor-protein pilot: extract proteins, then leave-one-species-out tests.
# ROOT is absolute: no BASH_SOURCE on SLURM.
ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/caax-receptor-test/results/2026-10-05_caax_receptor_test
PYENV=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
module load hmmer/3.4 mafft/7.505
export PATH=$PYENV:$PATH
cd $ROOT
python pilot_extract.py && python pilot_test.py
