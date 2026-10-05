#!/usr/bin/bash -l
#SBATCH -N 1 -c 8 --mem 32gb --time=0-02:00:00 -p short
#SBATCH --out /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/caax-receptor-test/results/2026-10-05_caax_receptor_test/slurm_%A_%a.log
#SBATCH --array=1-4
# Panel extension: 4 Agaricomycete genomes with a curated PR (B-locus) record.
# Usage: sbatch run_new_genomes.sh   (ROOT is absolute: no BASH_SOURCE on SLURM)
ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/caax-receptor-test/results/2026-10-05_caax_receptor_test
PYENV=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
export PATH=$PYENV:$PATH
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $ROOT/new_genomes.txt)
set -- $LINE
cd $ROOT
$PYENV/python scan_genome.py "$1" $ROOT/out_new "$2"
