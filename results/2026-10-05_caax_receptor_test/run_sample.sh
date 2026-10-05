#!/usr/bin/bash -l
#SBATCH -N 1 -c 8 --mem 32gb --time=0-02:00:00 -p short
#SBATCH --out /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/caax-receptor-test/results/2026-10-05_caax_receptor_test/slurm_sample_%A_%a.log
#SBATCH --array=1-11
# Chance-model scan: scan_genome.py on each sampled genome, 46 genomes per array task.
# ROOT is absolute: no BASH_SOURCE on SLURM.
ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/caax-receptor-test/results/2026-10-05_caax_receptor_test
PYENV=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
export PATH=$PYENV:$PATH
PER=46
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 2 ))
LAST=$(( FIRST + PER - 1 ))
cd $ROOT
mkdir -p out_sample
sed -n "${FIRST},${LAST}p" sample_genomes.tsv | cut -f1 | while read asm; do
    if [ -s out_sample/$asm.loci.tsv ]; then continue; fi
    $PYENV/python scan_genome.py "$asm" $ROOT/out_sample || echo "FAILED $asm"
done
