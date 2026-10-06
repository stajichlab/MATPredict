#!/usr/bin/bash -l
#SBATCH -N 1 -c 8 --mem 32gb --time=0-02:00:00 -p short
#SBATCH --array=1-62
#SBATCH --out /bigdata/stajichlab/jstajich/agari_work/slurm_%A_%a.log
# array_scan.py on every Agaricomycetes genome of the v0.6.0 run (scan==True in agari_genomes.tsv), 30 per task.
ROOT=/bigdata/stajichlab/jstajich/agari_work
PYENV=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
export PATH=$PYENV:$PATH
export SCRATCH=$ROOT/tmp; mkdir -p $SCRATCH $ROOT/out
PER=30
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 1 )); LAST=$(( FIRST + PER - 1 ))
cd $ROOT
awk -F'\t' 'NR>1 && $11=="True"{print $1}' agari_genomes.tsv | sed -n "${FIRST},${LAST}p" | while read asm; do
    [ -s out/$asm.chance.tsv ] && continue
    $PYENV/python array_scan.py "$asm" $ROOT/out || echo "FAILED $asm"
done
