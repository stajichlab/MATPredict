#!/usr/bin/bash -l
#SBATCH -N 1 -c 8 --mem 16gb --time=0-02:00:00 -p short
#SBATCH --array=1-62
#SBATCH --out /bigdata/stajichlab/jstajich/agari_work/prot_slurm_%A_%a.log
# extract_loci_proteins.py on every scanned Agaricomycetes genome (scan==True in agari_genomes.tsv), 30 per task.
ROOT=/bigdata/stajichlab/jstajich/agari_work
PYENV=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
export PATH=$PYENV:$PATH
export SCRATCH=$ROOT/tmp; mkdir -p $SCRATCH $ROOT/prot_out
PER=30
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * PER + 1 )); LAST=$(( FIRST + PER - 1 ))
cd $ROOT
awk -F'\t' 'NR>1 && $11=="True"{print $1}' agari_genomes.tsv | sed -n "${FIRST},${LAST}p" | while read asm; do
    [ -e prot_out/$asm.done ] && continue
    $PYENV/python extract_loci_proteins.py "$asm" $ROOT/prot_out && touch prot_out/$asm.done || echo "FAILED $asm"
done
