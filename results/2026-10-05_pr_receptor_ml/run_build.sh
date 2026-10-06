#!/usr/bin/bash -l
#SBATCH -N 1 -c 8 --mem 24gb --time=0-02:00:00 -p short
#SBATCH --array=1-23%6
#SBATCH --out /bigdata/stajichlab/jstajich/prml_work/slurm_build_%A_%a.log
# One task per genome in genomes.txt. W is absolute: no BASH_SOURCE on SLURM.
W=/bigdata/stajichlab/jstajich/prml_work
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
module load ncbi-blast/2.16.0+
export SCRATCH=$W/tmp
ASM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/genomes.txt)
python $W/build_genome.py $ASM $W/out
