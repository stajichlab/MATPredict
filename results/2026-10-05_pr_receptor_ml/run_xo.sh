#!/usr/bin/bash -l
#SBATCH -N 1 -c 8 --mem 16gb --time=0-01:00:00 -p short
#SBATCH --array=1-23%8
#SBATCH --out /bigdata/stajichlab/jstajich/prml_work/slurm_xo_%A_%a.log
W=/bigdata/stajichlab/jstajich/prml_work
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
ASM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/genomes.txt)
module load ncbi-blast/2.16.0+
python $W/recompute_xo.py $ASM $W/out
