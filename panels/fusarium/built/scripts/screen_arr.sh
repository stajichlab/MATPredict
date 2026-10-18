#!/bin/bash
#SBATCH -p short -A stajichlab -c 8 --mem 16G -t 2:00:00 -J mrx_screen
# usage: sbatch --array=1-N%K screen_arr.sh LISTFILE  (one run accession per line)
run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $1)
exec bash /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/screen.sh $run
