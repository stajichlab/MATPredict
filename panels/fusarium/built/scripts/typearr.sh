#!/bin/bash
#SBATCH -p short -A stajichlab -c 2 --mem 8G -t 2:00:00 -J mrx_type
# usage: sbatch --array=1-N%K typearr.sh LISTFILE ; each line: RUN PREFIX [PREFIX ...]
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $1)
exec bash /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm/typejob.sh $line
