#!/bin/bash
#SBATCH -p exfab -A exfab -c 1 --mem 4gb -t 3:00:00 -J fola-reads-type --array=0-14 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_type/out_v2/slurm_%A_%a.log
# matpredict reads-type (worktree reads-type) on the 148 Fola strains; panel = AB011379.2 (MAT1-1), AB011378.1 (MAT1-2).
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/reads-type/src
REF=/rhome/nicolel/bigdata/Fola/06_Align/Mating_Types/References
N=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/test/bin/python -m MATPredict reads-type --idiomorph MAT1-1=$REF/MAT-1.fasta --idiomorph MAT1-2=$REF/MAT-2.fasta --samples /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_type/chunks/c_$N.tsv --max-reads 8000000 --out /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_type/out_v2/types_$N.tsv
