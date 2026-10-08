#!/bin/bash
#SBATCH -p short -A stajichlab -c 1 --mem 4gb -t 2:00:00 -J fola-mixes --array=0-2 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_mixes/out/slurm_%A_%a.log
PAIRS=("VSP-1150 VSP-0798" "VSP-0777 VSP-2032" "VSP-1057 JCP360")
set -- ${PAIRS[$SLURM_ARRAY_TASK_ID]}
export PANEL_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/reads-type/panels/fusarium_oxysporum_fola READS_DIR=/rhome/nicolel/bigdata/Fola/05_Pre-Processing PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict
/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/reads-type/scripts/simulate_mixes.sh $1 $2 /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_mixes/out/mix_$1_$2.tsv
