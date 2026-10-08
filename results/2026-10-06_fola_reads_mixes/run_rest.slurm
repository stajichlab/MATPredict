#!/bin/bash
#SBATCH -p short -A stajichlab -c 1 --mem 4gb -t 2:00:00 -J fola-mixes2 --array=0-1 --exclude=c01 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_mixes/out/slurm2_%A_%a.log
# Remaining mix fractions after the first run stalled on node c01. Pair 2 stops at 80%: VSP-2032 has only ~7.0M R1 reads.
PAIRS=("VSP-0777 VSP-2032 20,50,80" "VSP-1057 JCP360 20,50,80,90,95,98,99,100")
set -- ${PAIRS[$SLURM_ARRAY_TASK_ID]}
export FRACS="${3//,/ }"
export PANEL_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/reads-type/panels/fusarium_oxysporum_fola READS_DIR=/rhome/nicolel/bigdata/Fola/05_Pre-Processing PROJ_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict
/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/reads-type/scripts/simulate_mixes.sh $1 $2 /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_reads_mixes/out/mixrest_$1_$2.tsv
