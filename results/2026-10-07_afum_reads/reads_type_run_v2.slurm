#!/bin/bash
#SBATCH -p short -A stajichlab -c 1 --mem 4gb -t 2:00:00 -J afum-reads-type --array=0-27 --exclude=c01 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/types/out_v2/slurm_%A_%a.log
# matpredict reads-type with the curated A1163 (MAT1-1) and Af293 (MAT1-2) locus segments; first 8M reads of each read set.
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/src
n=$(printf "%02d" $SLURM_ARRAY_TASK_ID)
/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/test/bin/python -m MATPredict reads-type --idiomorph MAT1-1=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/panels/aspergillus_fumigatus/MAT1-1.fasta --idiomorph MAT1-2=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/afum-reads-type/panels/aspergillus_fumigatus/MAT1-2.fasta --samples /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/chunks/s_$n.tsv --max-reads 8000000 --out /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-07_afum_reads/types/out_v2/types_$n.tsv
