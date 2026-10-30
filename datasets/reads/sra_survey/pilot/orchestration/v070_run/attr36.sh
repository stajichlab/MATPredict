#!/usr/bin/env bash
# usage: attr36.sh NAME SRCTREE LIST DETECT_ARGS [ENV=VAL]
P=/bigdata/stajichlab/jstajich/projects/MATPredict
A=$P/results/2026-10-26_v070_submit/attr
name=$1; tree=$2; list=$3; da=$4; shift 4
sbatch --parsable -J attr-$name -p highclock -c 8 --mem=32G -t 3:00:00 -o $P/logs/attr36_${name}_%j.log \
  --export=ALL,"$@",DETECT_ARGS="$da",REPO=$P,SRC=$tree/src,MATPREDICT_DB_ROOT=$tree/db,CLADE=attr36_$name,SAMPLE_LIST=$A/$list,OUT_DIR=$A/$name,GENOME_TIMEOUT=3600 \
  $P/.claude/worktrees/run-b317ece/scripts/run_clade_panel.slurm
