#!/bin/bash
#SBATCH -p batch -c 16 --mem 24G -t 01:30:00 -J dothi-sla2
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/dothi_sla2_%j.log
set -euo pipefail
W=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval
D=$W/results/2026-10-05_dothideo_sla2_test
cd $D
awk '/^>/{keep=($1==">KAI0195599.1"||$1==">KAI0195601.1"||$1==">KAI0195602.1"||$1==">KAI0195603.1")} keep' $W/results/2026-10-05_xylariales_nc1011_interval/synteny/q14.faa > queries.faa
grep -c '>' queries.faa
/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python sla2_distance.py $W/results/2026-10-03_ascomycota_v060 queries.faa sla2_distance.tsv 8
