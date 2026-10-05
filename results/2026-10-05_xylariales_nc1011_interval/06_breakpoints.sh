#!/bin/bash
#SBATCH -p batch -c 4 --mem 16G -t 01:00:00 -J xyl-bkpt
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/xyl_bkpt_%j.log
set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval/results/2026-10-05_xylariales_nc1011_interval/common.sh
set +eu; source /etc/profile 2>/dev/null; module load samtools/1.19.2 minimap2/2.30; set -eu
cd $D && /bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/python 06_breakpoints.py $D
