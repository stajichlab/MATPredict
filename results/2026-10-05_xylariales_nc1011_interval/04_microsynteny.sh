#!/bin/bash
#SBATCH -p batch -c 4 --mem 8G -t 02:00:00 -J xyl-synteny
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/xyl_synteny_%j.log
set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/xylariales-interval/results/2026-10-05_xylariales_nc1011_interval/common.sh
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd $D && python 04_microsynteny.py $D
