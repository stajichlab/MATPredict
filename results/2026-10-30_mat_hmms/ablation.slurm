#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 8G -t 03:00:00 -J mat-abl
#SBATCH -o /bigdata/stajichlab/jstajich/projects/MATPredict/logs/mat_abl_%j.log
set -euo pipefail
set +eu; source /etc/profile 2>/dev/null; module load mafft/7.505; set -eu; unset MAFFT_BINARIES
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-30_mat_hmms
python ablation.py
touch LOCO_DONE
