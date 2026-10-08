#!/bin/bash
#SBATCH -p exfab -A exfab -c 4 --mem 16gb -t 2:00:00 -J cs-tally --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_core_span_panel/tally_%j.log
set +eu; source /etc/profile; set -eu
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_core_span_panel
/usr/bin/python3.12 tally.py newdb/runs tally_newdb.tsv 10000
