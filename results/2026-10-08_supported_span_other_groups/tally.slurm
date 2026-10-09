#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 8gb -t 1:00:00 -J og-tally --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_supported_span_other_groups/tally_%j.log
set +eu; source /etc/profile; set -eu
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_supported_span_other_groups
/usr/bin/python3.12 tally.py floor33/runs tally_beyond_supported_floor33.tsv 10000 supported_span
