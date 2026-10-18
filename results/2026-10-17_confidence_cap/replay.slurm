#!/bin/bash -l
#SBATCH -A stajichlab -p short -c 4 --mem 16G -t 02:00:00 -J conf-replay -o /bigdata/stajichlab/jstajich/projects/conf/logs/conf-replay.out
W=/bigdata/stajichlab/jstajich/projects/conf
export PATH=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/default/bin:$PATH
PYTHONPATH=$W/src python $W/results/2026-10-17_confidence_cap/scan_replay.py --db $W/db \
  --reports "/bigdata/stajichlab/jstajich/projects/core-hmms/asco/wave_*/runs/*/detection_report.yaml" \
  --out $W/results/2026-10-17_confidence_cap/loci_features2.tsv --workers 4
