#!/usr/bin/bash
#SBATCH -N 1 -c 1 --mem 2gb --time 00:30:00 -J ss-validate
set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict
T=$R/results/2026-10-30_saccharomycopsis/tree
export PATH="$R/.pixi/envs/default/bin:$PATH"
cd $T
for r in 44080_gdd12_MATss_MATa-MATalpha 44080_gdd12_MATss_MATa 36022_joy008_MATpf_MATa 36022_tbrc-4498_MATpf_MATalpha 36038_pb4207_MATbz_MATa-MATalpha; do
  PYTHONPATH="$T/src" MATPREDICT_DB_ROOT="$T/db" python -m MATPredict curate-db validate --phylum Ascomycota --record-id $r > $T/$r.validate.json
  grep -n "taxonomy_current\|taxonomy_error\|sequence_match\|status" -A1 $T/$r.validate.json | head -12
done
