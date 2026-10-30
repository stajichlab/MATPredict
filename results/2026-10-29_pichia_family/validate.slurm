#!/usr/bin/bash
#SBATCH -N 1 -c 1 --mem 2gb --time 00:30:00 -J pk-validate
set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict
T=$R/results/2026-10-29_pichia_family/tree
export PATH="$R/.pixi/envs/default/bin:$PATH"
cd $T
for r in 4909_cbs-573_MATpk_MATa 4909_cbs-573_MATpk_MATalpha; do
  PYTHONPATH="$T/src" MATPREDICT_DB_ROOT="$T/db" python -m MATPredict curate-db validate --phylum Ascomycota --record-id $r > $T/$r.validate.json
  grep -n "taxonomy_current\|taxonomy_error\|sequence_match" -A1 $T/$r.validate.json | head -8
done
