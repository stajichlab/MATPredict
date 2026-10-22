#!/usr/bin/bash -l
# propose -> validate -> accept -> build-gff for every pilot record, then build-duckdb. Compute node.
set -uo pipefail
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/curate-pilot-test
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
P=$HOME/curate_pilot
export PATH=$ENV/bin:$PATH PYTHONPATH=$REPO/src MATPREDICT_DB_ROOT=$REPO/db PYTHONNOUSERSITE=1
export MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache MATPREDICT_HTML=0
: > $P/empty_config.toml; export MATPREDICT_CONFIG=$P/empty_config.toml   # no e-mail is sent to NCBI
cd $REPO
while IFS=$'\t' read -r order rid; do
  echo "=== $order $rid"
  python -m MATPredict curate-db propose --phylum Ascomycota --record-file $P/pilot_out/Ascomycota/$order/$rid/metadata.yaml || { echo "PROPOSE FAIL $rid"; continue; }
  python -m MATPredict curate-db validate --phylum Ascomycota --record-id $rid > $P/validate_$rid.json 2>&1 || echo "VALIDATE FAIL $rid"
  python - "$REPO/db/candidates/Ascomycota/$rid/metadata.yaml" <<'PY'
import sys, yaml
r = yaml.safe_load(open(sys.argv[1]))
v = r["validation"]
print("validation:", v["status"], v.get("accession_resolved"), v.get("accession_resolved_version"), v["sequence_match"]["status"], v["sequence_match"]["per_gene"], "taxonomy_current", v.get("taxonomy_current"), v.get("taxonomy_error"))
PY
  python -m MATPredict curate-db accept --phylum Ascomycota --order-or-family $order --record-id $rid || { echo "ACCEPT FAIL $rid"; continue; }
  python -m MATPredict curate-db build-gff --phylum Ascomycota --order-or-family $order --record-id $rid || echo "BUILDGFF FAIL $rid"
done < $P/pilot_records.tsv
python -m MATPredict curate-db build-duckdb --out $P/pilot.duckdb && echo DUCKDB_OK
