#!/bin/bash
#SBATCH --partition=short --time=2:00:00 --cpus-per-task=16 --mem=32gb
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-29_aligner_default/logs/%x_%j.log
set -u
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict
WT=$REPO/.claude/worktrees/run-algcand
export PYTHONPATH=$WT/src MATPREDICT_DB_ROOT=$WT/db MATPREDICT_CACHE_DIR=$REPO/.matpredict_cache
PY=$REPO/.pixi/envs/default/bin/python
export PATH=$REPO/.pixi/envs/default/bin:$PATH
OUT=$REPO/results/2026-09-29_aligner_default/lcg_runs
S=${SCRATCH:?}/lcg; mkdir -p $S
echo "src=$WT/src db=$WT/db list=$LIST"
run_one(){ org="$1"; g="$2"; out="$OUT/$org"
  [[ -s "$out/detection_report.yaml" ]] && return 0
  mkdir -p "$out"; fna="$S/$org.fna"; cp "$g" "$fna" || { echo "FAIL $org copy" >&2; return 0; }
  t0=$(date +%s)
  timeout 3600 "$PY" -m MATPredict detect --genome "$fna" --phylum Mucoromycota --out-dir "$out" \
     --evidence-diagnostics "$out/evidence_diagnostics.jsonl" > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; echo "$(( $(date +%s) - t0 ))" > "$out/wall_seconds"
  [[ $rc -ne 0 ]] && echo "FAIL $org rc=$rc: $(tail -1 "$out/stderr.log")" >&2
  rm -f "$fna"; }
export -f run_one; export OUT S PY
awk -F'\t' '{print $1" "$2}' "$LIST" | xargs -P 16 -L1 bash -c 'run_one "$0" "$1"'
echo "done $(find $OUT -name detection_report.yaml -size +0 | wc -l)"
