#!/bin/bash
#SBATCH --partition=short --time=2:00:00 --cpus-per-task=16 --mem=48gb --job-name=rt-check
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_runtime_check/logs/%x_%j.log
set -u
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict
D=$REPO/results/2026-10-01_runtime_check
PY=$REPO/.pixi/envs/default/bin/python
export PATH=$REPO/.pixi/envs/default/bin:$PATH MATPREDICT_CACHE_DIR=$REPO/.matpredict_cache
S=${SCRATCH:?}/rt; mkdir -p $S
echo "node=$(hostname) start=$(date)"
run_one(){ c="$1"; org="$2"; g="$3"; mode="$4"
  WT=$REPO/.claude/worktrees/run-$c; out=$D/$mode/$c/$org; mkdir -p "$out"
  fna=$S/${c}_${org}.fna; cp "$g" "$fna"
  pre=""; [[ $mode == prof ]] && pre="-m cProfile -o $out/profile.pstats"
  t0=$(date +%s.%N)
  PYTHONPATH=$WT/src MATPREDICT_DB_ROOT=$WT/db timeout 3600 "$PY" $pre -m MATPredict detect --genome "$fna" --phylum Mucoromycota --out-dir "$out" \
     --evidence-diagnostics "$out/evidence_diagnostics.jsonl" > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; t1=$(date +%s.%N)
  echo -e "$c\t$org\t$mode\t$rc\t$(echo "$t1 - $t0" | bc)" >> $D/$mode/timings.tsv
  rm -f "$fna"; }
export -f run_one; export D S PY REPO
mkdir -p $D/time $D/prof
awk -F'\t' '{print $1" "$2" "$3" time"}' $D/tasks.tsv | xargs -P 8 -L1 bash -c 'run_one "$0" "$1" "$2" "$3"'
echo "timing done $(date)"
awk -F'\t' '{print $1" "$2" "$3" prof"}' $D/prof_tasks.tsv | xargs -P 6 -L1 bash -c 'run_one "$0" "$1" "$2" "$3"'
echo "end=$(date)"
