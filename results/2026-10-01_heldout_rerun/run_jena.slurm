#!/usr/bin/bash
#SBATCH -p short -N 1 -n 1 -c 16 --mem 32gb --time 2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_heldout_rerun/logs/%x_%j.log
set -u
W=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-52b3ff9
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default
O=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_heldout_rerun
export PATH="$E/bin:$PATH" PYTHONPATH="$W/src" MATPREDICT_DB_ROOT="$W/db" MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache W E O INPUT
S="${SCRATCH:?}"; export S
col=2; [ "$INPUT" = contigs ] && col=3
run_one(){ strain="$1"; fa="$2"
  [ -z "$fa" ] && { echo "SKIP $strain no $INPUT file" >&2; return 0; }
  out="$O/jena_$INPUT/$strain"; [ -s "$out/detection_report.yaml" ] && return 0
  mkdir -p "$out"; cp "$fa" "$S/$strain.$INPUT.fna"
  start=$(date +%s)
  timeout 3600 python -m MATPredict detect --genome "$S/$strain.$INPUT.fna" --phylum Mucoromycota \
     --out-dir "$out" --evidence-diagnostics "$out/evidence_diagnostics.jsonl" > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; echo $(( $(date +%s) - start )) > "$out/wall_seconds"
  rm -f "$S/$strain.$INPUT.fna"
  [ $rc -ne 0 ] && echo "FAIL $strain rc=$rc: $(tail -1 $out/stderr.log)" >&2; return 0; }
export -f run_one
awk -F'\t' -v c=$col '{print $1"\t"$c}' "$O/jena_inputs.tsv" | xargs -P 8 -L1 -d '\n' bash -c 'IFS=$'"'"'\t'"'"' read s f <<< "$0"; run_one "$s" "$f"'
echo "done $INPUT: $(ls $O/jena_$INPUT/*/detection_report.yaml 2>/dev/null | wc -l)"
