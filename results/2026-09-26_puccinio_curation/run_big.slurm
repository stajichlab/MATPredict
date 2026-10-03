#!/usr/bin/bash
#SBATCH --partition=epyc --time=10:00:00 --cpus-per-task=16 --mem=48gb --job-name=puccinio-big
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/puccinio_big_%j.log
# Large (>=750 Mb) Pucciniales genomes, after arm, 4 h per-genome limit.
set -u
: "${SCRATCH:?}"
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default
export PATH="$E/bin:$PATH" PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-293640d/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-293640d/db
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
OUT=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_puccinio_curation/after_big
mkdir -p "$OUT" "$SCRATCH/fasta"
run_one() { a="$1"; t="$2"; o="$OUT/runs/$a"; [[ -s "$o/detection_report.yaml" ]] && return 0
  mkdir -p "$o"; f="$SCRATCH/fasta/$a.fna"; zcat "$LIB/$a.fa.gz" > "$f"
  t0=$(date +%s)
  timeout 14400 python -m MATPredict detect --genome "$f" --taxid "$t" --out-dir "$o" \
     --evidence-diagnostics "$o/evidence_diagnostics.jsonl" > "$o/stdout.log" 2> "$o/stderr.log"
  rc=$?; echo "$(( $(date +%s) - t0 ))" > "$o/wall_seconds"
  [[ $rc -ne 0 ]] && echo "FAIL $a rc=$rc: $(tail -1 $o/stderr.log)" >&2; rm -f "$f"; }
export -f run_one; export OUT LIB SCRATCH
xargs -P 4 -L1 bash -c 'run_one "$0" "$1"' < /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_puccinio_curation/pilot_big.tsv
echo done
