#!/bin/bash
#SBATCH --partition=short --time=2:00:00 --cpus-per-task=16 --mem=48gb -J mucoro_mat
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/logs/%x_%A_%a.log
#SBATCH --array=0-3
# Mucoromycotina MAT campaign (2026-10-03): BFD Mucorales+Umbelopsidales, LCG, Jena.
# All genomes: detect --phylum Mucoromycota --emit-cds-fasta, no taxid (uniform mode;
# LCG/Jena have no taxid). Code + db: frozen worktree run-7c7ed99 (= PR #10 head).
set -u
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict
WT=$REPO/.claude/worktrees/run-7c7ed99
[[ -f $WT/src/MATPredict/__init__.py && -f $WT/db/Mucoromycota/order.yml ]] || { echo "worktree missing" >&2; exit 1; }
export PYTHONPATH=$WT/src MATPREDICT_DB_ROOT=$WT/db MATPREDICT_CACHE_DIR=$REPO/.matpredict_cache
PY=$REPO/.pixi/envs/default/bin/python
export PATH=$REPO/.pixi/envs/default/bin:$PATH
C=$REPO/results/2026-10-03_mucoromycotina_mat
OUT=$C/runs
LIST=$C/chunk_${SLURM_ARRAY_TASK_ID}.tsv
S=${SCRATCH:?}/mm; mkdir -p $S
echo "src=$WT/src db=$WT/db list=$LIST $(wc -l < $LIST) genomes host=$(hostname)"
run_one(){ id="$1"; g="$2"; out="$OUT/$id"
  [[ -s "$out/detection_report.yaml" ]] && return 0
  mkdir -p "$out"; fna="$S/$id.fna"
  if [[ "$g" == *.gz ]]; then zcat "$g" > "$fna"; else cp "$g" "$fna"; fi || { echo "FAIL $id copy" >&2; return 0; }
  t0=$(date +%s)
  timeout 3600 "$PY" -m MATPredict detect --genome "$fna" --phylum Mucoromycota --out-dir "$out" \
     --emit-cds-fasta --evidence-diagnostics "$out/evidence_diagnostics.jsonl" > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; echo "$(( $(date +%s) - t0 ))" > "$out/wall_seconds"
  [[ $rc -ne 0 ]] && echo "FAIL $id rc=$rc: $(tail -1 "$out/stderr.log")" >&2
  rm -f "$fna"; }
export -f run_one; export OUT S PY
awk -F'\t' '{print $1" "$2}' "$LIST" | xargs -P 16 -L1 bash -c 'run_one "$0" "$1"'
echo "done $(find $OUT -name detection_report.yaml -size +0 | wc -l) total reports"
