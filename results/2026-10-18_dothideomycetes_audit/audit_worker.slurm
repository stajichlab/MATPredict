#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J doth-audit
#SBATCH -p short,highclock,exfab
# One array task = one chunk file (chunks/chunk_NNN.tsv: ASMID<TAB>TAXID). Runs `matpredict detect` on current main
# (worktree run-c9c0ec6 = origin/main c9c0ec6) with --partial-core-policy coverage_profile, all other defaults,
# lineage routing from the taxid. Resumable: a genome with a report is skipped.
set -uo pipefail
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-c9c0ec6
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
OUT=$REPO/results/2026-10-18_dothideomycetes_audit
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=$ENV/bin:$PATH
export PYTHONPATH=$REPO/src
export MATPREDICT_DB_ROOT=$REPO/db
export MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
export MATPREDICT_HTML=0
export PYTHONNOUSERSITE=1
N=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
LIST=$OUT/chunks/chunk_$N.tsv
W=${SCRATCH:?}/doth_audit_${SLURM_JOB_ID}_$N
mkdir -p "$W" "$OUT/runs"
where=$(python -c 'import MATPredict,sys;print(MATPredict.__file__)')
case "$where" in "$REPO"/src/*) ;; *) echo "WRONG CODE: $where" >&2; exit 1;; esac
echo "code=$where host=$(hostname) chunk=$N"
while IFS=$'\t' read -r asm taxid; do
  out=$OUT/runs/$asm
  [[ -s $out/detection_report.yaml ]] && continue
  mkdir -p "$out"
  fna=$W/$asm.fna
  zcat "$LIB/$asm.fa.gz" > "$fna" || { echo "FAIL $asm decompress"; continue; }
  t0=$(date +%s)
  timeout 3600 python -m MATPredict detect --genome "$fna" --taxid "$taxid" --out-dir "$out" \
     --partial-core-policy coverage_profile \
     --evidence-diagnostics "$out/evidence_diagnostics.jsonl" > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?
  echo $(( $(date +%s) - t0 )) > "$out/wall_seconds"
  [[ $rc -ne 0 ]] && echo "FAIL $asm rc=$rc: $(tail -1 "$out/stderr.log")"
  rm -f "$fna"
done < "$LIST"
rm -rf "$W"
echo done
