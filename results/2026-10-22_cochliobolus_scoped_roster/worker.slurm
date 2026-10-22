#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J cochl-scoped -p short
# One array task = one chunk ($WD/chunks_$CHUNKS/chunk_NNN.tsv: ASMID<TAB>TAXID).
# Env: VARIANT (output subdir), CODE (worktree), FLAGS (extra detect flags), CHUNKS (chunk set name).
# Clean-run settings on every variant: --partial-core-policy coverage_profile --confidence-cap-policy off
# --strong-single-gene-policy off (unless FLAGS overrides the last), no --homothallic-policy. Resumable.
set -uo pipefail
WD=/bigdata/stajichlab/jstajich/projects/cochl_scoped
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
CODE=${CODE:?}
export PATH=$ENV/bin:$PATH PYTHONPATH=$CODE/src MATPREDICT_DB_ROOT=$CODE/db PYTHONNOUSERSITE=1
export MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache MATPREDICT_HTML=0
N=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
LIST=$WD/chunks_${CHUNKS:?}/chunk_$N.tsv
OUT=$WD/runs/${VARIANT:?}
W=${SCRATCH:?}/cochl_${SLURM_JOB_ID}_$N
mkdir -p "$W" "$OUT"
where=$(python -c 'import MATPredict;print(MATPredict.__file__)')
case "$where" in "$CODE"/src/*) ;; *) echo "WRONG CODE: $where" >&2; exit 1;; esac
echo "code=$where host=$(hostname) chunk=$N variant=$VARIANT flags=${FLAGS:-}"
while IFS=$'\t' read -r asm taxid; do
  out=$OUT/$asm
  [[ -s $out/detection_report.yaml ]] && continue
  mkdir -p "$out"; fna=$W/$asm.fna
  zcat "$LIB/$asm.fa.gz" > "$fna" || { echo "FAIL $asm decompress"; continue; }
  t0=$(date +%s)
  timeout 3600 python -m MATPredict detect --genome "$fna" --taxid "$taxid" --out-dir "$out" \
     --partial-core-policy coverage_profile --confidence-cap-policy off --strong-single-gene-policy off \
     ${FLAGS:-} > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; echo $(( $(date +%s) - t0 )) > "$out/wall_seconds"
  [[ $rc -ne 0 ]] && echo "FAIL $asm rc=$rc: $(tail -1 "$out/stderr.log")"
  rm -f "$fna" "$out/_reference.faa"
done < "$LIST"
rm -rf "$W"; echo done
