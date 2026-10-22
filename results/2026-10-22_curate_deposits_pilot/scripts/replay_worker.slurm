#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J cp-replay -p short
# One array task = one chunk (chunks/chunk_NNN: ASM<TAB>TAXID). `detect` with production defaults plus the clean-run
# settings (--partial-core-policy coverage_profile --confidence-cap-policy off --strong-single-gene-policy off; no
# --homothallic-policy). REPO selects the code+db (base = origin/main, pilot = curate-deposits-pilot). Resumable.
set -uo pipefail
REPO=${CP_REPO:?}
SUB=${CP_SUB:?}
P=/bigdata/stajichlab/jstajich/projects/MATPredict/curate_pilot_replay
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=$ENV/bin:$PATH PYTHONPATH=$REPO/src MATPREDICT_DB_ROOT=$REPO/db PYTHONNOUSERSITE=1
export MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache MATPREDICT_HTML=0
export MATPREDICT_CONFIG=$P/empty_config.toml   # no e-mail is sent to NCBI
N=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
LIST=$P/chunks/chunk_$N
OUT=$P/runs/$SUB
W=${SCRATCH:-/tmp}/cp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$W" "$OUT"
where=$(python -c 'import MATPredict;print(MATPredict.__file__)')
case "$where" in "$REPO"/src/*) ;; *) echo "WRONG CODE: $where" >&2; exit 1;; esac
while IFS=$'\t' read -r asm taxid; do
  out=$OUT/$asm
  [[ -s $out/detection_report.yaml ]] && continue
  mkdir -p "$out"; fna=$W/$asm.fna
  zcat "$LIB/$asm.fa.gz" > "$fna" || { echo "FAIL $asm decompress"; continue; }
  /usr/bin/time -f "%e s %M KB" -o "$out/time.txt" timeout 3600 python -m MATPredict detect --genome "$fna" --taxid "$taxid" \
     --out-dir "$out" --partial-core-policy coverage_profile --confidence-cap-policy off \
     --strong-single-gene-policy off > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; [[ $rc -ne 0 ]] && echo "FAIL $asm rc=$rc: $(tail -1 "$out/stderr.log")"
  rm -f "$fna"
done < "$LIST"
rm -rf "$W"; echo done
