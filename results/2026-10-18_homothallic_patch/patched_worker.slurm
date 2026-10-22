#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J hl-run -p short
# One array task = one chunk (chunks/chunk_NNN.tsv: ASMID<TAB>TAXID). `detect` of the homothallic-loci-patch worktree with
# --homothallic-policy both_genes and --partial-core-policy coverage_profile (the audit's setting); every other default.
# Resumable. The report's `homothallic_note.previous_locus_class` carries the legacy class of every re-classed locus.
set -uo pipefail
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/hl-patch
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
OUT=$REPO/results/2026-10-18_homothallic_patch/runs_out
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=$ENV/bin:$PATH PYTHONPATH=$REPO/src MATPREDICT_DB_ROOT=$REPO/db PYTHONNOUSERSITE=1
export MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache MATPREDICT_HTML=0
FLAGS=${HL_FLAGS:---homothallic-policy both_genes}
SUB=${HL_SUB:-runs}
N=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
LIST=$REPO/results/2026-10-18_homothallic_patch/${HL_CHUNKS:-chunks}/chunk_$N.tsv
W=${SCRATCH:?}/hl_${SLURM_JOB_ID}_$N
mkdir -p "$W" "$OUT/$SUB"
where=$(python -c 'import MATPredict;print(MATPredict.__file__)')
case "$where" in "$REPO"/src/*) ;; *) echo "WRONG CODE: $where" >&2; exit 1;; esac
while IFS=$'\t' read -r asm taxid; do
  out=$OUT/$SUB/$asm
  [[ -s $out/detection_report.yaml ]] && continue
  mkdir -p "$out"; fna=$W/$asm.fna
  zcat "$LIB/$asm.fa.gz" > "$fna" || { echo "FAIL $asm decompress"; continue; }
  timeout 3600 python -m MATPredict detect --genome "$fna" --taxid "$taxid" --out-dir "$out" \
     --partial-core-policy coverage_profile $FLAGS > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; [[ $rc -ne 0 ]] && echo "FAIL $asm rc=$rc: $(tail -1 "$out/stderr.log")"
  rm -f "$fna"
done < "$LIST"
rm -rf "$W"; echo done
