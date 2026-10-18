#!/bin/bash
#SBATCH -p short -A stajichlab -c 3 --mem 12gb -t 2:00:00 -J cocci-full-detect --array=1-558%20 --output=/bigdata/stajichlab/jstajich/projects/MATPredict_cocci/full_out/detect/slurm_%A_%a.log
# matpredict detect on the freeze assembly of one strain (row N of samples.tsv, header excluded).
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/datasets/reads/coccidioides/pilot/env.sh
export OUT=/bigdata/stajichlab/jstajich/projects/MATPredict_cocci/full_out PILOT=/bigdata/stajichlab/jstajich/projects/MATPredict_cocci/datasets/reads/coccidioides/full/samples.tsv
# awk, not `read`: tab is IFS whitespace, so an empty centre field would shift the columns
row=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$PILOT")
f() { echo "$row" | awk -F'\t' -v n=$1 '{print $n}'; }
strain=$(f 1); species=$(f 3); r1=$(f 6); r2=$(f 7); assembly=$(f 9)
[[ -n "$assembly" ]] || { echo "skip $strain (no assembly)"; exit 0; }
case "$species" in *immitis) taxid=5501;; *posadasii) taxid=199306;; esac
python -m MATPredict detect --genome "$assembly" --taxid $taxid --out-dir "$OUT/detect/runs/$strain" > "$OUT/detect/runs/$strain.log" 2>&1 && echo "ok $strain" || echo "FAIL $strain"
