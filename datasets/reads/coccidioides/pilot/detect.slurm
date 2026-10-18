#!/bin/bash
#SBATCH -p short -A stajichlab -c 3 --mem 12gb -t 2:00:00 -J cocci-detect --array=1-24 --output=/bigdata/stajichlab/jstajich/projects/MATPredict_cocci/pilot_out/detect/slurm_%A_%a.log
# matpredict detect on the freeze assembly of one pilot strain (row N of pilot_samples.tsv, header excluded).
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/datasets/reads/coccidioides/pilot/env.sh
# awk, not `read`: tab is IFS whitespace, so an empty centre field would shift the columns
row=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$PILOT")
f() { echo "$row" | awk -F'\t' -v n=$1 '{print $n}'; }
strain=$(f 1); species=$(f 3); r1=$(f 6); r2=$(f 7); assembly=$(f 9)
case "$species" in *immitis) taxid=5501;; *posadasii) taxid=199306;; esac
python -m MATPredict detect --genome "$assembly" --taxid $taxid --out-dir "$OUT/detect/runs/$strain" > "$OUT/detect/runs/$strain.log" 2>&1 && echo "ok $strain" || echo "FAIL $strain"
