#!/bin/bash
#SBATCH -p short -A stajichlab -c 1 --mem 6gb -t 2:00:00 -J cocci-reads --array=1-24 --output=/bigdata/stajichlab/jstajich/projects/MATPredict_cocci/pilot_out/reads/slurm_%A_%a.log
# matpredict reads-type of one pilot strain against each panel in $PANEL_LIST (default: the trimmed cocci panel; first 8M reads).
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/datasets/reads/coccidioides/pilot/env.sh
# awk, not `read`: tab is IFS whitespace, so an empty centre field would shift the columns
row=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$PILOT")
f() { echo "$row" | awk -F'\t' -v n=$1 '{print $n}'; }
strain=$(f 1); species=$(f 3); r1=$(f 6); r2=$(f 7); assembly=$(f 9)
IN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/input
printf 'sample\treads\n%s\t%s,%s\n' "$strain" "$IN/$r1" "$IN/$r2" > "$OUT/reads/samples_$strain.tsv"
for panel in ${PANEL_LIST:-cocci}; do
  python -m MATPredict reads-type --idiomorph MAT1-1=$PANELS/${panel}_MAT1-1.fasta --idiomorph MAT1-2=$PANELS/${panel}_MAT1-2.fasta \
    --samples "$OUT/reads/samples_$strain.tsv" --max-reads 8000000 --out "$OUT/reads/${strain}.${panel}_panel.tsv"
done
