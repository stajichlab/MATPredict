#!/bin/bash
#SBATCH -p short -A stajichlab -c 1 --mem 4gb -t 1:00:00 -J cocci-profile --output=/bigdata/stajichlab/jstajich/projects/MATPredict_cocci/pilot_out/profile_%j.log
# Per-window depth of the trimmed MAT1-1 panel in a MAT1-2 strain (B3286) and a MAT1-1 strain (B11299), 2M reads each.
source /bigdata/stajichlab/jstajich/projects/MATPredict_cocci/datasets/reads/coccidioides/pilot/env.sh
IN=/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/input
for s in B3286 B11299; do
  row=$(awk -F'\t' -v s=$s '$1==s' "$PILOT"); r1=$(echo "$row" | cut -f6); r2=$(echo "$row" | cut -f7)
  echo "## $s MAT1-1 panel"; python $CLONE/datasets/reads/coccidioides/pilot/profile_panel.py $PANELS/cocci_MAT1-1.fasta $IN/$r1 $IN/$r2 2000000 500
done
