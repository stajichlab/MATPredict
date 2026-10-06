#!/usr/bin/bash
#SBATCH --job-name=sexMP-final-ml
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_sexMP_final_ml/trees_%j.log
set -euo pipefail
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_sexMP_final_ml
cd "$D"
module load iqtree/3.0.1 raxml-ng/2.0.2
iqtree3 -s selected_ids.afa -m MFP -B 1000 -alrt 1000 -T AUTO --threads-max 8 -pre iq -redo
BEST=$(grep -m1 'Best-fit model' iq.iqtree | sed 's/.*: //; s/ chosen.*//')
echo "IQ-TREE best model: $BEST"
M="$BEST"
if ! raxml-ng --parse --msa selected_ids.afa --model "$M" --prefix rx_check --redo >/dev/null 2>&1; then
  M=$(echo "$BEST" | sed -E 's/^[^+]+/LG/'); echo "RAxML-NG cannot use $BEST; trying $M"
  raxml-ng --parse --msa selected_ids.afa --model "$M" --prefix rx_check --redo >/dev/null 2>&1 || M="LG+G4+F"
fi
echo "RAxML-NG model: $M"
raxml-ng --all --msa selected_ids.afa --model "$M" --bs-trees autoMRE{500} --threads 8 --seed 2026 --prefix rx --redo
echo done
