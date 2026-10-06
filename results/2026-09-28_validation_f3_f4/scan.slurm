#!/bin/bash
#SBATCH -p short -c 16 --mem 24G -t 2:00:00
set -u
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
export PATH="$E:$PATH"
module load miniprot 2>/dev/null || true
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-28_validation_f3_f4
SCAN=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_pheromone_positional/scan_genome.py
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_pheromone_positional
while read a; do [ -z "$a" ] && continue; [ -s "$D/out_neg/$a.loci.tsv" ] && continue
  "$E/python" "$SCAN" "$a" "$D/out_neg" || echo "FAIL $a" >&2
done < "$1"
