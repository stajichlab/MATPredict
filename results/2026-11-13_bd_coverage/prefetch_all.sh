#!/usr/bin/bash
# Downloads every run in bd_runs.tsv with NCBI prefetch (login node). One copy only.
module load sratoolkit/3.2.0
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_bd_coverage/sra
mkdir -p "$D"
tail -n +2 /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_bd_coverage/bd_runs.tsv | cut -f2 | while read -r r; do
  [[ -s "$D/$r/$r.sra" ]] && continue
  prefetch "$r" -O "$D" --max-size 30G > "$D/prefetch_$r.log" 2>&1
done
echo finished > "$D/prefetch_done.txt"
