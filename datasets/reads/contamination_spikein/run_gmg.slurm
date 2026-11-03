#!/usr/bin/env bash
# Screen the mixed read sets with GmG at each depth. --max-reads counts reads PER FILE, so N is N pairs.
# Usage: sbatch --export=ALL,OUT=/path/spikein run_gmg.slurm      (DEPTHS="100000 500000 2000000 6000000")
#SBATCH -p short -c 16 --mem 48gb -t 2:00:00 -J spikein_gmg -o %x.%j.log
set -euo pipefail
: "${OUT:?output folder with samples.tsv}"
GMG=${GMG:-/bigdata/stajichlab/shared/projects/GmG_runs/bin/gmg}
IDX=${IDX:-/bigdata/stajichlab/shared/projects/GmG_runs/release/2026.10.5/bfd-2026.10.5-s10000.gmgi}
for N in ${DEPTHS:-100000 500000 2000000 6000000}; do
  "$GMG" run --sample-sheet "$OUT/samples.tsv" --index "$IDX" -o "$OUT/gmg/depth_$N" \
    --max-reads "$N" --threads "${SLURM_CPUS_PER_TASK:-16}" --resume
done
echo "gmg spike-in screen complete"
