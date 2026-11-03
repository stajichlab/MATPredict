#!/usr/bin/env bash
# One array task = one condition of the spike-in: build its mixed read files (see spikein.py mix-one).
# Usage: sbatch --array=1-38%12 --export=ALL,SRC=sources.tsv,OUT=/path/spikein run_mix.slurm
#SBATCH -p short -c 2 --mem 8gb -t 2:00:00 -J spikein_mix -o %x.%A_%a.log
set -euo pipefail
: "${SRC:?sources.tsv}" "${OUT:?output folder with conditions.tsv}"
HERE=${HERE:-$(cd "$(dirname "$0")" && pwd)}
python3 "$HERE/spikein.py" mix-one --sources "$SRC" --outdir "$OUT" --index "$SLURM_ARRAY_TASK_ID" --pairs "${PAIRS:-6000000}"
