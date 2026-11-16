#!/usr/bin/bash
#SBATCH --job-name=hd-allele
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/hdallele_%A.log
set -uo pipefail
: "${SRC:?}" "${BASE:?}" "${OUT:?}"
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
python3 "$SRC/scripts/hd_allele_identity.py" --base "$BASE" --out "$OUT"
