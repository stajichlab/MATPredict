#!/usr/bin/bash
#SBATCH --job-name=hd1-hmm
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/hd1hmm_%A.log
# Env: SRC, BASE, OUT. Builds the AMF HD1-like HMM and searches all genomes in BASE/list_all.txt.
set -uo pipefail
: "${SRC:?}" "${BASE:?}" "${OUT:?}"
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
module load hmmer/3.4
python3 "$SRC/scripts/hd1_hmm.py" build --base "$BASE" --out "$OUT" || exit 1
python3 "$SRC/scripts/hd1_hmm.py" search --base "$BASE" --out "$OUT" --lists "$BASE/list_all.txt"
echo done
