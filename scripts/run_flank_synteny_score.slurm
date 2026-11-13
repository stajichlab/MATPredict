#!/usr/bin/bash
#SBATCH --job-name=ffs-synt
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/ffs2_%A.log
# Flank-first scan, stage 2. Env: SRC, STAGE1 (per_genome dir), LIST, OUT; optional GROUP_BY (genus|all), REF.
set -uo pipefail
: "${SRC:?}" "${STAGE1:?}" "${LIST:?}" "${OUT:?}"
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
module load hmmer/3.4
ARGS=(--stage1 "$STAGE1" --list "$LIST" --out "$OUT" --group-by "${GROUP_BY:-genus}" --threads 8)
[[ -n "${REF:-}" ]] && ARGS+=(--ref "$REF")
python3 "$SRC/scripts/flank_synteny_score.py" "${ARGS[@]}"
echo done
