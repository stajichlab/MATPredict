#!/usr/bin/bash
#SBATCH --job-name=ffs-scan
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/ffs_%A.log
# Flank-first scan, stage 1. Env: SRC (worktree), OUT (results folder), LIST (genome dir names, one per line),
# CONTROLS (control FASTA), WINDOW (genes per side, default 5).
set -uo pipefail
: "${SRC:?}" "${OUT:?}" "${LIST:?}" "${CONTROLS:?}"
WINDOW="${WINDOW:-5}"
ANN=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation
PIXI=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
module load hmmer/3.4
export PATH="$PIXI:$PATH"
W="${SCRATCH:?}/ffs"; mkdir -p "$W" "$OUT/per_genome"
diamond makedb --in "$CONTROLS" -d "$W/controls" --quiet
HMM="$OUT/pfam/hmg_homeo.hmm"
while read -r g; do
  [[ -z "$g" ]] && continue
  P="$ANN/$g/predict_results"
  [[ -s "$OUT/per_genome/$g.loci.tsv" ]] && continue
  if [[ ! -s "$P/$g.proteins.fa.gz" || ! -s "$P/$g.gff3.gz" ]]; then echo "MISSING $g" >&2; continue; fi
  python3 "$SRC/scripts/flank_first_scan.py" --name "$g" --proteins "$P/$g.proteins.fa.gz" --gff3 "$P/$g.gff3.gz" \
    --pfam-hmm "$HMM" --controls-dmnd "$W/controls.dmnd" --controls-faa "$CONTROLS" --out "$OUT/per_genome" \
    --window "$WINDOW" --threads 4 || echo "FAIL $g" >&2
done < "$LIST"
echo done
