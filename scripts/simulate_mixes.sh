#!/bin/bash
# Simulated mixed samples for `matpredict reads-type`: N reads in total, FRAC percent from the MAT1-2 strain.
# usage: simulate_mixes.sh MAT11_STRAIN MAT12_STRAIN OUT_TSV   (env: PANEL_DIR, READS_DIR, PROJ_ROOT)
# No temporary FASTQ: each mixture is a process substitution. Absolute paths only (SLURM-safe).
set -eu
A=$1; B=$2; OUT=$3
N=${N_READS:-8000000}
: "${PANEL_DIR:?}" "${READS_DIR:?}" "${PROJ_ROOT:?}"
export PYTHONPATH=$PROJ_ROOT/.claude/worktrees/reads-type/src
PY=$PROJ_ROOT/.pixi/envs/test/bin/python
rm -f "$OUT.part"
for f in 0 1 2 5 10 20 50 80 90 95 98 99 100; do
  nb=$(( N * f / 100 )); na=$(( N - nb ))
  $PY -m MATPredict reads-type --idiomorph MAT1-1=$PANEL_DIR/MAT1-1.fasta --idiomorph MAT1-2=$PANEL_DIR/MAT1-2.fasta \
    --reads <(zcat $READS_DIR/${A}_R1_trimmed.fastq.gz | head -n $(( na * 4 ))) <(zcat $READS_DIR/${B}_R1_trimmed.fastq.gz | head -n $(( nb * 4 ))) \
    --sample ${A}+${B}_pct${f}_MAT1-2 --out "$OUT.$f.tsv" 2>/dev/null
  if [ ! -f "$OUT.part" ]; then head -1 "$OUT.$f.tsv" > "$OUT.part"; fi
  tail -n +2 "$OUT.$f.tsv" >> "$OUT.part"; rm -f "$OUT.$f.tsv"
done
mv "$OUT.part" "$OUT"
