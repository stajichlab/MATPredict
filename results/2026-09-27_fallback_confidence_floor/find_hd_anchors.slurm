#!/usr/bin/bash
#SBATCH --job-name=hd-anchor
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=24gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/hd_anchor_%j.log
# Locate MIP1 and beta-fg (Agaricus H97 proteins) in each genome with a
# phylum_fallback HD call: top tblastn hits, for the adjacency test.
set -uo pipefail
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_fallback_confidence_floor
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
W="${SCRATCH:?}/hd_anchor"; mkdir -p "$W" "$D/hd_hits"
one() {
  g="$1"; out="$D/hd_hits/$g.tsv"; [[ -s "$out" ]] && return 0
  zcat "$LIB/$g.fa.gz" > "$W/$g.fna" 2>/dev/null || { echo "FAIL $g unzip" >&2; return 0; }
  makeblastdb -in "$W/$g.fna" -dbtype nucl -out "$W/$g.db" >/dev/null 2>&1
  tblastn -query "$D/hd_anchor_queries.faa" -db "$W/$g.db" -evalue 1e-10 -seg no \
    -outfmt "6 qseqid sseqid pident length sstart send evalue bitscore" -max_target_seqs 5 \
    > "$W/$g.tmp" 2>/dev/null && mv "$W/$g.tmp" "$out"
  [[ -s "$out" ]] || echo -e "none\tnone\t0\t0\t0\t0\t1\t0" > "$out"
  rm -f "$W/$g.fna" "$W/$g.db".*
}
export -f one; export D LIB W
xargs -P 16 -I{} bash -c 'one {}' < "$D/hd_genomes.txt"
echo "done $(ls $D/hd_hits | wc -l) / $(wc -l < $D/hd_genomes.txt)"
