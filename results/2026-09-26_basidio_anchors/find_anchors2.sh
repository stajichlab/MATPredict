#!/usr/bin/bash
#SBATCH --job-name=basidio-anchors2
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=24gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/basidio_anchors2_%A.log
# Positional test for Basidiomycota MAT anchors: tblastn MIP1, beta-fg, ICMT,
# and every curated HD-class and receptor-class protein against each pilot
# genome, so HD/receptor hits can be checked for adjacency to the anchors.
set -uo pipefail
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_basidio_anchors
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
W="${SCRATCH:?}/basidio_anchors"; mkdir -p "$W" "$D/hits2"
one() {
  g="$1"; out="$D/hits2/$g.tsv"; [[ -s "$out" ]] && return 0
  zcat "$LIB/$g.fa.gz" > "$W/$g.fna" 2>/dev/null || { echo "FAIL $g unzip" >&2; return 0; }
  makeblastdb -in "$W/$g.fna" -dbtype nucl -out "$W/$g.db" >/dev/null 2>&1
  tblastn -query "$D/queries2.faa" -db "$W/$g.db" -evalue 1e-5 -seg no -num_threads 2 \
    -outfmt "6 qseqid sseqid pident length sstart send evalue bitscore" -max_target_seqs 20 \
    > "$W/$g.tmp" 2>/dev/null && mv "$W/$g.tmp" "$out"
  [[ -s "$out" ]] || echo -e "none\tnone\t0\t0\t0\t0\t1\t0" > "$out"
  awk '/^>/{n=substr($1,2)} !/^>/{l[n]+=length($0)} END{for(k in l) print k"\t"l[k]}' "$W/$g.fna" > "$D/hits2/$g.len"
  rm -f "$W/$g.fna" "$W/$g.db".*
}
export -f one; export D LIB W
cut -f1 "$D/pilot.tsv" | xargs -P 8 -I{} bash -c 'one {}'
echo "done $(ls $D/hits2/*.tsv | wc -l) / $(wc -l < $D/pilot.tsv)"
