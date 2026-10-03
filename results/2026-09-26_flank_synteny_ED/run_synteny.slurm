#!/usr/bin/bash
#SBATCH --job-name=flank-synteny-ED
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/flank_synteny_ED_%A.log
# Forward: all curated Mucoromycota MAT-locus proteins (tblastn) against each genome.
# Reverse: every flank-gene region's translated protein (blastp) against three
# Mucorales proteomes, for a reciprocal-best-hit ortholog call.
set -uo pipefail
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_flank_synteny_ED
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
W="${SCRATCH:?}/flank_synteny"; mkdir -p "$W" "$D/hits"
one() {
  g="$1"; out="$D/hits/$g.done"; [[ -e "$out" ]] && return 0
  zcat "$LIB/$g.fa.gz" > "$W/$g.fna" 2>/dev/null || { echo "FAIL $g unzip" >&2; return 0; }
  [[ -s "$W/$g.fna" ]] || { echo "FAIL $g empty" >&2; rm -f "$W/$g.fna"; return 0; }
  makeblastdb -in "$W/$g.fna" -dbtype nucl -out "$W/$g.db" >/dev/null 2>&1
  tblastn -query "$D/queries.faa" -db "$W/$g.db" -evalue 10 -seg no -num_threads 2 \
    -max_target_seqs 50 \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sseq" \
    > "$W/$g.tbn" 2>/dev/null || { echo "FAIL $g tblastn" >&2; rm -f "$W/$g".*; return 0; }
  python3 "$D/regions.py" "$W/$g.tbn" "$g" "$W/$g.regions.faa"
  if [[ -s "$W/$g.regions.faa" ]]; then
    blastp -query "$W/$g.regions.faa" -db "$D/ref/mucorales3" -evalue 1e-3 -num_threads 2 \
      -max_target_seqs 5 -outfmt "6 qseqid sseqid pident length evalue bitscore" \
      > "$W/$g.back" 2>/dev/null
  else : > "$W/$g.back"; fi
  cut -f1-10 "$W/$g.tbn" | zstd -q -c > "$D/hits/$g.tbn.zst"
  zstd -q -c "$W/$g.back" > "$D/hits/$g.back.zst"
  [[ -s "$W/$g.regions.faa" ]] && zstd -q -c "$W/$g.regions.faa" > "$D/hits/$g.regions.faa.zst"
  touch "$out"
  rm -f "$W/$g".*
}
export -f one; export D LIB W
cut -f1 "$D/genomes.tsv" | xargs -P 4 -I{} bash -c 'one {}'
echo "done $(ls $D/hits/*.done | wc -l) / $(wc -l < $D/genomes.tsv)"
