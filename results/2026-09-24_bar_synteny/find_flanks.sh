#!/usr/bin/bash
#SBATCH --job-name=bar-synteny
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=48gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/bar_synteny_%A.log
# Independent positional evidence for the bar question: locate each genome's
# SLA2 and APN2 orthologs (top tblastn hit) so a MAT candidate can be checked
# for adjacency. Pezizomycotina MAT sits between SLA2 and APN2.
set -uo pipefail
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-24_bar_synteny
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
W="${SCRATCH:?}/bar_synteny"; mkdir -p "$W" "$D/hits"
cat /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-24_pilots/lists/*.tsv | cut -f1 | sort -u > "$W/genomes.txt"
one() {
  g="$1"; out="$D/hits/$g.tsv"; [[ -s "$out" ]] && return 0
  zcat "$LIB/$g.fa.gz" > "$W/$g.fna" 2>/dev/null || { echo "FAIL $g unzip" >&2; return 0; }
  [[ -s "$W/$g.fna" ]] || { echo "FAIL $g empty" >&2; rm -f "$W/$g.fna"; return 0; }
  makeblastdb -in "$W/$g.fna" -dbtype nucl -out "$W/$g.db" >/dev/null 2>&1
  tblastn -query "$D/flank_queries.faa" -db "$W/$g.db" -evalue 1e-10 -seg no \
    -outfmt "6 qseqid sseqid pident length sstart send evalue bitscore" -max_target_seqs 5 \
    > "$W/$g.tmp" 2>/dev/null && mv "$W/$g.tmp" "$out"
  [[ -s "$out" ]] || echo -e "none\tnone\t0\t0\t0\t0\t1\t0" > "$out"
  rm -f "$W/$g.fna" "$W/$g.db".*
}
export -f one; export D LIB W
xargs -P 16 -I{} bash -c 'one {}' < "$W/genomes.txt"
echo "done $(ls $D/hits | wc -l) / $(wc -l < $W/genomes.txt)"
