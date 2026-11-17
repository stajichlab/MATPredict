#!/usr/bin/bash -l
#SBATCH --job-name=orf-nr
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=96gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/orfnr_%j.log
# DIAMOND blastp of candidate ORF proteins against the local NCBI nr cluster database. Read-only discovery job.
# Env: QUERY (FASTA), OUT (folder), DB (diamond .dmnd; default the 2026-01-28 nr_cluster_seq.dmnd).
set -uo pipefail
: "${QUERY:?}" "${OUT:?}"
DB="${DB:-/bigdata/operations/pkgadmin/srv/projects/db/ncbi/diamond/20260128/nr_cluster_seq.dmnd}"
module load diamond/2.1.21
mkdir -p "$OUT"
diamond version > "$OUT/diamond_version.txt"
diamond blastp -d "$DB" -q "$QUERY" -o "$OUT/nr_hits.tsv" --more-sensitive --evalue 1e-3 --max-target-seqs 10 --threads ${SLURM_CPUS_PER_TASK:-16} \
  --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle --tmpdir "${SCRATCH:?}" > "$OUT/diamond.log" 2>&1
echo "done: $(wc -l < "$OUT/nr_hits.tsv") hits"
