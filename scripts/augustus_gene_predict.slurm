#!/usr/bin/bash -l
#SBATCH --job-name=aug-pred
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=8:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/aug_%A.log
# Gene prediction with Augustus for assemblies that BFD skipped ("skipped_too_small"), so they can enter the flank-first scan.
# Env: SRC (worktree with this script), LIST (TSV: ASMID<TAB>name), OUT (folder), SPECIES (Augustus model, default rhizopus_oryzae).
# Output: $OUT/<name>/predict_results/<name>.proteins.fa.gz and <name>.gff3.gz (headers rewritten g1.t1 -> g1-T1 so
# flank_first_scan.py maps protein to gene). The model is not trained on these lineages: gene models are approximate.
set -uo pipefail
: "${SRC:?}" "${LIST:?}" "${OUT:?}"
SPECIES="${SPECIES:-rhizopus_oryzae}"
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
# `module load augustus` prints a conda error in a batch shell; set the paths directly (checked on the login node 2026-11-14)
AUG=/opt/linux/rocky/8.x/x86_64/pkgs/augustus/3.5.0/env
export PATH="$AUG/bin:$PATH"
export AUGUSTUS_CONFIG_PATH=/opt/linux/centos/8.x/x86_64/pkgs/augustus/3.5.0/env/config/
W="${SCRATCH:?}/aug"; mkdir -p "$W"
while IFS=$'\t' read -r asm name; do
  [[ -z "$asm" ]] && continue
  D="$OUT/$name/predict_results"; mkdir -p "$D"
  [[ -s "$D/$name.gff3.gz" ]] && continue
  zcat "$LIB/$asm.fa.gz" > "$W/$name.fa"
  # one Augustus run per contig chunk would be faster; the genomes here are 11-24 Mb so one run is enough
  augustus --species="$SPECIES" --gff3=on --protein=on --codingseq=off --softmasking=0 "$W/$name.fa" > "$W/$name.aug.gff3" 2> "$W/$name.aug.err" || { echo "FAIL $name"; tail -5 "$W/$name.aug.err"; continue; }
  python3 "$SRC/scripts/augustus_proteins.py" "$W/$name.aug.gff3" "$D/$name.proteins.fa.gz"
  awk -F'\t' '/^#/{next} $3=="gene"{print}' "$W/$name.aug.gff3" | gzip > "$D/$name.gff3.gz"
  echo "$name genes=$(zcat "$D/$name.gff3.gz" | wc -l) proteins=$(zcat "$D/$name.proteins.fa.gz" | grep -c '>')"
  rm -f "$W/$name".*
done < "$LIST"
echo done
