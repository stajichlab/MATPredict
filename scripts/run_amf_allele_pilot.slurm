#!/usr/bin/bash
#SBATCH --job-name=amf-allele
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=24gb
#SBATCH --time=6:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/amfallele_%A.log
# Allele read counts at the AMF HD locus from long reads. Env: SRC, BASE, SRA (prefetch dir), OUT.
# RUNS: space-separated SRR:platform:label (platform pb|ont).
set -uo pipefail
: "${SRC:?}" "${BASE:?}" "${SRA:?}" "${OUT:?}" "${RUNS:?}"
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
module load sratoolkit/3.2.0
module load minimap2/2.30
module load samtools/1.22.1
W="${SCRATCH:?}/amfallele"; mkdir -p "$W" "$OUT"
python3 "$SRC/scripts/amf_allele_refs.py" "$BASE" "$OUT/locus_refs.fa" \
  Rhizophagus_irregularis_A4:A4 Rhizophagus_irregularis_A5:A5 Rhizophagus_irregularis_C2:C2 \
  Rhizophagus_irregularis_A1:A1 Rhizophagus_irregularis_DAOM-197198:DAOM197198 \
  Rhizophagus_irregularis_DAOM_240409:DAOM240409 Rhizophagus_irregularis_B3:B3 | tee "$OUT/refs.txt"
for spec in $RUNS; do
  IFS=: read -r run plat label <<< "$spec"
  fasterq-dump "$SRA/$run/$run.sra" -e 8 -O "$W" --temp "$W" -f > "$OUT/fasterq_$run.log" 2>&1 || { echo "FAIL fasterq $run"; continue; }
  fq=$(ls "$W"/${run}*.fastq | head -1)
  nreads=$(( $(wc -l < "$fq") / 4 ))
  preset=map-pb; [[ "$plat" == ont ]] && preset=map-ont
  minimap2 -ax $preset -t 8 --secondary=no "$OUT/locus_refs.fa" "$fq" 2> "$OUT/minimap_$run.log" | samtools sort -@ 4 -o "$W/$run.bam" - 
  samtools index "$W/$run.bam"
  python3 "$SRC/scripts/amf_allele_read_count.py" "$W/$run.bam" 0.8 > "$OUT/counts_${label}_${run}.tsv"
  python3 "$SRC/scripts/amf_allele_read_count.py" "$W/$run.bam" 0.5 > "$OUT/counts50_${label}_${run}.tsv"
  samtools view -F 0x904 "$W/$run.bam" | awk '{print $1"\t"$3"\t"$4"\t"$6"\t"$12"\t"$13}' | cut -c1-200 > "$OUT/primary_${label}_${run}.tsv"
  echo "$run $label reads=$nreads mapped_primary=$(samtools view -c -F 0x904 "$W/$run.bam")" | tee -a "$OUT/summary.txt"
  [[ "${KEEP:-0}" == 1 ]] && cp "$W/$run.bam" "$W/$run.bam.bai" "$OUT/" ; rm -f "$W"/${run}*.fastq "$W/$run.bam"*
done
echo done
