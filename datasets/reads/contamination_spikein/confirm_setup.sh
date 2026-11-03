#!/bin/bash -l
# Build the bowtie2 index of the two Coccidioides reference genomes (C. immitis RS and C. posadasii C735) used by confirm_flagged.slurm.
#SBATCH -p short -c 8 --mem 16gb -t 1:00:00 -J confirm_setup -o %x.%j.log
set -euo pipefail
module load bowtie2/2.5.3 2>/dev/null || true
G=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
OUT=${OUT:?output folder}; mkdir -p "$OUT/ref"
zcat "$G/GCF_000149335.2_ASM14933v2.fa.gz" "$G/GCF_000151335.2_JCVI-cpa1-1.0.fa.gz" > "$OUT/ref/cocci.fa"
bowtie2-build --threads "${SLURM_CPUS_PER_TASK:-8}" "$OUT/ref/cocci.fa" "$OUT/ref/cocci" > "$OUT/ref/build.log" 2>&1
ls -la "$OUT/ref" | head
