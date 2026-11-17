#!/usr/bin/bash -l
#SBATCH --job-name=bd-multiref
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/bd_multiref_%j.log
# Align four more Bd assemblies to CLFT044 (same preset as job 32089404) and call lineage-diagnostic divergent sites.
set -euo pipefail
module load minimap2/2.30
SRC=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/bd-coverage
OUT=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-15_bd_multiref
REF=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-13_bd_coverage/ref
G=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
W="${SCRATCH:?}"
cd "$W"
declare -A FA=( [CLFT067]=GCA_036289345.1_UCR_BdCLFT067_1.0.fa.gz [CLFT071]=GCA_029704095.1_UCR_BdCLFT071_1.0.fa.gz [JAM81]=GCF_000203795.2_v1.1.fa.gz [RTP6]=GCA_003595275.1_ASM359527v1.fa.gz )
for n in "${!FA[@]}"; do
  zcat "$G/${FA[$n]}" > $n.fa
  minimap2 -x asm10 -c --cs=short -t ${SLURM_CPUS_PER_TASK} $REF/CLFT044.fa $n.fa 2> $OUT/$n.minimap2.log | gzip > $OUT/${n}_vs_CLFT044.paf.gz
done
python3 "$SRC/scripts/bd_ref_lineage_sites.py" --sites $REF/divergent_sites.tsv \
  --paf JEL423=$REF/JEL423_vs_CLFT044.paf.gz --paf CLFT067=$OUT/CLFT067_vs_CLFT044.paf.gz --paf CLFT071=$OUT/CLFT071_vs_CLFT044.paf.gz \
  --paf JAM81=$OUT/JAM81_vs_CLFT044.paf.gz --paf RTP6=$OUT/RTP6_vs_CLFT044.paf.gz \
  --brazil CLFT067 CLFT071 --other JEL423 JAM81 RTP6 --out $OUT/sites
echo done
