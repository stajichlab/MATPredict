#!/usr/bin/bash
#SBATCH --job-name=sexMP-tree
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/logs/sexMP_tree_%j.log
# IQ-TREE on the HMM-guided HMG-box alignment: ModelFinder, 1000 UFBoot, 1000 SH-aLRT.
set -euo pipefail
module load iqtree/3.0.1
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_sexMP_phylogeny
cd "$D"
iqtree3 -s aln_hmm.ids.fa -st AA -m MFP -B 1000 -alrt 1000 -T ${SLURM_CPUS_PER_TASK:-4} --prefix tree_hmm -redo
