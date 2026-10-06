#!/bin/bash
#SBATCH -p exfab -A exfab -c 4 --mem 16gb -t 2:00:00 -J afum-hyb --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-05_afum_hybrids_detect/slurm_%j.log
# MATPredict detect on the 3 A. fumigatus putative hybrids (Lofgren et al. 2022), main @ 2569348.
set -e
W=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-2569348
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH PYTHONPATH=$W/src MATPREDICT_DB_ROOT=$W/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
for s in IFM_59359 DMC_AF100-1_3 IFM_61407; do
  python -m MATPredict detect --genome /bigdata/stajichlab/shared/projects/Afumigatus_pangenome/genomes/$s.scaffolds.fasta --taxid 746128 --out-dir /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-05_afum_hybrids_detect/$s > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-05_afum_hybrids_detect/$s.log 2>&1 || echo "FAIL $s"
done
echo done
