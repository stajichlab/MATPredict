#!/bin/bash
#SBATCH -p exfab -A exfab -c 4 --mem 8gb -t 2:00:00 -J offtax --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-04_offline_taxonomy_check/slurm_%j.log
# Offline taxonomy check: detect with --taxid, MATPREDICT_OFFLINE=1, empty NCBI cache, no network use.
set -e
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/offline-taxonomy/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/offline-taxonomy/db
export MATPREDICT_TAXONOMY=/bigdata/stajichlab/jstajich/projects/MATPredict/results/ncbi_taxonomy_2026-10-01.tsv.zst MATPREDICT_OFFLINE=1
export MATPREDICT_CACHE_DIR=${SCRATCH:?}/emptycache; mkdir -p $MATPREDICT_CACHE_DIR
zcat /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/GCA_001661245.1_Pacta1_2.fa.gz > $SCRATCH/g.fna
/usr/bin/time -f '%e s %M KB' python -m MATPredict detect --genome $SCRATCH/g.fna --taxid 669874 --out-dir /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-04_offline_taxonomy_check/pacta > /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-04_offline_taxonomy_check/stdout.log 2> /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-04_offline_taxonomy_check/stderr.log
echo "cache files after run: $(ls $MATPREDICT_CACHE_DIR | wc -l)"
