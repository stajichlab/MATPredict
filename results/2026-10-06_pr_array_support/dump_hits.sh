#!/usr/bin/bash -l
#SBATCH -p short -c 4 --mem 16G -t 1:00:00 -J dumphits
R=/bigdata/stajichlab/jstajich/projects/MATPredict
W=$R/.claude/worktrees/run-pr-array-cand
O=$R/results/2026-10-06_pr_array_support/dump
mkdir -p $O
export PATH=$R/.pixi/envs/default/bin:$PATH PYTHONPATH=$W/src MATPREDICT_DB_ROOT=$W/db MATPREDICT_CACHE_DIR=$R/.matpredict_cache
for a in "GCF_000271585.1_Trametes_versicolor_v1.0 717944" "GCA_016772295.1_ASM1677229v1 1132390"; do set -- $a
 zcat /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/$1.fa.gz > $SCRATCH/$1.fna
 DUMP=$O/$1.pkl python $R/results/2026-10-06_pr_array_support/dump_hits.py --genome $SCRATCH/$1.fna --taxid $2 --out-dir $O/$1
done
