#!/usr/bin/bash
# HPCC: v0.8.1 on the 134 Ogataea genomes with and without the six candidate MATog records (<= 40 CPUs in all).
set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict
W=$R/.claude/worktrees/run-fix228
O=$R/results/2026-11-04_matog_candidates/loo
S=$W/src; SLURM=$W/scripts/run_clade_panel.slurm
sub() {  # tag db list cpus extra-args
  DETECT_ARGS="--preset v0.8.1 $5" sbatch --parsable -p short --time=2:00:00 --mem=24gb --cpus-per-task=$4 -J og-$1 -o $R/logs/og_$1_%j.log \
   --export=ALL,REPO=$R,SRC=$S,MATPREDICT_DB_ROOT=$2,CLADE=og-$1,SAMPLE_LIST=$3,OUT_DIR=$O/$1,GENOME_TIMEOUT=3600 $SLURM
}
SRC_ACC="GCA_030570775.1 GCA_030559035.1 GCA_030568015.1 GCA_030574315.1 GCA_002008015.1 GCA_002008025.1"
RID=(51923_nrrl-y-17664_MATog_MATa-MATalpha 1249076_nrrl-yb-4943_MATog_MATa-MATalpha 169506_nrrl-y-27599_MATog_MATa-MATalpha 246385_nrrl-y-27963_MATog_MATa-MATalpha 5477_ndk27a1_MATog_MATalpha 5477_uniss-cb60_MATog_MATa)
grep -v -E "$(echo $SRC_ACC | tr ' ' '|')" $O/lists/all.tsv > $O/lists/others.tsv
wc -l $O/lists/others.tsv
sub base $R/.claude/worktrees/run-a2933f7/db $O/lists/all.tsv 8 ""
sub cand $O/db_cand $O/lists/others.tsv 8 ""
i=0
for a in $SRC_ACC; do
  grep "$a" $O/lists/all.tsv > $O/lists/self_$i.tsv
  sub self$i $O/db_cand $O/lists/self_$i.tsv 2 "--exclude-records ${RID[$i]}"
  i=$((i+1))
done
