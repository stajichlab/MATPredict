set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict; O=$R/results/2026-11-04_148_hyphopichia
cd $R && git fetch -q origin floor-1b-experiment
W=$R/.claude/worktrees/run-floor1b; [ -d $W ] || git worktree add -q $W origin/floor-1b-experiment
cd $W && git checkout -q --detach origin/floor-1b-experiment; rm -rf $W/db; ln -s $R/.claude/worktrees/run-a2933f7/db $W/db
for tag in base1b:$R/.claude/worktrees/run-a2933f7/db cand1b:$O/testdb; do n=${tag%%:*}; d=${tag#*:}; rm -rf $O/$n
DETECT_ARGS="--preset v0.8.1" sbatch --parsable -p short --time=2:00:00 --mem=24gb --cpus-per-task=8 -J h148-$n -o $R/logs/h148_${n}_%j.log \
 --export=ALL,MATPREDICT_FLOOR_EXPERIMENT=1b,REPO=$R,SRC=$W/src,MATPREDICT_DB_ROOT=$d,CLADE=h148-$n,SAMPLE_LIST=$O/lists/s400.tsv,OUT_DIR=$O/$n,GENOME_TIMEOUT=3600 $W/scripts/run_clade_panel.slurm; done
