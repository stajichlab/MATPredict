#!/usr/bin/bash
# Candidate-only regression panel: reuses an existing base run, submits the candidate
# side (all panel groups plus Zygo 23) and a dependent diff. Same job shape as
# scripts/run_regression_panel.sh. The environment at submit time (for example
# MATPREDICT_MMSEQS_EVALUE, MMSEQS_BIN) passes to the jobs through --export=ALL.
#   CAND_WT=<frozen worktree> BASE_OUT=<existing panel dir with base/ and lists/> \
#   OUT=<new results dir> TAG=<short job-name tag> bash submit_cand_only.sh
set -euo pipefail
REPO="${REPO:-/bigdata/stajichlab/jstajich/projects/MATPredict}"
: "${CAND_WT:?}"; : "${BASE_OUT:?}"; : "${OUT:?}"; : "${TAG:?}"
PARTITION="${PARTITION:-short}"; TIME="${TIME:-2:00:00}"; MEM="${MEM:-32gb}"; CPUS="${CPUS:-16}"
[[ "$CAND_WT" == /bigdata/* && -f "$CAND_WT/src/MATPredict/__init__.py" ]] || { echo "bad CAND_WT" >&2; exit 1; }
[[ -d "$BASE_OUT/base" && -d "$BASE_OUT/lists" ]] || { echo "bad BASE_OUT" >&2; exit 1; }
mkdir -p "$OUT/lists" "$REPO/logs"
cp "$BASE_OUT"/lists/*.tsv "$OUT/lists/"; cp "$BASE_OUT/regression_panel.tsv" "$OUT/"
groups=$(ls "$OUT/lists" | sed 's/\.tsv$//')
jobids=()
for g in $groups; do
    j=$(sbatch --parsable --partition="$PARTITION" --time="$TIME" --mem="$MEM" --cpus-per-task="$CPUS" \
        --job-name="reg-$TAG-$g" --output="$REPO/logs/regression_%j.log" \
        --export=ALL,REPO="$REPO",SRC="$CAND_WT/src",MATPREDICT_DB_ROOT="$CAND_WT/db",CLADE="reg_${TAG}_$g",SAMPLE_LIST="$OUT/lists/$g.tsv",OUT_DIR="$OUT/cand/$g" \
        "$CAND_WT/scripts/run_clade_panel.slurm")
    jobids+=("$j"); echo "$j	cand	$g" | tee -a "$OUT/jobs.tsv"
done
j=$(sbatch --parsable --partition="$PARTITION" --time="$TIME" --mem="$MEM" --cpus-per-task="$CPUS" \
    --job-name="reg-$TAG-zygo23" --output="$REPO/logs/regression_%j.log" \
    --export=ALL,SRC="$CAND_WT/src",OUT="$OUT/cand/zygo23" "$CAND_WT/scripts/run_zygo_regression.slurm")
jobids+=("$j"); echo "$j	cand	zygo23" | tee -a "$OUT/jobs.tsv"
deps=$(IFS=:; echo "${jobids[*]}")
pairs=""; for g in $groups zygo23; do pairs+=" --pair $g $BASE_OUT/base/$g $OUT/cand/$g"; done
j=$(sbatch --parsable --partition="$PARTITION" --time=0:30:00 --mem=8gb --cpus-per-task=1 \
    --job-name="reg-$TAG-diff" --dependency="afterany:$deps" --output="$REPO/logs/regression_%j.log" \
    --wrap="set -e; export PYTHONPATH=$CAND_WT/src; P=$REPO/.pixi/envs/default/bin/python; \
\$P $CAND_WT/scripts/regression_check.py diff --title 'regression $TAG vs base' --out $OUT/diff $pairs; \
\$P $CAND_WT/scripts/check_record_selfcall.py --db $CAND_WT/db --reports $OUT/cand/*/runs --out $OUT/diff/record_selfcall_candidate.tsv || true")
echo "$j	diff	all" | tee -a "$OUT/jobs.tsv"
