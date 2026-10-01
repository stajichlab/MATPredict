#!/usr/bin/bash
# Pre-sign-off regression check (curator ruling, J. Stajich, 2026-09-30).
#
# Runs the fixed regression panel (testset/regression_panel.tsv) plus Zygo 23 on
# two FROZEN worktrees -- baseline and candidate -- then diffs them with
# scripts/regression_check.py. Every new record or classifier rebuild must
# attach the resulting regression_summary.md before sign-off.
# See docs/regression-check.md.
#
#   BASE_WT=/bigdata/.../.claude/worktrees/run-<base-sha> \
#   CAND_WT=/bigdata/.../.claude/worktrees/run-<cand-sha> \
#   OUT=/bigdata/.../results/<date>_regression_<name> \
#   bash scripts/run_regression_panel.sh
#
# Both worktrees must be on /bigdata (/scratch is node-local) and must not be
# edited while the jobs run. Paths come from the environment, never from
# BASH_SOURCE (it does not resolve inside SLURM work directories).
# Job sizing: each panel group is one `short` job (<= 2 h). Measured medians:
# Mucoromycota ~130 s/genome, Ascomycota MAT ~620 s, Basidiomycota ~150 s, at
# 16-way parallelism -> ~0.2-0.7 h per group.
set -euo pipefail

REPO="${REPO:-/bigdata/stajichlab/jstajich/projects/MATPredict}"
: "${BASE_WT:?set BASE_WT to the baseline frozen worktree}"
: "${CAND_WT:?set CAND_WT to the candidate frozen worktree}"
: "${OUT:?set OUT to the results directory}"
PANEL="${PANEL:-$CAND_WT/testset/regression_panel.tsv}"
PARTITION="${PARTITION:-short}"
TIME="${TIME:-2:00:00}"
MEM="${MEM:-32gb}"
CPUS="${CPUS:-16}"

for wt in "$BASE_WT" "$CAND_WT"; do
    [[ -f "$wt/src/MATPredict/__init__.py" && -f "$wt/db/Ascomycota/order.yml" ]] \
        || { echo "not a MATPredict worktree with src and db: $wt" >&2; exit 1; }
    [[ "$wt" == /bigdata/* ]] || { echo "worktree must be on /bigdata: $wt" >&2; exit 1; }
done
[[ -s "$PANEL" ]] || { echo "panel missing: $PANEL" >&2; exit 1; }

mkdir -p "$OUT/lists" "$REPO/logs"
cp "$PANEL" "$OUT/regression_panel.tsv"
groups=$(awk -F'\t' '!/^#/ && $1!="asmid" {print $3}' "$PANEL" | sort -u)
for g in $groups; do
    awk -F'\t' -v g="$g" '!/^#/ && $3==g {print $1"\t"$2}' "$PANEL" > "$OUT/lists/$g.tsv"
done

jobids=()
for side in base cand; do
    wt="$BASE_WT"; [[ "$side" == cand ]] && wt="$CAND_WT"
    for g in $groups; do
        j=$(sbatch --parsable --partition="$PARTITION" --time="$TIME" --mem="$MEM" \
            --cpus-per-task="$CPUS" --job-name="reg-$side-$g" \
            --output="$REPO/logs/regression_%j.log" \
            --export=ALL,REPO="$REPO",SRC="$wt/src",MATPREDICT_DB_ROOT="$wt/db",CLADE="reg_${side}_$g",SAMPLE_LIST="$OUT/lists/$g.tsv",OUT_DIR="$OUT/$side/$g" \
            "$wt/scripts/run_clade_panel.slurm")
        jobids+=("$j"); echo "$j	$side	$g" | tee -a "$OUT/jobs.tsv"
    done
    j=$(sbatch --parsable --partition="$PARTITION" --time="$TIME" --mem="$MEM" \
        --cpus-per-task="$CPUS" --job-name="reg-$side-zygo23" \
        --output="$REPO/logs/regression_%j.log" \
        --export=ALL,SRC="$wt/src",OUT="$OUT/$side/zygo23" \
        "$wt/scripts/run_zygo_regression.slurm")
    jobids+=("$j"); echo "$j	$side	zygo23" | tee -a "$OUT/jobs.tsv"
done

deps=$(IFS=:; echo "${jobids[*]}")
pairs=""
for g in $groups zygo23; do pairs+=" --pair $g $OUT/base/$g $OUT/cand/$g"; done
j=$(sbatch --parsable --partition="$PARTITION" --time=0:30:00 --mem=8gb --cpus-per-task=1 \
    --job-name=reg-diff --dependency="afterany:$deps" \
    --output="$REPO/logs/regression_%j.log" \
    --wrap="set -e; export PYTHONPATH=$CAND_WT/src; P=$REPO/.pixi/envs/default/bin/python; \
\$P $CAND_WT/scripts/regression_check.py diff --title 'regression $(basename "$CAND_WT") vs $(basename "$BASE_WT")' --out $OUT/diff $pairs; \
\$P $CAND_WT/scripts/check_record_selfcall.py --db $CAND_WT/db --reports $OUT/cand/*/runs --out $OUT/diff/record_selfcall_candidate.tsv || true")
echo "$j	diff	all" | tee -a "$OUT/jobs.tsv"
echo "summary will be written to $OUT/diff/summary.md"
