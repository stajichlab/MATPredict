#!/bin/bash
# Top up the thin-evidence species (species.txt) to the cap of 700 runs. Prepared 2026-11-10 during the HPCC shutdown; run ONCE, on the HPCC login node,
# from the root of a frozen worktree of main (git fetch; git worktree add --detach .claude/worktrees/run-<sha> FETCH_HEAD; cd there).
#
#   bash datasets/reads/sra_survey/pilot/orchestration/thin_topup_2026-11-10/run.sh refresh   # 1. new availability snapshot (about 2 minutes)
#   bash datasets/reads/sra_survey/pilot/orchestration/thin_topup_2026-11-10/run.sh lists     # 2. run lists; read R/expand_summary.tsv
#   bash datasets/reads/sra_survey/pilot/orchestration/thin_topup_2026-11-10/run.sh submit    # 3. queue stage thin_topup per wave
#   bash datasets/reads/sra_survey/pilot/orchestration/thin_topup_2026-11-10/run.sh check     # 4. after the summaries: runs to rerun
#   5. in a checkout of a new branch, per wave: append_stage.py --wave-dir datasets/reads/sra_survey/pilot/<wave> --stage $W/<wave>/thin_topup/summary/results_all.tsv
#      (appends the new runs, decides the topped-up species again with their thin-evidence cap); rebuild the reads site; open a PR.
#
# Choices (2026-11-10): cap 700 as round 3; no per-study limit (--per-study 700: the curator had the runs the 25-per-study limit left out added, expand300b);
# SKETCH=1 SCREEN=1 and 8 GB as round 3; preempt first, no retry arrays (rerun what `check` lists on stajichlab).
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/MATPredict
REPO=$(pwd)
HERE=$REPO/datasets/reads/sra_survey/pilot/orchestration/thin_topup_2026-11-10
W=$P/results/sra_reads_typing/waves
R=${R:-$P/results/2026-11-10_thin_topup}
PY=$P/.pixi/envs/default/bin/python
DRV=$REPO/datasets/reads/sra_survey/pilot/wave_driver.py
[ -f "$DRV" ] || { echo "run from the root of a MATPredict checkout" >&2; exit 2; }
mkdir -p "$R"
case "${1:-}" in
  refresh)
    bash "$REPO/datasets/reads/availability/refresh.sh" 2>&1 | tee "$R/refresh.log" ;;
  lists)
    SNAP=$(ls -d "$REPO"/datasets/reads/availability/snapshots/20* | sort | tail -1)
    echo "snapshot $SNAP" | tee "$R/snapshot.txt"
    "$PY" "$REPO/datasets/reads/availability/expand_runs.py" --snapshot "$SNAP" --waves-root "$W" --out "$R" --cap 700 --per-study 700 --species-file "$HERE/species.txt"
    column -t -s$'\t' "$R/expand_summary.tsv" ;;
  submit)
    [ -e "$R/submitted" ] && { echo "already submitted:"; cat "$R/submitted"; exit 1; }
    for f in "$R"/runs_wave*.tsv; do
      w=$(basename "$f" .tsv); w=${w#runs_}
      "$PY" "$DRV" submit-list --wave "$W/$w" --runs "$f" --stage-name thin_topup --parts preempt:preempt,stajichlab:stajichlab --retry-part none \
        --mem 8gb --env SKETCH=1 --env SCREEN=1 | tee -a "$R/submitted"
    done ;;
  check)
    for f in "$R"/runs_wave*.tsv; do w=$(basename "$f" .tsv); "$PY" "$DRV" incomplete --wave "$W/${w#runs_}" --stage-name thin_topup; done ;;
  *) sed -n 2,13p "$0"; exit 2 ;;
esac
