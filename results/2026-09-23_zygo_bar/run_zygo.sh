#!/usr/bin/bash
set -u
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default
SRC=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/src
[[ -f "$SRC/MATPredict/__init__.py" ]] || { echo "SRC missing on $(hostname)" >&2; exit 1; }
export PATH="$E/bin:$PATH"
D=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-23_zygo_bar
run_one(){ org="$1"; fna="$2"
  out="$D/runs/$org"; [[ -s "$out/detection_report.yaml" ]] && return 0
  mkdir -p "$out"
  PYTHONPATH="$SRC" timeout 3600 "$E/bin/python" -m MATPredict detect \
     --genome "$fna" --phylum Mucoromycota --out-dir "$out" \
     > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; [[ $rc -ne 0 ]] && echo "FAIL $org rc=$rc: $(tail -1 $out/stderr.log)" >&2; return 0; }
export -f run_one; export D E SRC
awk -F'\t' '{print $1" "$6}' "$D/zygo_truth.tsv" | sort -u | xargs -P 6 -L1 bash -c 'run_one "$0" "$1"'
echo "zygo done: $(ls $D/runs/*/detection_report.yaml 2>/dev/null | wc -l)"
