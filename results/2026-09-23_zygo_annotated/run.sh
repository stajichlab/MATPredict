#!/usr/bin/bash
set -u
R=/bigdata/stajichlab/jstajich/projects/MATPredict
E=$R/.pixi/envs/default
SRC=$R/.claude/worktrees/polish-scope-cuts/src
[[ -f "$SRC/MATPredict/__init__.py" ]] || { echo "SRC missing on $(hostname)" >&2; exit 1; }
export PATH="$E/bin:$PATH"; export PYTHONPATH="$SRC"
A=/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/annotate
D=$R/results/2026-09-23_zygo_annotated
W="${SCRATCH:?}/zygo_annot"; mkdir -p "$W/fasta"

# Step 1: the annotation converter -- the proteome path REQUIRES deflines
# carrying 'contig:start-end:strand', which the raw .proteins.fa lacks.
mapfile -t ORGS < <(cut -f1 $R/results/2026-09-23_zygo_bar/zygo_truth.tsv | sort -u)
for org in "${ORGS[@]}"; do
  gbk=$(ls "$A/$org/annotate_results/"*.gbk 2>/dev/null | head -1)
  [[ -s "$gbk" ]] && echo "$gbk"
done > "$W/gbk.list"
echo "converting $(wc -l < "$W/gbk.list") annotations"
xargs -a "$W/gbk.list" "$E/bin/python" "$R/.claude/worktrees/polish-scope-cuts/scripts/convert_annotations.py" --out-dir "$W/fasta" || true
ls "$W/fasta" | head -4

# Step 2: detect, annotated path
run_one(){
  org="$1"
  stem=$(basename "$(ls $A/$org/annotate_results/*.gbk 2>/dev/null | head -1)" .gbk)
  fna="$W/fasta/$stem.fna"; faa="$W/fasta/$stem.faa"
  [[ -s "$fna" && -s "$faa" ]] || { echo "SKIP $org (no converted pair)" >&2; return 0; }
  out="$D/runs/$org"; [[ -s "$out/detection_report.yaml" ]] && return 0
  mkdir -p "$out"
  timeout 3600 "$E/bin/python" -m MATPredict detect \
     --genome "$fna" --proteins "$faa" --phylum Mucoromycota --out-dir "$out" \
     > "$out/stdout.log" 2> "$out/stderr.log"
  rc=$?; [[ $rc -ne 0 ]] && echo "FAIL $org rc=$rc: $(tail -1 $out/stderr.log)" >&2; return 0; }
export -f run_one; export A D E W SRC
printf '%s\n' "${ORGS[@]}" | xargs -P 6 -L1 bash -c 'run_one "$0"'
echo "annotated zygo done: $(ls $D/runs/*/detection_report.yaml 2>/dev/null | wc -l)"
