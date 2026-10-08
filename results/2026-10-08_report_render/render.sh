#!/bin/bash
# Render each staged real campaign report to HTML and PDF with the pinned environment (WeasyPrint 69); time each.
# usage: render.sh   (run after: pixi install -e test --frozen in /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da)
set -u
PY=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/.pixi/envs/test/bin/python
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_report_render
mkdir -p out; : > render_timing.tsv
$PY -c "import weasyprint,sys; print('weasyprint',weasyprint.__version__)" | tee weasyprint_version.txt
echo -e "report\thtml_s\tpdf_s\thtml_bytes\tpdf_bytes" >> render_timing.tsv
for d in inputs/*/; do
  n=$(basename $d)
  s=$(date +%s.%N); $PY -m MATPredict report genome --run $d --sample $n --out out/$n.html > out/$n.html.log 2>&1; e=$(date +%s.%N)
  s2=$(date +%s.%N); $PY -m MATPredict report genome --run $d --sample $n --out out/$n.html --pdf out/$n.pdf > out/$n.pdf.log 2>&1; e2=$(date +%s.%N)
  echo -e "$n\t$(echo "$e - $s" | bc)\t$(echo "$e2 - $s2" | bc)\t$(stat -c %s out/$n.html)\t$(stat -c %s out/$n.pdf 2>/dev/null || echo 0)" >> render_timing.tsv
done
column -t -s$'\t' render_timing.tsv
