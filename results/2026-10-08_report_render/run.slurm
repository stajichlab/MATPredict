#!/bin/bash
#SBATCH -p exfab -A exfab -c 4 --mem 16gb -t 2:00:00 -J mp-deepdive --array=0-3 --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_report_render/out/slurm_%A_%a.log
# Deep dive on 4 real campaign genomes: (1) old-format campaign YAML -> HTML+PDF; (2) fresh detect with main e9060da (timed, writes report.html) -> PDF.
set +eu; source /etc/profile; set -eu
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/src MATPREDICT_DB_ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/db MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/.pixi/envs/test/bin:$PATH
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_report_render
line=$(sed -n "$((SLURM_ARRAY_TASK_ID+2))p" genomes.tsv); n=$(echo "$line" | cut -f1); g=$(echo "$line" | cut -f2); t=$(echo "$line" | cut -f3); x=$(echo "$line" | cut -f4)
echo "$n on $(hostname)"
s=$(date +%s.%N); python -m MATPredict report genome --run inputs/$n --sample $n --out out/${n}_campaign.html --pdf out/${n}_campaign.pdf; echo -e "$n\told_yaml_html_pdf_s\t$(echo "$(date +%s.%N) - $s" | bc)" > out/$n.timing.tsv
s=$(date +%s.%N); python -m MATPredict detect --genome $g --taxid $t $x --sample $n --out-dir runs/$n > runs/$n.log 2>&1; rc=$?
echo -e "$n\tdetect_fresh_s\t$(echo "$(date +%s.%N) - $s" | bc)\texit $rc" >> out/$n.timing.tsv
s=$(date +%s.%N); python -m MATPredict report genome --run runs/$n --out out/${n}_fresh.html --pdf out/${n}_fresh.pdf; echo -e "$n\tfresh_html_pdf_s\t$(echo "$(date +%s.%N) - $s" | bc)" >> out/$n.timing.tsv
