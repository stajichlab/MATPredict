#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 8gb -t 1:00:00 -J ss-bridge --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_supported_span/bridging_%j.log
set +eu; source /etc/profile; set -eu
export PYTHONPATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-supported/src PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-e9060da/.pixi/envs/test/bin:$PATH
cd /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-08_supported_span
python bridging.py /bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-supported/db bridging.tsv /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_basidiomycota_v060/reports_all.tar.zst /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_ascomycota_v060/reports_all.tar.zst /bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-03_mucoromycotina_mat/reports_all.tar.zst
