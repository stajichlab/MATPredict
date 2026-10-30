#!/bin/bash
#SBATCH --job-name=v070-pack
#SBATCH --partition=highclock
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/pack_%j.log
set -euo pipefail
R=/bigdata/stajichlab/jstajich/projects/MATPredict/results
export PATH=$(dirname $(ls /bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/*/bin/zstd 2>/dev/null | head -1)):$PATH
cd $R/2026-10-26_v070_ascomycota && find wave_* \( -name detection_report.yaml -o -name wall_seconds \) | sort | tar -cf - -T - | zstd -19 -T8 -q -o reports_all.tar.zst
cd $R/2026-10-26_v070_basidiomycota && find runs \( -name detection_report.yaml -o -name wall_seconds \) | sort | tar -cf - -T - | zstd -19 -T8 -q -o reports_all.tar.zst
cd $R/2026-10-26_v070_mucoromycotina && find runs \( -name detection_report.yaml -o -name wall_seconds \) | sort | tar -cf - -T - | zstd -19 -T8 -q -o reports_all.tar.zst
ls -l $R/2026-10-26_v070_*/reports_all.tar.zst
