#!/bin/bash
#SBATCH --job-name=v070-tables
#SBATCH --partition=highclock
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/tables_%x_%j.log
# usage: sbatch --export=ALL,CAMP=asco|basidio|mucoro tables_v070.slurm
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/MATPredict
R=$P/results
SAMPLES=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv
PY=$P/.claude/worktrees/run-b317ece
cd $PY
PYBIN=$(ls $P/.pixi/envs/default/bin/python)
export PYTHONPATH=$PY/src
case $CAMP in
 asco)    C=$R/2026-10-26_v070_ascomycota
          $PYBIN scripts/build_cleanrun_tables.py --campaign $C --lists $C/lists/wave_*.tsv --baseline $R/2026-10-03_ascomycota_v060/reports_all.tar.zst --samples $SAMPLES --out $C/tables --jobs 8 ;;
 basidio) C=$R/2026-10-26_v070_basidiomycota
          $PYBIN scripts/build_cleanrun_tables.py --campaign $C --baseline $R/2026-10-03_basidiomycota_v060/reports_all.tar.zst --samples $SAMPLES --out $C/tables --jobs 8 ;;
 mucoro)  C=$R/2026-10-26_v070_mucoromycotina
          $PYBIN scripts/build_cleanrun_tables.py --layout mucoro --campaign $C --baseline-calls $R/2026-10-03_mucoromycotina_mat/calls.tsv --baseline $R/2026-10-03_mucoromycotina_mat/reports_all.tar.zst --samples $SAMPLES --out $C/tables --jobs 8 ;;
esac
cp $C/MANIFEST.txt $C/tables/MANIFEST.txt
