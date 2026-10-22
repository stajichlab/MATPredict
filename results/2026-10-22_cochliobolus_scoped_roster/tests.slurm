#!/usr/bin/bash -l
#SBATCH -N 1 -c 4 --mem 16gb --time 02:00:00 -J cochl-tests -p short
set -uo pipefail
CODE=${CODE:?}
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
export PATH=$ENV/bin:$PATH PYTHONPATH=$CODE/src PYTHONNOUSERSITE=1 MATPREDICT_HTML=0
export MATPREDICT_CACHE_DIR=/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache
cd $CODE
python -c 'import MATPredict;print(MATPredict.__file__)'
git log --oneline -1
python -m pytest -q -p no:cacheprovider ${TESTS:-tests/db tests/detect} --no-header -rf 2>&1 | tail -30 || true
