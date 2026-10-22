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
# The pixi env carries an editable MATPredict install (mat-prof/src) that can shadow PYTHONPATH:
# put this worktree's src first and refuse to run if another checkout is imported.
python - ${TESTS:-tests/db tests/detect} <<'PY'
import os, sys
sys.path.insert(0, os.path.join(os.getcwd(), "src"))
import MATPredict
assert MATPredict.__file__.startswith(os.path.join(os.getcwd(), "src")), MATPredict.__file__
print("testing", MATPredict.__file__)
import pytest
sys.exit(pytest.main(["-q", "-p", "no:cacheprovider", "--no-header", "-rf", *sys.argv[1:]]))
PY
