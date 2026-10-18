# Shared environment for the Coccidioides pilot (source me). CLONE is a clone of this branch on hpcc.
export CLONE=${CLONE:-/bigdata/stajichlab/jstajich/projects/MATPredict_cocci}
export TOOLS=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin   # python + exonerate, miniprot, diamond, hmmer
export PATH=$TOOLS:$PATH
export PYTHONPATH=$CLONE/src
export MATPREDICT_DB_ROOT=$CLONE/db
export MATPREDICT_CACHE_DIR=${MATPREDICT_CACHE_DIR:-/bigdata/stajichlab/jstajich/projects/MATPredict/.matpredict_cache}
export OUT=${OUT:-$CLONE/pilot_out}
export PILOT=$CLONE/datasets/reads/coccidioides/pilot/pilot_samples.tsv
export PANELS=$CLONE/panels/coccidioides
