#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J cochl-proto
#SBATCH -p short
#SBATCH --array=0-9%5
# Prototype db (GAP1/BGL1 roster genes) detect run. LIST (env) = tsv genome<TAB>taxid; OUTNAME (env) = result dir name; EXTRA (env) = extra detect args.
# One task = every 10th genome. Resumable.
set -uo pipefail
WD=/bigdata/stajichlab/jstajich/projects/cochl_recall
CODE=$WD/code
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=$ENV/bin:$PATH PYTHONPATH=$CODE/src MATPREDICT_DB_ROOT=${DBROOT:-$WD/proto_db} MATPREDICT_CACHE_DIR=$WD/cache_proto MATPREDICT_HTML=0 PYTHONNOUSERSITE=1
OUT=$WD/${OUTNAME:?}; mkdir -p $OUT
W=${SCRATCH:?}/cochl_proto_${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W
n=0
while IFS=$'\t' read -r g taxid; do
  if [ $(( n % 10 )) -eq $SLURM_ARRAY_TASK_ID ] && [ ! -s $OUT/$g/detection_report.yaml ]; then
    zcat $LIB/$g.fa.gz > $W/$g.fna
    mkdir -p $OUT/$g
    timeout 3000 python -m MATPredict detect --genome $W/$g.fna --taxid $taxid --out-dir $OUT/$g \
      --partial-core-policy coverage_profile ${EXTRA:-} \
      --evidence-diagnostics $OUT/$g/evidence_diagnostics.jsonl > $OUT/$g/stdout.log 2> $OUT/$g/stderr.log || echo "FAIL $g"
    rm -f $W/$g.fna
  fi
  n=$((n+1))
done < $WD/${LIST:?}
rm -rf $W
