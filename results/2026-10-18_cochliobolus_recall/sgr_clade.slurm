#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J cochl-sgr
#SBATCH -p short
#SBATCH --array=0-7%6
# matpredict detect on each clade genome with --strong-single-gene-policy on (+ the adjacent-gene prototype refs added when ADJ=1).
# One task = every 8th genome of clade_genome_names_clade.txt. Resumable.
set -uo pipefail
WD=/bigdata/stajichlab/jstajich/projects/cochl_recall
CODE=$WD/code
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=$ENV/bin:$PATH PYTHONPATH=$CODE/src MATPREDICT_DB_ROOT=$CODE/db MATPREDICT_CACHE_DIR=$WD/cache MATPREDICT_HTML=0 PYTHONNOUSERSITE=1
OUT=$WD/sgr_clade; mkdir -p $OUT $WD/cache
W=${SCRATCH:?}/cochl_sgr_${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W
n=0
while IFS=$'\t' read -r g taxid; do
  if [ $(( n % 8 )) -eq $SLURM_ARRAY_TASK_ID ] && [ ! -s $OUT/$g/detection_report.yaml ]; then
    zcat $LIB/$g.fa.gz > $W/$g.fna
    mkdir -p $OUT/$g
    timeout 3000 python -m MATPredict detect --genome $W/$g.fna --taxid $taxid --out-dir $OUT/$g \
      --partial-core-policy coverage_profile --strong-single-gene-policy on \
      --evidence-diagnostics $OUT/$g/evidence_diagnostics.jsonl > $OUT/$g/stdout.log 2> $OUT/$g/stderr.log || echo "FAIL $g"
    rm -f $W/$g.fna
  fi
  n=$((n+1))
done < $WD/clade_names_taxid.tsv
rm -rf $W
