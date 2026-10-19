#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J cochl-confirm
#SBATCH -p short,highclock,exfab
# Independent confirmation on the four C4/C5 genomes + detect with --strong-single-gene-policy on (analysis only).
set -uo pipefail
WD=/bigdata/stajichlab/jstajich/projects/cochl_recall
CODE=$WD/code
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
PFAM=/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/run-c9c0ec6/results/2026-10-18_dothideomycetes_audit/pfam_subset.hmm
export PATH=$ENV/bin:$PATH
export PYTHONPATH=$CODE/src
export MATPREDICT_DB_ROOT=$CODE/db
export DBROOT=$CODE/db
export MATPREDICT_CACHE_DIR=$WD/cache
export MATPREDICT_HTML=0 PYTHONNOUSERSITE=1
mkdir -p $WD/cache $WD/confirm $WD/sgr
W=${SCRATCH:?}/cochl_confirm_$SLURM_JOB_ID; mkdir -p $W
for asm in GCA_028891585.1_CocheC5_4m GCF_028858645.1_CocheC4_6 GCF_000354255.1_CocheC4_1 GCA_000338975.1_CocheC5_3; do
  zcat $LIB/$asm.fa.gz > $W/$asm.fna
  echo "=== $asm"
  python $CODE/../confirm_locus.py $W/$asm.fna $WD/confirm/$asm $PFAM 2>&1 | tee $WD/confirm/$asm.log
  taxid=5016; mkdir -p $WD/sgr/$asm
  python -m MATPredict detect --genome $W/$asm.fna --taxid $taxid --out-dir $WD/sgr/$asm \
     --partial-core-policy coverage_profile --strong-single-gene-policy on \
     --evidence-diagnostics $WD/sgr/$asm/evidence_diagnostics.jsonl > $WD/sgr/$asm.stdout 2> $WD/sgr/$asm.stderr
  echo "detect rc=$?"
done
rm -rf $W
