#!/usr/bin/bash -l
#SBATCH -N 1 -c 2 --mem 8gb --time 02:00:00 -J cochl-flank
#SBATCH -p short
#SBATCH --array=0-7%6
# one task = every 8th genome of clade_genome_names.txt (offset = array index)
set -uo pipefail
WD=/bigdata/stajichlab/jstajich/projects/cochl_recall
ENV=/bigdata/stajichlab/jstajich/projects/mat-prof/.pixi/envs/test
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
export PATH=$ENV/bin:$PATH
W=${SCRATCH:?}/cochl_flank_${SLURM_JOB_ID}_$SLURM_ARRAY_TASK_ID; mkdir -p $W $WD/flank
n=0
while read g; do
  if [ $(( n % 8 )) -eq $SLURM_ARRAY_TASK_ID ]; then
    if [ ! -s $WD/flank/$g.tsv ]; then
      zcat $LIB/$g.fa.gz > $W/$g.fna
      python $WD/flank_map.py $W/$g.fna $WD/_reference.faa $WD/flank/$g.tsv $g || echo "FAIL $g"
      rm -f $W/$g.fna
    fi
  fi
  n=$((n+1))
done < $WD/clade_genome_names_clade.txt
rm -rf $W
