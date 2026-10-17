# sourced by every job
R=/bigdata/stajichlab/jstajich/projects/MATPredict_cocci/reads_asm
REPO=/bigdata/stajichlab/jstajich/projects/MATPredict_cocci
S=${SCRATCH:-/scratch/$USER}/reads_asm
mkdir -p $S
for m in diamond/2.1.12 bwa-mem2/2.3 spades/4.2.0 megahit/1.2.9 minimap2/2.30 ncbi-blast/2.17.0+ samtools/1.22.1 seqtk/1.4 miniprot/0.2 exonerate/2.4.0; do module load $m >/dev/null 2>&1 || echo "module load $m failed"; done
module load aria2/1.36.0 >/dev/null 2>&1
module load apptainer >/dev/null 2>&1
SIF=/bigdata/stajichlab/shared/singularity_cache/AAFTF-0.7.0-beta.4-18-g9ce56a1.sif   # $NXF_APPTAINER_CACHEDIR; AAFTF 0.7.0-beta.4, bwa 0.7.19, samtools 1.24, SPAdes 4.3.0 inside
AP="apptainer exec --env MPLCONFIGDIR=/tmp -B $R -B $S -B $REPO $SIF"
export PATH=/rhome/jstajich/bin:$PATH
# fetch RUN: download both mates into $S/RUN, verify md5
fetch() {
  local run=$1; mkdir -p $S/$run; cd $S/$run
  read u1 u2 m1 m2 < <(awk -F'\t' -v r=$run '$1==r{print $5,$6,$7,$8}' $R/runs.tsv)
  for pair in "1 $u1 $m1" "2 $u2 $m2"; do set -- $pair
    f=${run}_$1.fastq.gz
    if [ ! -s $f ] || [ "$(md5sum $f | cut -d' ' -f1)" != "$3" ]; then
      rm -f $f $f.aria2
      for try in 1 2 3 4; do
        aria2c -q -x 8 -s 8 -k 4M --retry-wait=5 --max-tries=5 -o $f https://$2 && break
        echo "retry $try $f"; sleep 5
      done
    fi
    [ "$(md5sum $f | cut -d' ' -f1)" == "$3" ] && echo "md5 ok $f" || { echo "MD5 FAIL $f"; return 1; }
  done
}
# subsample RUN NPAIRS_TARGET: every Nth pair (N from read_count), writes $S/RUN/sub_{1,2}.fastq.gz
subsample() {
  local run=$1 target=$2 rc N; rc=$(awk -F'\t' -v r=$run '$1==r{print $4}' $R/runs.tsv)
  N=$(( rc / target )); [ $N -lt 1 ] && N=1
  for m in 1 2; do
    zcat $S/$run/${run}_$m.fastq.gz | awk -v N=$N '{ if (int((NR-1)/4) % N == 0) print }' | gzip -1 > $S/$run/sub_$m.fastq.gz
  done
  echo "subsample $run: every ${N}th pair of $rc"
}
