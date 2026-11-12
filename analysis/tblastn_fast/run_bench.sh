#!/usr/bin/bash
# Compare tblastn default vs -task tblastn-fast: wall time and HSPs, per genome.
# Usage: run_bench.sh OUTDIR   (needs SCRATCH; run on one node, 4 threads)
set -eu
module load ncbi-blast/2.17.0+
ROOT=/bigdata/stajichlab/jstajich/projects/MATPredict
OUT=${1:?outdir}; W=${SCRATCH:?}/tbfast; mkdir -p "$OUT" "$W"
LIB=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
PY=$ROOT/.pixi/envs/default/bin/python
PYTHONPATH=$ROOT/src $PY -c "
from pathlib import Path
from MATPredict.detect.reference_fasta import build_reference_fasta
build_reference_fasta(Path('$ROOT/db'), Path('$W/ref.faa'))"
grep -c '>' $W/ref.faa
FMT='6 qseqid sseqid pident length qlen qcovhsp evalue bitscore sstart send sframe'
$PY - <<P > $W/asm.txt
import json
d=json.load(open('$ROOT/results/truth_assemblies.json'))
print('\n'.join(sorted({v['bfd_asmid'] for v in d.values() if v['bfd_asmid']})))
P
echo -e "asm\tmode\tseconds\thsps" > $OUT/timing.tsv
while read a; do
  zcat $LIB/$a.fa.gz > $W/g.fa
  makeblastdb -in $W/g.fa -dbtype nucl -out $W/db -logfile /dev/null
  for rep in 1 2; do
   for mode in tblastn tblastn-fast; do
    s=$(date +%s.%N)
    tblastn -task $mode -query $W/ref.faa -db $W/db -seg no -num_threads 4 -outfmt "$FMT" > $W/o.tsv
    e=$(date +%s.%N)
    echo -e "$a\t$mode\t$(echo "$e - $s" | bc)\t$(wc -l < $W/o.tsv)" >> $OUT/timing.tsv
    [ $rep = 1 ] && zstd -q -f $W/o.tsv -o $OUT/$a.$mode.tsv.zst
   done
  done
done < $W/asm.txt
