#!/bin/bash
# HPCC: windows, then default and relaxed miniprot of refs.faa on each.
export PATH=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin:$PATH
cd "$(dirname "$0")"
python3 windows.py win
for f in win/*.fa; do b=${f%.fa}; miniprot -t4 --gff -P MP "$f" refs.faa > $b.miniprot.gff 2> $b.miniprot.log; miniprot -t4 -k5 -l2 --outs=0.3 --outc=0.2 -N 5 --gff -P MR "$f" refs.faa > $b.relaxed.gff 2>> $b.miniprot.log; done
ls win | wc -l
