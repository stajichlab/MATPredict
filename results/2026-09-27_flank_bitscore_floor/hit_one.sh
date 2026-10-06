#!/bin/bash
# tblastn (detect settings: -seg no, evalue 10) of every core protein of the families
# flagged for this genome, genome-wide, with the report's genetic code. Records genome length.
set -euo pipefail
g=$1; rd=$2; od=$3; fams=$4
E=/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin
lib=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
[ -s "$od/$g.tsv.zst" ] && exit 0
w=$TMPD/$g; mkdir -p $w
zcat $lib/$g.fa.gz > $w/g.fna
grep -v '>' $w/g.fna | tr -d '\n' | wc -c > $od/$g.len
gc=$(python3 -c "import yaml;print(yaml.safe_load(open('$rd/detection_report.yaml')).get('genetic_code') or 1)")
python3 - "$rd/_reference.faa" "$w/core.faa" "$fams" <<'PY'
import sys
FL={"SLA2","APN2","COX13","sla2","apn2","cox13","PAP1","OBP1","PIK1","tptA","rnhA","glrA","algA","btbA","MIP1","beta_fg","RibL6","BAP31","CAF1","STE20"}
tags=[f"_{f.split(':')[1]}_" for f in sys.argv[3].split(',')]
out=open(sys.argv[2],'w'); keep=False
for l in open(sys.argv[1]):
    if l.startswith('>'):
        h=l[1:].strip(); rec=h.split('|')[0]; gene=h.split('|')[-1]
        keep=(gene not in FL) and any(t in rec+'_' for t in tags)
    if keep: out.write(l)
PY
[ -s $w/core.faa ] || { echo "no core for $g $fams" >&2; rm -rf $w; exit 0; }
$E/makeblastdb -in $w/g.fna -dbtype nucl -out $w/db >/dev/null
$E/tblastn -query $w/core.faa -db $w/db -db_gencode $gc -evalue 10 -seg no -num_threads 1 -max_target_seqs 200 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen' | zstd -q -o $od/$g.tsv.zst
rm -rf $w
