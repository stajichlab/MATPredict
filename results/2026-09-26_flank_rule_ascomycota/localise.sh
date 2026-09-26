#!/bin/bash
# For each genome with a changed call: tblastn the family's CORE reference proteins
# (from the run's _reference.faa) genome-wide with the report's genetic code.
# usage: localise.sh GENOME RUNDIR OUTDIR   (run with $TMPD set to a scratch dir)
set -euo pipefail
g=$1; rd=$2; od=$3
lib=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes
[ -s "$od/$g.tsv.zst" ] && exit 0
w=$TMPD/$g; mkdir -p $w
zcat $lib/$g.fa.gz > $w/g.fna
gc=$(python3 -c "import yaml;print(yaml.safe_load(open('$rd/detection_report.yaml')).get('genetic_code') or 1)")
python3 - "$rd/_reference.faa" "$w/core.faa" <<'PY'
import sys,re
flank=re.compile(r'^(SLA2|APN2|COX13|sla2|apn2|cox13|PAP1|OBP1|PIK1|tptA|rnhA|glrA|algA|CIMG_.*)$')
out=open(sys.argv[2],'w'); keep=False
for l in open(sys.argv[1]):
    if l.startswith('>'): keep=not flank.match(l.strip().split('|')[-1])
    if keep: out.write(l)
PY
makeblastdb -in $w/g.fna -dbtype nucl -out $w/db >/dev/null
tblastn -query $w/core.faa -db $w/db -db_gencode $gc -evalue 1e-3 -num_threads 2 -max_target_seqs 50 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen' | zstd -q -o $od/$g.tsv.zst
rm -rf $w
