#!/bin/bash
#SBATCH -p exfab -A exfab -c 2 --mem 2gb -t 0:10:00 -J probe-net
echo host=$(hostname); for t in curl pigz gzip awk gawk mawk md5sum python3 bc; do echo "$t: $(command -v $t || echo MISSING)"; done
awk --version 2>&1 | head -1
curl -sS -m 20 -o /dev/null -w "ena_http=%{http_code} t=%{time_total}\n" -I https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR058/SRR058918/SRR058918_1.fastq.gz
s=$(date +%s.%N); curl -sS -m 60 https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR058/SRR058918/SRR058918_1.fastq.gz | head -c 200000000 | wc -c; e=$(date +%s.%N); echo "200MB-ish in $(echo "$e - $s" | bc) s"
echo $SCRATCH; df -h $SCRATCH | tail -1
