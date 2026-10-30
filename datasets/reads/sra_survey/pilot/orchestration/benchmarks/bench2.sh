#!/bin/bash
#SBATCH -p exfab -A exfab -c 4 --mem 8gb -t 0:40:00 -J bench-sra2
set -uo pipefail
cd "$SCRATCH" 2>/dev/null || cd /tmp
echo host=$(hostname) date=$(date -Is)
t() { local s=$(date +%s.%N); "$@"; local rc=$?; local e=$(date +%s.%N); printf "  secs=%.1f rc=%s\n" $(echo "$e - $s" | bc) $rc; }
ena() { local run=$1 n=${#1}; local d=""; [ $n -gt 9 ] && d="/0$(printf %02d $((10#${run: -2})))" ; echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${run:0:6}${d}/${run}/${run}_1.fastq.gz"; }
run=SRR1587944
echo "== ENA url for $run: $(ena $run)"
echo "-- ENA 60 s cap"; t bash -c "curl -sS -L -m 60 -w '%{http_code} ' $(ena $run) | wc -c"
echo "-- NCBI FTP .sra 60 s cap"; t bash -c "curl -sS -L -m 60 -w '%{http_code} %{url_effective}\n' -o /dev/null https://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR158/$run/$run.sra"
echo "-- NCBI FTP bytes in 30 s"; t bash -c "curl -sS -L -m 30 https://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR158/$run/$run.sra | wc -c"
echo "-- ENA bytes in 30 s"; t bash -c "curl -sS -L -m 30 $(ena $run) | wc -c"
echo "-- AWS 4 parallel streams, 20 s cap each"; s=$(date +%s.%N); for i in 1 2 3 4; do (curl -sS -m 20 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR28962295/SRR28962295 | wc -c > par$i.txt) & done; wait; e=$(date +%s.%N); echo "  total bytes $(cat par*.txt | paste -sd+ | bc) in $(echo "$e - $s" | bc) s"; rm -f par*.txt
module load sratoolkit/3.2.0 2>&1 | head -2; which fastq-dump fasterq-dump
echo "-- fastq-dump local on an AWS-downloaded .sra"
t curl -sS -m 120 -o $run.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run; ls -l $run.sra
t bash -c "fastq-dump -X 500000 -Z --split-spot --ncbi_error_report never $run.sra 2>&1 | head -c 400 | head -2; fastq-dump -X 500000 -Z --split-spot --ncbi_error_report never $run.sra 2>/dev/null | wc -l"
echo "-- windows: spots 5,000,000-5,250,000 and 100,000-350,000"; t bash -c "fastq-dump -N 5000000 -X 5250000 -Z --split-spot --ncbi_error_report never $run.sra 2>/dev/null | wc -l"
rm -f $run.sra; echo done
