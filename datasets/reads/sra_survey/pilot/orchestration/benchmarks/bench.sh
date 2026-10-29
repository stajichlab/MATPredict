#!/bin/bash
#SBATCH -p exfab -A exfab -c 4 --mem 8gb -t 0:40:00 -J bench-sra
# Throughput of four ways to get reads of one run onto a compute node (60 s cap each, partial reads only; nothing is kept).
set -uo pipefail
cd "$SCRATCH" 2>/dev/null || cd /tmp
echo host=$(hostname) date=$(date -Is)
t() { local s=$(date +%s.%N); "$@"; local rc=$?; local e=$(date +%s.%N); printf "  secs=%.1f rc=%s\n" $(echo "$e - $s" | bc) $rc; }
for run in SRR5087068 SRR1587944 SRR28962295; do
  p6=${run:0:6}; pre=${run:0:3}
  echo "== $run"
  echo "-- ENA https (60 s cap)"; t bash -c "curl -sS -m 60 https://ftp.sra.ebi.ac.uk/vol1/fastq/${p6}/$( [ ${#run} -eq 10 ] && echo 0$(echo ${run:9:1}) )/${run}/${run}_1.fastq.gz 2>/dev/null | wc -c" 2>&1 | head -3
  echo "-- NCBI FTP sra-instant .sra (60 s cap)"; t bash -c "curl -sS -m 60 https://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/${pre}/${p6}/${run}/${run}.sra | wc -c"
  echo "-- AWS open-data .sra (60 s cap)"; t bash -c "curl -sS -m 60 https://sra-pub-run-odp.s3.amazonaws.com/sra/${run}/${run} | wc -c"
done
echo "-- sratoolkit remote random access, first 500000 spots"
module load sratoolkit/3.2.0 2>&1 | head -2
export NCBI_SETTINGS=/dev/null
for run in SRR5087068 SRR1587944; do
  echo "== $run"; t bash -c "timeout 120 fastq-dump -X 500000 -Z --split-spot --ncbi_error_report never $run 2>&1 | wc -c"
done
echo done
