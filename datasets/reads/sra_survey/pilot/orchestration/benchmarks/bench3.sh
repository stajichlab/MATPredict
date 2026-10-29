#!/bin/bash -l
#SBATCH -p exfab -A exfab -c 4 --mem 8gb -t 0:30:00 -J bench-sra3
cd "$SCRATCH" 2>/dev/null || cd /tmp
echo host=$(hostname) date=$(date -Is)
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/sratoolkit/3.2.0/bin:$PATH; which fastq-dump fasterq-dump
t() { local s=$(date +%s.%N); "$@"; local rc=$?; local e=$(date +%s.%N); printf "  secs=%.1f rc=%s\n" $(echo "$e - $s" | bc) $rc; }
run=SRR1587944
t curl -sS -m 120 -o $run.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/$run/$run; ls -l $run.sra
export VDB_CONFIG=/dev/null
echo "-- first 500000 spots, split-spot to stdout"
t bash -c "fastq-dump -X 500000 -Z --split-spot --ncbi_error_report never $run.sra 2>err1.txt | wc -l"; head -c 300 err1.txt
echo "-- window 5,000,000-5,500,000"
t bash -c "fastq-dump -N 5000000 -X 5500000 -Z --split-spot --ncbi_error_report never $run.sra 2>err2.txt | wc -l"; head -c 300 err2.txt
echo "-- format check (first record, 8 lines)"; fastq-dump -X 2 -Z --split-spot --ncbi_error_report never $run.sra 2>/dev/null | head -8
rm -f $run.sra; echo done
