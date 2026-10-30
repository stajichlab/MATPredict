#!/bin/bash
# Fetch a window around each replay MTL/MATsc call (coordinates from results/2026-10-29_near_tie/flank_synteny.tsv)
# from NCBI nuccore. Window = call +/- 25 kb. Writes data/<assembly>_<contig>_<start>-<stop>.fa
cd "$(dirname "$0")"; mkdir -p data
while IFS=$'\t' read -r asm contig start end; do
  s=$((start-25000)); [ $s -lt 1 ] && s=1; e=$((end+25000))
  out="data/${asm}_${contig}_${s}-${e}.fa"
  [ -s "$out" ] || curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${contig}&rettype=fasta&seq_start=${s}&seq_stop=${e}" > "$out"
  echo "$out $(grep -v '>' "$out" | tr -d '\n' | wc -c)"; sleep 0.4
done <<'T'
GCA_023627795.1	JAMAIE010000005.1	614225	630853
GCA_003706655.3	PPLX03000005.1	120524	137156
GCA_003123585.1	QBLI01000035.1	9748	26373
GCA_003123885.1	QBLJ01000026.1	125869	142498
GCA_049863505.1	CM113344.1	1130885	1139023
T
