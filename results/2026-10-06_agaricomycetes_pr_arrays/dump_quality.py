#!/usr/bin/env python3
"""Run on HPCC (/usr/bin/python3.12 has pyarrow): dump assembly quality for all BFD genomes.
Usage: dump_quality.py OUT.tsv   -> ASMID, complete_pct (BUSCO fungi_odb12), N50_bp, total_length_bp, contig_count"""
import sys
import pyarrow.parquet as pq
T = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/tables/"
b = {r["ASMID"]: r["complete_pct"] for r in pq.read_table(T + "busco_genome.parquet", columns=["ASMID", "complete_pct"]).to_pylist()}
with open(sys.argv[1], "w") as fo:
    fo.write("genome\tbusco_complete_pct\tn50_bp\ttotal_length_bp\tcontig_count\n")
    for r in pq.read_table(T + "asm_stats.parquet", columns=["ASMID", "N50_bp", "total_length_bp", "contig_count"]).to_pylist():
        fo.write(f"{r['ASMID']}\t{b.get(r['ASMID'], '')}\t{r['N50_bp']}\t{r['total_length_bp']}\t{r['contig_count']}\n")
