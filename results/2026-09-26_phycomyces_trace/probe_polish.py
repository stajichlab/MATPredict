"""Call detect's polish functions directly on the NRRL_1554 locus windows.

Usage: probe_polish.py <genome.fa> <contig> <start> <end> <reference.faa>
"""
import sys
from pathlib import Path
from MATPredict.detect.family_registry import load_all_families, load_record_families
from MATPredict.detect.search import polish_with_miniprot, polish_with_exonerate
import os
db = Path(os.environ["MATPREDICT_DB_ROOT"])
fam = [f for f in load_all_families(db) if str(f.key) .startswith("Mucoromycota") or getattr(f.key, "phylum", "") == "Mucoromycota"]
fam = [f for f in fam if "MAT" in str(f.key)][0]
rec = load_record_families(db)
genome, contig, s, e, ref = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
for g in ("sexP", "rnhA"):
    for name, fn in (("miniprot", polish_with_miniprot), ("exonerate", polish_with_exonerate)):
        m = fn(genome_fasta=Path(genome), family=fam, gene_name=g, reference_fasta=Path(ref),
               record_families=rec, window=(contig, s, e), genetic_code=1)
        print(g, name, None if m is None else (m.contig, m.start, m.end, getattr(m, "record_id", None), getattr(m, "identity", None)))
