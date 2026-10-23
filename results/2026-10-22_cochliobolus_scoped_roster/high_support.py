"""Classify every `high` Ascomycota:MAT call in a variant: 'supported' when a SLA2/APN2/COX13 flank is
modelled (status polished_*) over >= 50% of its reference; otherwise 'tie_break' (high only through an
unmodelled or short fragment). usage: high_support.py RUNS/VARIANT GENOME_LIST > out.tsv"""
import sys, yaml
from pathlib import Path
print('genome\tcontig\tidiomorph\tclass\tflank_evidence')
for g in open(sys.argv[2]).read().split():
    p = Path(sys.argv[1]) / g / 'detection_report.yaml'
    if not p.exists() or not p.stat().st_size: continue
    d = yaml.load(p.read_text(), Loader=yaml.CSafeLoader)
    for r in d.get('detected') or []:
        if r.get('family') != 'Ascomycota:MAT' or r.get('confidence') != 'high': continue
        fl = [e for e in r.get('gene_evidence') or [] if e.get('gene') in ('SLA2', 'APN2', 'COX13')]
        ok = any(str(e.get('status', '')).startswith('polished') and (e.get('coverage') or 0) >= 50 for e in fl)
        desc = ';'.join(f"{e['gene']}:{e.get('identity')}/{e.get('coverage')}/{e.get('status')}" for e in fl) or 'none'
        print(f"{g}\t{r['contig']}\t{r['idiomorph']}\t{'supported' if ok else 'tie_break'}\t{desc}")
