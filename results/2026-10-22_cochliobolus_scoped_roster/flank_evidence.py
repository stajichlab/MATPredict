"""Every GAP1/BGL1 gene_evidence row in detected loci of a variant: identity, coverage, distance to nearest core gene.
usage: flank_evidence.py RUNS/VARIANT > out.tsv"""
import sys
from pathlib import Path
import yaml
print('genome\tcontig\tlocus_idiomorph\tconfidence\tgene\tstatus\tidentity\tcoverage\tdist_to_core\tcore_best_identity\tcore_best_cov')
for d in sorted(Path(sys.argv[1]).iterdir()):
    p = d / 'detection_report.yaml'
    if not p.exists() or not p.stat().st_size: continue
    doc = yaml.load(p.read_text(), Loader=yaml.CSafeLoader)
    for r in doc.get('detected') or []:
        ev = r.get('gene_evidence') or []
        core = [e for e in ev if e.get('gene', '').startswith('MAT1-') and not e.get('lost_to')]
        for e in ev:
            if e.get('gene') not in ('GAP1', 'BGL1'): continue
            dist = min((max(0, max(c['start'], e['start']) - min(c['end'], e['end'])) for c in core), default='')
            cb = max(core, key=lambda c: c.get('identity') or 0) if core else {}
            print('\t'.join(str(x) for x in (d.name, r['contig'], r['idiomorph'], r['confidence'], e['gene'], e.get('status'), e.get('identity'), e.get('coverage'), dist, cb.get('identity', ''), cb.get('coverage', ''))))
