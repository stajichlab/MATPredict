"""Collect one row per genome per variant from detection_report.yaml.
usage: collect.py RUNS_DIR VARIANT [VARIANT ...] > out.tsv
Columns: variant genome status report_sha (sha256 of the report minus its `run` block) calls (';'-joined locus keys)
A locus key = family|contig|start|end|idiomorph|confidence|locus_class|genes_found(comma)."""
import hashlib, json, sys
from pathlib import Path
import yaml
try:
    Loader = yaml.CSafeLoader
except AttributeError:
    Loader = yaml.SafeLoader
runs = Path(sys.argv[1])
print('variant\tgenome\tstatus\treport_sha\tcalls_sha\tn_calls\tcalls\twall')
for v in sys.argv[2:]:
    for d in sorted((runs / v).iterdir()):
        rep = d / 'detection_report.yaml'
        wall = (d / 'wall_seconds').read_text().strip() if (d / 'wall_seconds').exists() else ''
        if not rep.exists() or rep.stat().st_size == 0:
            print(f'{v}\t{d.name}\tno_report\t\t\t0\t\t{wall}'); continue
        doc = yaml.load(rep.read_text(), Loader=Loader)
        doc.pop('run', None)
        sha = hashlib.sha256(json.dumps(doc, sort_keys=True, default=str).encode()).hexdigest()[:16]
        calls = []
        for r in doc.get('detected') or []:
            calls.append('|'.join(str(x) for x in (r.get('family'), r.get('contig'), r.get('start'), r.get('end'), r.get('idiomorph'),
                                                   r.get('confidence'), r.get('locus_class'), ','.join(r.get('genes_found') or []))))
        calls.sort()
        csha = hashlib.sha256(';'.join(calls).encode()).hexdigest()[:16]
        print(f"{v}\t{d.name}\t{'called' if calls else 'not_detected'}\t{sha}\t{csha}\t{len(calls)}\t{';'.join(calls)}\t{wall}")
