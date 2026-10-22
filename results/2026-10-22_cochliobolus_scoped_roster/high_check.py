"""For given genomes, print each detected MAT call's genes with identity/coverage/status (proto2)."""
import sys, yaml
for g in sys.argv[2:]:
    d = yaml.load(open(f'{sys.argv[1]}/{g}/detection_report.yaml'), Loader=yaml.CSafeLoader)
    for r in d.get('detected') or []:
        ev = ', '.join(f"{e['gene']} {e.get('identity')}/{e.get('coverage')}/{e.get('status')}{'/lost_to='+str(e['lost_to']) if e.get('lost_to') else ''}" for e in r.get('gene_evidence') or [])
        print(g, r['idiomorph'], r['confidence'], r['locus_class'], '|', ev)
