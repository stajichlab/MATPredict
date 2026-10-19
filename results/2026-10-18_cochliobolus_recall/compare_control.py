"""Baseline (audit run, current main) vs prototype roster on the control genomes (called, outside the clade).
usage: compare_control.py audit_calls.tsv proto_control.tsv out.tsv"""
import csv, sys
from collections import defaultdict, Counter
def load(p):
    d = defaultdict(list)
    for r in csv.DictReader(open(p), delimiter='\t'):
        if r['n_detected'] != '0':
            d[r['genome']].append(r)
        else:
            d.setdefault(r['genome'], [])
    return d
base = load(sys.argv[1]); prot = load(sys.argv[2])
rows = []; res = Counter()
for g in prot:
    b = base.get(g, []); p = prot[g]
    bi = sorted(r['idiomorph'] for r in b); pi = sorted(r['idiomorph'] for r in p)
    bc = ','.join(sorted(r['confidence'] for r in b)); pc = ','.join(sorted(r['confidence'] for r in p))
    if not p and b: o = 'LOST_CALL'
    elif bi != pi: o = 'IDIOMORPH_SET_CHANGED'
    elif bc != pc: o = 'confidence_changed'
    else: o = 'same'
    # same locus? compare coordinates of first call
    res[o] += 1
    rows.append(dict(genome=g, base_idiomorph=','.join(bi), proto_idiomorph=','.join(pi), base_conf=bc, proto_conf=pc, base_genes=';'.join(r['genes_found'] for r in b),
                     proto_genes=';'.join(r['genes_found'] for r in p), outcome=o))
print(res, 'n', len(rows))
with open(sys.argv[3], 'w') as f:
    w = csv.DictWriter(f, list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
