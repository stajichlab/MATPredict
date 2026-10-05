"""Per-gene and per-gap RNA-seq depth over the NC1011 SLA2-APN2 block (SRR8861595, hisat2, unstranded)."""
import csv, statistics as st, sys
D = sys.argv[1]
dep = {}
for l in open(f'{D}/work/rna/SRR8861595_region_depth.tsv'):
    c, p, d = l.split(); dep[int(p)] = int(d)
def stats(a, b):
    v = [dep.get(i, 0) for i in range(a, b + 1)]
    return (b - a + 1, round(sum(v) / len(v), 1), st.median(v), min(v), max(v), round(sum(x >= 5 for x in v) / len(v), 2))
rows = []
for r in csv.DictReader(open(f'{D}/work/genes.tsv'), delimiter='\t'):
    a, b = int(r['start']), int(r['end'])
    rows.append(('gene', r['protein'], r['product'][:40], a, b) + stats(a, b))
gaps = (('gap', 'prev gene - SLA2', 243655, 243983), ('gap', 'SLA2 - COX13', 247734, 248208),
        ('gap', 'COX13 - APN2', 249502, 249619), ('gap', 'APN2 - H604', 252217, 252999))
for k, n, a, b in gaps: rows.append((k, n, '', a, b) + stats(a, b))
with open(f'{D}/work/rna/coverage_summary.tsv', 'w') as o:
    o.write('type\tname\tproduct\tstart\tend\tlength\tmean_depth\tmedian_depth\tmin\tmax\tfrac_ge5x\n')
    for r in rows: o.write('\t'.join(str(x) for x in r) + '\n')
print('rows', len(rows))
