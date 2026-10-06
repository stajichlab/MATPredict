#!/usr/bin/env python3
"""For failures on good assemblies, what do the withheld clusters in suppressed_loci contain?
usage: 09_withheld_pairs.py REPORT_DIR  (reads failure_classes.tsv; writes withheld_pairs.tsv)
Counts per genome: any withheld cluster with an HD pair (HD1+HD2 or bE+bW), with receptor+precursor,
with >=2 distinct genes of one family, and the largest polished_genes seen in any withheld cluster."""
import csv, glob, os, sys, yaml, collections
rep = sys.argv[1]
paths = {os.path.basename(os.path.dirname(p)): p for p in glob.glob(rep + '/*/runs/*/detection_report.yaml')}
rows = [r for r in csv.DictReader(open('failure_classes.tsv'), delimiter='\t') if r['good_asm'] == '1']
REC = {'pheromone_receptor', 'STE3', 'bar3', 'bbr2', 'pra1'}
out = []
for r in rows:
    d = yaml.safe_load(open(paths[r['genome']])) or {}
    sup = d.get('suppressed_loci') or []
    hdpair = recpre = 0
    maxpol = 0
    for s in sup:
        g = set(s.get('genes_found') or [])
        maxpol = max(maxpol, s.get('polished_genes') or 0)
        if ({'HD1', 'HD2'} <= g) or ({'bE', 'bW'} <= g) or ({'SXI1', 'SXI2'} <= g):
            hdpair += 1
        if (g & REC) and any(x.startswith(('pheromone_B', 'fungal_mating', 'MF', 'mfa', 'bap', 'bbp', 'caax')) for x in g):
            recpre += 1
    out.append([r['genome'], r['order'], r['cls'], r['fail_class'], len(sup), hdpair, recpre, maxpol])
with open('withheld_pairs.tsv', 'w') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(['genome', 'order', 'cls', 'fail_class', 'n_suppressed', 'withheld_hd_pair_clusters', 'withheld_receptor_precursor_clusters', 'max_polished_in_withheld'])
    w.writerows(out)
agg = collections.defaultdict(lambda: [0, 0, 0, 0])
for g, o, c, fc, n, h, p, m in out:
    a = agg[(o, c)]
    a[0] += 1; a[1] += h > 0; a[2] += p > 0; a[3] += (h > 0 or p > 0)
print('order\tcls\tgood_n\thd_pair_withheld\treceptor_precursor_withheld\teither')
for k, v in sorted(agg.items()):
    print(k[0], k[1], *v, sep='\t')
