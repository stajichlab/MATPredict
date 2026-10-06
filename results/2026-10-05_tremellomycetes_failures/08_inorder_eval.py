#!/usr/bin/env python3
"""Evaluate the in-order reference test. usage: 08_inorder_eval.py INORDER_DIR
A KN hit = miniprot alignment of an in-order KN fragment, identity >= 0.5, query coverage >= 0.5;
STE3 hit likewise (identity >= 0.5). Loci merged at 3 kb. Linked = KN and STE3 on one contig within 150 kb."""
import csv, collections, os, sys
d = sys.argv[1]
C = {r['genome']: r for r in csv.DictReader(open('failure_classes.tsv'), delimiter='\t')}
J = {r['genome']: r for r in csv.DictReader(open('joined.tsv'), delimiter='\t')}


def loci(rows):
    rows.sort(); out = []
    for c, s, e in rows:
        if out and out[-1][0] == c and s <= out[-1][2] + 3000:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return out


res = []
for g, j in J.items():
    p = os.path.join(d, g + '.inorder.tsv')
    if not os.path.exists(p):
        continue
    kn, st = [], []
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if float(f[6]) >= 0.5 and float(f[7]) >= 0.5:
            (kn if f[4] == 'KN' else st).append((f[0], int(f[2]), int(f[3])))
    K, S = loci(kn), loci(st)
    dm = None
    for a in K:
        for b in S:
            if a[0] == b[0]:
                x = 0 if (a[1] <= b[2] and b[1] <= a[2]) else min(abs(a[1] - b[2]), abs(b[1] - a[2]))
                dm = x if dm is None else min(dm, x)
    res.append((g, j['order'], j['cls'], j['good_asm'], len(K), len(S), dm))
with open('inorder_test.tsv', 'w') as fo:
    w = csv.writer(fo, delimiter='\t'); w.writerow(['genome', 'order', 'cls', 'good_asm', 'n_kn', 'n_ste3', 'd_kn_ste3'])
    w.writerows(res)
agg = collections.defaultdict(lambda: [0, 0, 0, 0])
for g, o, cl, good, k, s, dm in res:
    if good != '1':
        continue
    a = agg[(o, cl)]
    a[0] += 1; a[1] += k > 0; a[2] += s > 0; a[3] += dm is not None and dm <= 150000
print('order\tcls\tn\tKN>0\tSTE3>0\tlinked<=150kb')
for k, v in sorted(agg.items()):
    print(k[0], k[1], *v, sep='\t')
