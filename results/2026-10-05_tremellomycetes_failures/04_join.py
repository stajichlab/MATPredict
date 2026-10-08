#!/usr/bin/env python3
"""Join report_summary.tsv and scan_summary.tsv; print per-order tables and write joined.tsv."""
import csv, collections, statistics
R = {r['genome']: r for r in csv.DictReader(open('report_summary.tsv'), delimiter='\t')}
S = {r['genome']: r for r in csv.DictReader(open('scan_summary.tsv'), delimiter='\t')}
rows = []
def num(x):
    return int(x) if x not in ('', None) else None
for g, r in R.items():
    s = S.get(g, {})
    j = dict(r); j.update({k: v for k, v in s.items() if k != 'genome'})
    called = r['status'] == 'called'
    prcalled = r['called_families'] == 'PR'
    j['cls'] = 'uncalled' if not called else ('PR_only' if prcalled else 'called_other')
    rows.append(j)
with open('joined.tsv', 'w') as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
def med(v):
    v = [x for x in v if x is not None]
    return statistics.median(v) if v else None
for o in ['Tremellales', 'Trichosporonales', 'Filobasidiales', 'Cystofilobasidiales']:
    print('==', o)
    for cls in ['called_other', 'PR_only', 'uncalled']:
        for good in ['1', '0']:
            x = [r for r in rows if r['order'] == o and r['cls'] == cls and r['good_asm'] == good and r.get('status') != 'missing']
            if not x: continue
            n = len(x)
            ste = sum(num(r['n_ste3']) > 0 for r in x)
            kn = sum(num(r['n_kn']) > 0 for r in x)
            hdb = sum(num(r['n_hd_mp_strong']) > 0 for r in x)
            both = sum(num(r['n_ste3']) > 0 and (num(r['n_kn']) > 0 or num(r['n_hd_mp_strong']) > 0) for r in x)
            dk = [num(r['d_kn_ste3']) for r in x]
            dh = [num(r['d_hdstrong_ste3']) for r in x]
            same = sum(d is not None for d in dk)
            sames = sum(d is not None for d in dh)
            near = sum(d is not None and d <= 20000 for d in dh)
            near2 = sum(d is not None and d <= 20000 for d in dk)
            print(f'{cls:13s} good={good} n={n:3d} STE3>0:{ste:3d} KN>0:{kn:3d} HDstrong>0:{hdb:3d} STE3&HD:{both:3d} '
                  f'sameContig(KN,STE3):{same:3d} (HDstrong):{sames:3d} <=20kb KN:{near2:3d} HDs:{near:3d} medSTE3={med([num(r["n_ste3"]) for r in x])}')
