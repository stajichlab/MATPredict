#!/usr/bin/env python3
"""Classify failures (uncalled, or called with PR only) on joined.tsv.

usage: 06_classify.py [LINK_BP]   (run in the results directory)

Classes (scan evidence; only good assemblies are classified by gene arrangement):
  d_assembly       BUSCO < 70 or N50 < 20 kb or contigs > 5000
  e_absent         no STE3-like locus and no homeodomain locus found by the scan
  c_linked         HD1-class gene (Pfam PF05920 Homeobox_KN, or miniprot >= 45% identity to a curated HD) and STE3 on one contig
                   within LINK bp: a linked or fused HD+PR arrangement; no single database family
                   expresses it as one locus
  c_split          HD and STE3 both present but on different contigs or farther apart than LINK
  b_receptor_only  STE3 present, no HD locus found
  a_hd_only        HD present, no STE3 found
"""
import csv, collections, sys
LINK = int(sys.argv[1]) if len(sys.argv) > 1 else 150000
R = list(csv.DictReader(open('joined.tsv'), delimiter='\t'))


def num(x):
    return int(x) if x not in ('', None) else None


KEEP = ['genome', 'order', 'species', 'cls', 'good_asm', 'fail_class', 'd_kn_ste3', 'n_ste3', 'n_ste3_caax',
        'n_kn', 'n_hd_mp_strong', 'busco', 'n50', 'contigs']
out = []
tab = collections.Counter()
for r in R:
    if r['cls'] == 'called_other':
        continue
    if r['good_asm'] != '1':
        c = 'd_assembly'
    else:
        ste = num(r['n_ste3']) or 0
        hd = (num(r['n_kn']) or 0) + (num(r['n_hd_mp_strong']) or 0)
        ds = [x for x in (num(r['d_kn_ste3']), num(r['d_hdstrong_ste3'])) if x is not None]
        d = min(ds) if ds else None
        if ste == 0 and hd == 0:
            c = 'e_absent'
        elif ste and hd and d is not None and d <= LINK:
            c = 'c_linked'
        elif ste and hd:
            c = 'c_split'
        elif ste:
            c = 'b_receptor_only'
        else:
            c = 'a_hd_only'
    r['fail_class'] = c
    out.append(r)
    tab[(r['order'], r['cls'], c)] += 1
with open('failure_classes.tsv', 'w') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(KEEP)
    for r in out:
        w.writerow([r[k] for k in KEEP])
cl = ['c_linked', 'c_split', 'b_receptor_only', 'a_hd_only', 'e_absent', 'd_assembly']
print('LINK', LINK)
print('order\tcls\t' + '\t'.join(cl) + '\ttotal')
for o in ['Tremellales', 'Trichosporonales', 'Filobasidiales', 'Cystofilobasidiales']:
    for k in ['uncalled', 'PR_only']:
        v = [tab[(o, k, c)] for c in cl]
        if sum(v):
            print(o, k, *v, sum(v), sep='\t')
