#!/usr/bin/env python3
"""Classify failures (uncalled, or called with PR only) on joined.tsv.

usage: 06_classify.py [LINK_BP]   (run in the results directory; default 150000)
writes failure_classes.tsv (at LINK_BP) and prints the class table.

Classes (good assemblies; poor assemblies are d_assembly):
  d_assembly            BUSCO < 70, N50 < 20 kb or contigs > 5000 (or no BUSCO value)
  e_absent              no STE3-like locus and no HD1-class or homeobox locus found
  linked_hd1            HD1-class gene (Pfam PF05920 KN, or miniprot >= 45% identity to a curated HD)
                        on one contig with a STE3 locus within LINK bp
  linked_homeobox_only  no HD1-class gene within LINK of STE3, but a homeodomain (PF00046) locus is
                        (this is where an HD2-type allele would fall)
  hd1_found_not_linked  HD1-class gene present, but not within LINK of STE3 (other contig or farther)
  no_hd1_domain_found   STE3 present; no HD1-class (KN) gene found anywhere and no homeodomain
                        locus within LINK of STE3 (not a statement that HD is absent)
"""
import csv, collections, sys
LINK = int(sys.argv[1]) if len(sys.argv) > 1 else 150000
R = list(csv.DictReader(open('joined.tsv'), delimiter='\t'))


def num(x):
    return int(x) if x not in ('', None) else None


def classify(r, link):
    if r['good_asm'] != '1':
        return 'd_assembly'
    ste = num(r['n_ste3']) or 0
    hd1 = (num(r['n_kn']) or 0) + (num(r['n_hd_mp_strong']) or 0)
    d1 = [x for x in (num(r['d_kn_ste3']), num(r['d_hdstrong_ste3'])) if x is not None]
    d1 = min(d1) if d1 else None
    dhb = num(r['d_hb_ste3'])
    if ste == 0 and hd1 == 0 and (num(r['n_hdbox']) or 0) == 0:
        return 'e_absent'
    if ste and d1 is not None and d1 <= link:
        return 'linked_hd1'
    if ste and dhb is not None and dhb <= link:
        return 'linked_homeobox_only'
    if ste and hd1:
        return 'hd1_found_not_linked'
    return 'no_hd1_domain_found'


KEEP = ['genome', 'order', 'species', 'cls', 'good_asm', 'fail_class', 'd_kn_ste3', 'd_hdstrong_ste3', 'd_hb_ste3',
        'd_hdany_ste3', 'n_ste3', 'n_ste3_caax', 'n_kn', 'n_hdbox', 'n_hd_mp_strong', 'busco', 'n50', 'contigs']
CL = ['linked_hd1', 'linked_homeobox_only', 'hd1_found_not_linked', 'no_hd1_domain_found', 'e_absent', 'd_assembly']


def table(link):
    tab = collections.Counter()
    out = []
    for r in R:
        if r['cls'] == 'called_other':
            continue
        r = dict(r)
        r['fail_class'] = classify(r, link)
        out.append(r)
        tab[(r['order'], r['cls'], r['fail_class'])] += 1
    return out, tab


def show(tab, link):
    print('LINK', link)
    print('order\tcls\t' + '\t'.join(CL) + '\ttotal')
    tot = [0] * len(CL)
    for o in ['Tremellales', 'Trichosporonales', 'Filobasidiales', 'Cystofilobasidiales']:
        for k in ['uncalled', 'PR_only']:
            v = [tab[(o, k, c)] for c in CL]
            if sum(v):
                print(o, k, *v, sum(v), sep='\t')
                tot = [a + b for a, b in zip(tot, v)]
    print('total', '', *tot, sum(tot), sep='\t')


if __name__ == '__main__':
    out, tab = table(LINK)
    with open('failure_classes.tsv', 'w') as fo:
        w = csv.writer(fo, delimiter='\t')
        w.writerow(KEEP)
        for r in out:
            w.writerow([r[k] for k in KEEP])
    show(tab, LINK)
    print()
    for lk in (50000, 100000, 200000, 300000):
        show(table(lk)[1], lk)
