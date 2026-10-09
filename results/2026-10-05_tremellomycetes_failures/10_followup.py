#!/usr/bin/env python3
"""Review follow-ups, from existing tables (joined.tsv).
1. withheld-cluster medians by order and group (all genomes, and the uncalled-good subset).
2. HD allele-type cross-tab: HD1-type (KN or strong miniprot hit) versus HD2-type (homeobox PF00046
   within 150 kb of STE3, no KN) versus neither, by call status, within Trichosporonales and by genus.
Usage: 10_followup.py > followup_summary.txt
"""
import csv, collections, statistics as st
R = list(csv.DictReader(open('joined.tsv'), delimiter='\t'))


def num(x):
    return int(x) if x not in ('', None) else None


print('== withheld clusters (n_suppressed), median [n]')
for o in ['Tremellales', 'Trichosporonales', 'Filobasidiales', 'Cystofilobasidiales']:
    allg = [int(r['n_suppressed']) for r in R if r['order'] == o]
    print(o, 'all genomes', st.median(allg), [len(allg)])
    for c in ['called_other', 'PR_only', 'uncalled']:
        for good in ['1', '0']:
            x = [int(r['n_suppressed']) for r in R if r['order'] == o and r['cls'] == c and r['good_asm'] == good]
            if x:
                print('  ', c, 'good' if good == '1' else 'poor', st.median(x), [len(x)])

print('\n== CAAX near STE3: share of genomes with a strict-CAAX ORF within 10 kb of a STE3 locus')
for o in ['Trichosporonales', 'Filobasidiales', 'Cystofilobasidiales']:
    for c in ['called_other', 'PR_only', 'uncalled']:
        x = [r for r in R if r['order'] == o and r['cls'] == c and r['good_asm'] == '1']
        if x:
            print(o, c, sum((num(r['n_ste3_caax']) or 0) > 0 for r in x), '/', len(x))


def allele(r, link=150000):
    kn = (num(r['n_kn']) or 0) + (num(r['n_hd_mp_strong']) or 0) > 0
    dhb = num(r['d_hb_ste3'])
    if kn:
        return 'HD1type(KN)'
    if dhb is not None and dhb <= link:
        return 'HD2type(homeobox near STE3)'
    return 'neither'


print('\n== Trichosporonales, good assemblies: allele type x call status')
T = [r for r in R if r['order'] == 'Trichosporonales' and r['good_asm'] == '1']
ct = collections.Counter((r['cls'], allele(r)) for r in T)
for c in ['called_other', 'PR_only', 'uncalled']:
    print(c, {a: ct[(c, a)] for a in ['HD1type(KN)', 'HD2type(homeobox near STE3)', 'neither']})
print('\nallele x CAAX-near-STE3 x call status (KN, CAAX):')
cc = collections.Counter((r['cls'], (num(r['n_kn']) or 0) > 0, (num(r['n_ste3_caax']) or 0) > 0) for r in T)
for c in ['called_other', 'PR_only', 'uncalled']:
    print(c, 'KN+/CAAX-', cc[(c, True, False)], 'KN+/CAAX+', cc[(c, True, True)], 'KN-/CAAX+', cc[(c, False, True)], 'KN-/CAAX-', cc[(c, False, False)])

print('\n== by genus (Trichosporonales, good): genus, status -> KN+ / n ; allele classes')
by = collections.defaultdict(list)
for r in T:
    by[r['species'].split()[0]].append(r)
print('genus\tstatus\tn\tKN+\tHD2type\tneither\tCAAX+')
for g, x in sorted(by.items(), key=lambda kv: -len(kv[1])):
    for c in ['called_other', 'PR_only', 'uncalled']:
        y = [r for r in x if r['cls'] == c]
        if y:
            a = collections.Counter(allele(r) for r in y)
            print(g, c, len(y), a['HD1type(KN)'], a['HD2type(homeobox near STE3)'], a['neither'], sum((num(r['n_ste3_caax']) or 0) > 0 for r in y), sep='\t')

print('\n== species with both a called/PR-only and an uncalled genome (Trichosporonales, good): allele types')
sp = collections.defaultdict(list)
for r in T:
    sp[r['species']].append(r)
n = 0
for s, x in sorted(sp.items()):
    cl = {r['cls'] for r in x}
    if 'uncalled' in cl and len(cl) > 1:
        n += 1
        print(s, [(r['cls'], allele(r).split('(')[0]) for r in x])
print('species with both:', n)
