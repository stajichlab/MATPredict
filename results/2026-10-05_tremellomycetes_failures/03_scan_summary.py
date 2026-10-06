#!/usr/bin/env python3
"""Per-genome summary of the pipeline-independent scan (scan_ste3_hd.py).

usage: 03_scan_summary.py SCAN_DIR IDS_TXT OUT_TSV
SCAN_DIR holds <genome>.genes.tsv and <genome>.contigs.tsv.

Per genome: number of STE3-like loci (Pfam PF02076 or miniprot of 1,095 STE3
queries), how many have a strict-CAAX ORF within 10 kb, number of homeodomain
loci (Pfam PF05920 = Homeobox_KN, the HD1 class; PF00046 = Homeobox), and the
smallest same-contig distance between an HD locus and an STE3 locus.
"""
import csv
import os
import sys
from collections import defaultdict

scan, ids, out = sys.argv[1:4]
GAP = 3000
CAAX_W = 10000


def merge(rows, gap):
    """rows: (contig, start, end, tag); returns merged loci dicts per contig."""
    by = defaultdict(list)
    for r in rows:
        by[r[0]].append(r)
    loci = []
    for ctg, lst in by.items():
        lst.sort(key=lambda r: r[1])
        cur = None
        for c, s, e, tag in lst:
            if cur and s <= cur['end'] + gap:
                cur['end'] = max(cur['end'], e)
                cur['tags'].add(tag)
            else:
                if cur:
                    loci.append(cur)
                cur = {'contig': c, 'start': s, 'end': e, 'tags': {tag}}
        if cur:
            loci.append(cur)
    return loci


def dist(a, b):
    if a['contig'] != b['contig']:
        return None
    if a['start'] <= b['end'] and b['start'] <= a['end']:
        return 0
    return min(abs(a['start'] - b['end']), abs(b['start'] - a['end']))


def mind(A, B):
    best = None
    for a in A:
        for b in B:
            d = dist(a, b)
            if d is not None and (best is None or d < best):
                best = d
    return best


cols = ['genome', 'n_contigs_scanned', 'n_ste3', 'n_ste3_hmm', 'n_ste3_caax', 'n_hd_loci', 'n_kn', 'n_hdbox',
        'n_hd_mp_strong', 'n_kn_hd_pairs_5kb', 'd_kn_ste3', 'd_kn_ste3caax', 'd_hdany_ste3', 'd_hdany_ste3caax',
        'd_hdstrong_ste3', 'd_hb_ste3', 'status']
with open(out, 'w') as fo:
    w = csv.writer(fo, delimiter='\t')
    w.writerow(cols)
    for g in open(ids).read().split():
        gp = os.path.join(scan, g + '.genes.tsv')
        if not os.path.exists(gp):
            w.writerow([g] + [''] * (len(cols) - 2) + ['missing'])
            continue
        clen = {}
        for line in open(os.path.join(scan, g + '.contigs.tsv')):
            k, v = line.split()
            clen[k] = int(v)
        ste3, hd, caax, hdmp = [], [], [], []
        for r in csv.DictReader(open(gp), delimiter='\t'):
            s, e = int(r['start']), int(r['end'])
            src = r['source']
            if r['kind'] == 'STE3':
                ste3.append((r['contig'], s, e, 'hmm' if src.startswith('hmm') else 'mp'))
            elif r['kind'] == 'HD':
                if src.startswith('hmm'):
                    hd.append((r['contig'], s, e, 'KN' if src.endswith('PF05920') else 'HB'))
                else:
                    name, ident, cov = r['detail'].rsplit(':', 2)
                    # strong = identity 0.45+ to a curated HD protein over 0.5+ of the query
                    if float(ident) >= 0.45 and float(cov) >= 0.5:
                        hdmp.append((r['contig'], s, e, 'mpstrong'))
                    if float(ident) >= 0.35:
                        hd.append((r['contig'], s, e, 'mp'))
            elif r['kind'] == 'CAAX':
                caax.append((r['contig'], s, e, 'caax'))
        S = merge(ste3, GAP)
        H = merge(hd, GAP)
        K = [h for h in H if 'KN' in h['tags']]
        B = [h for h in H if 'HB' in h['tags']]
        C = merge(caax, 1)
        cx = defaultdict(list)
        for c in C:
            cx[c['contig']].append(c)
        Sc = []
        for s in S:
            near = [c for c in cx.get(s['contig'], []) if c['start'] <= s['end'] + CAAX_W and c['end'] >= s['start'] - CAAX_W]
            if near:
                Sc.append(s)
        pairs = 0
        for k in K:
            for b in B:
                d = dist(k, b)
                if d is not None and d <= 5000:
                    pairs += 1
        HS = merge(hdmp, GAP)
        w.writerow([g, len(clen), len(S), sum('hmm' in s['tags'] for s in S), len(Sc), len(H), len(K), len(B),
                    len(HS), pairs, mind(K, S), mind(K, Sc), mind(H, S), mind(H, Sc), mind(HS, S), mind(B, S), 'ok'])
