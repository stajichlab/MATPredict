#!/usr/bin/env python3
"""Summarise v0.6.0 detection reports for the four Tremellomycetes orders.

usage: 01_report_summary.py REPORT_DIR GENOMES_TSV LOCI_TSV QUALITY_TSV OUT_TSV
REPORT_DIR holds <chunk>/runs/<genome>/detection_report.yaml (from reports_all.tar.zst).
"""
import sys, csv, glob, os, collections, yaml
rep, gtsv, ltsv, qtsv, out = sys.argv[1:6]
ORD = ('Tremellales', 'Trichosporonales', 'Filobasidiales', 'Cystofilobasidiales')
G = [r for r in csv.DictReader(open(gtsv), delimiter='\t') if r['order'] in ORD]
Q = {r['genome']: r for r in csv.DictReader(open(qtsv), delimiter='\t')}
loci = collections.defaultdict(list)
for r in csv.DictReader(open(ltsv), delimiter='\t'):
    loci[r['genome']].append(r)
paths = {os.path.basename(os.path.dirname(p)): p for p in glob.glob(rep + '/*/runs/*/detection_report.yaml')}
cols = ['genome', 'order', 'family', 'species', 'status', 'busco', 'n50', 'contigs', 'good_asm',
        'called_families', 'call_classes', 'call_genes', 'verification',
        'n_suppressed', 'sup_best_genes', 'sup_best_family', 'sup_best_reason', 'sup_reasons',
        'nd_reason', 'nd_found', 'nd_missing', 'nd_notsearchable', 'nd_best_fraction']
with open(out, 'w') as o:
    w = csv.writer(o, delimiter='\t'); w.writerow(cols)
    for g in G:
        n = g['genome']; q = Q.get(n, {})
        d = yaml.safe_load(open(paths[n])) if n in paths else {}
        busco = q.get('busco', ''); n50 = int(q['n50']) if q.get('n50') else ''; ctg = int(q['contigs']) if q.get('contigs') else ''
        good = (busco != '' and float(busco) >= 70 and n50 != '' and n50 >= 20000 and ctg != '' and ctg <= 5000)
        L = loci.get(n, [])
        sup = d.get('suppressed_loci') or []
        best = max(sup, key=lambda s: len(s.get('genes_found') or []), default={})
        nd = (d.get('not_detected') or [])
        ndb = max(nd, key=lambda s: s.get('best_fraction_found') or 0, default={})
        w.writerow([n, g['order'], g['family'], g['species'], g['status'], busco, n50, ctg, int(good),
                    '|'.join(sorted({r['family_called'].split(':')[-1] for r in L})),
                    '|'.join(sorted({r['locus_class'] for r in L})),
                    ';'.join(r['genes_found'] for r in L),
                    '|'.join(sorted({('unverified' if 'unverified' in r['verification'] else 'ok') for r in L})),
                    len(sup), '|'.join(best.get('genes_found') or []), (best.get('family') or '').split(':')[-1],
                    best.get('withheld_reason', ''),
                    ';'.join('%s=%d' % kv for kv in collections.Counter(s.get('withheld_reason') for s in sup).items()),
                    (ndb.get('reason') or '')[:120], '|'.join(ndb.get('genes_found') or []),
                    '|'.join(ndb.get('genes_missing') or []), '|'.join(ndb.get('genes_not_searchable') or []),
                    ndb.get('best_fraction_found', '')])
