"""One row per detected locus (or one 'none' row) for every run dir under argv[1]. Output TSV on stdout.
columns: genome, n_detected, contig, start, end, idiomorph, confidence, locus_class, genes_found, polished_genes, flags (strong_single_gene / split_locus ...), core_identity"""
import sys, os, yaml
root = sys.argv[1]
print('\t'.join(['genome', 'n_detected', 'contig', 'start', 'end', 'idiomorph', 'confidence', 'locus_class', 'genes_found', 'polished_genes', 'ssg_identity', 'ssg_coverage', 'two_idiomorphs', 'n_suppressed']))
for g in sorted(os.listdir(root)):
    p = f'{root}/{g}/detection_report.yaml'
    if not os.path.exists(p): continue
    y = yaml.safe_load(open(p))
    det = y.get('detected') or []
    two = ','.join(str(x) for x in (y.get('two_idiomorphs') or []))[:60]
    if not det:
        print('\t'.join([g, '0', '', '', '', '', '', '', '', '', '', '', two, str(y.get('suppressed_unpolished'))]))
    for d in det:
        s = d.get('strong_single_gene') or {}
        print('\t'.join(str(x) for x in [g, len(det), d.get('contig'), d.get('start'), d.get('end'), d.get('idiomorph'), d.get('confidence'), d.get('locus_class'),
                                         ','.join(d.get('genes_found') or []), d.get('polished_genes'), s.get('core_identity', ''), s.get('core_coverage', ''), two, y.get('suppressed_unpolished')]))
