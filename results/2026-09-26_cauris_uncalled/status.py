import yaml,os,csv,sys
scan=sys.argv[1]
rows=list(csv.reader(open('2026-09-26_cauris_uncalled/cauris_samples.csv')))
for r in rows:
    a=r[0]; p=f'{scan}/runs/{a}/detection_report.yaml'
    if not os.path.exists(p): print(a,'NO_REPORT',r[2]); continue
    y=yaml.safe_load(open(p))
    det=y.get('detected') or []; sup=y.get('suppressed_loci') or []
    if not det:
        print(a, r[2], y.get('routing_mode'), 'UNCALLED', [(s['contig'],s['start'],s['end'],s['idiomorph'],s['polished_genes'],s['genes_found']) for s in sup])
