import json,sys
from collections import Counter
for f in sys.argv[1:]:
    rows=[json.loads(l) for l in open(f)]
    ev=[r for r in rows if r['kind']=='evidence']
    print(f.split('/')[0], dict(Counter(r['kind'] for r in rows)), 'clusters',len(ev),'admitted',sum(r['admitted'] for r in ev),
          'core>=80:',sum(1 for r in ev if 'core_MAT' in r['roles'] and (r['best_identity'] or 0)>=80))
    for r in sorted(ev,key=lambda r:-(r['best_identity'] or 0))[:4]:
        print('   ',r['contig'],r['cluster_start'],r['cluster_end'],'genes',r['gene_count'],r['roles'],r['best_identity'],'admitted',r['admitted'])
