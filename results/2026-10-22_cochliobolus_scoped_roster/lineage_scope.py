"""For each genome taxid in all.tsv: is Pleosporaceae (28556) in its NCBI lineage (MATPredict's own resolver)?"""
import sys
from MATPredict.db.taxonomy import default_lineage_taxids
out = open(sys.argv[2], 'w')
out.write('genome\ttaxid\tin_pleosporaceae\tlineage_ok\n')
cache = {}
for line in open(sys.argv[1]):
    g, t = line.split()
    t = int(t)
    if t not in cache:
        try:
            lin = set(default_lineage_taxids(t)) | {t}; ok = True
        except Exception as e:
            lin = set(); ok = False
        cache[t] = (28556 in lin, ok)
    a, ok = cache[t]
    out.write(f'{g}\t{t}\t{a}\t{ok}\n')
