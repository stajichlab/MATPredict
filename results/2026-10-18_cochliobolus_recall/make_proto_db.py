"""THROWAWAY prototype of remedy (d'): a curated roster that includes the genes that sit next to the MAT idiomorph in Bipolaris.
Copies nothing in the repo; edits a copy of db/ (path argv[1]) in place:
  * db/Ascomycota/order.yml: two new flanking_conserved roster genes, GAP1 and BGL1 (B. maydis proteins AAB84002.1, AAB84005.1);
  * the two Bipolaris maydis records (5016_c4_MAT_MAT1-2, 5016_c5_MAT_MAT1-1) get GAP1 and BGL1 in proteins.faa (gene_index 1, 2).
usage: make_proto_db.py DB_COPY ADJ_REFS.faa"""
import re, sys
from pathlib import Path
db = Path(sys.argv[1])
adj = {}
k = None
for l in open(sys.argv[2]):
    if l.startswith('>'):
        k = l[1:].split('|')[2].strip(); adj[k] = []
    else:
        adj[k].append(l.strip())
adj = {k: ''.join(v) for k, v in adj.items()}
oy = db / 'Ascomycota/order.yml'
t = oy.read_text()
anchor = "      - {name: COX13, role: flanking_conserved}\n"
assert anchor in t
t = t.replace(anchor, anchor + "      - {name: GAP1, role: flanking_conserved}\n      - {name: BGL1, role: flanking_conserved}\n", 1)
oy.write_text(t)
for rec in ('5016_c4_MAT_MAT1-2', '5016_c5_MAT_MAT1-1'):
    f = db / 'Ascomycota/Pleosporales' / rec / 'proteins.faa'
    s = f.read_text().rstrip('\n') + '\n'
    for i, g in ((1, 'GAP1'), (2, 'BGL1')):
        s += f'>{rec}|gene_index={i}|name={g}|role=flanking_conserved\n{adj[g]}\n'
    f.write_text(s)
print('prototype db written', db)
