#!/usr/bin/env python3
"""Build hd_queries.faa: curated Basidiomycota HD-class proteins from db/ plus the
Phaffia rhodozyma HD1/HD2 deposits (KU315762, KU315770, KU315771, KU315777, KU315779).

usage: make_hd_queries.py DB_BASIDIO_DIR PHAFFIA_FAA OUT_FAA
"""
import glob
import re
import sys

db, phaffia, outp = sys.argv[1:4]
names = {'HD1', 'HD2', 'SXI1', 'SXI2', 'bE', 'bW', 'Y', 'Z'}
out = []
seen = set()
for f in sorted(glob.glob(db + '/*/*/proteins.faa')):
    rid = f.split('/')[-2]
    name = None
    seq = []

    def flush():
        if name in names and seq:
            s = ''.join(seq)
            if (name, s) not in seen:
                seen.add((name, s))
                out.append(f">{rid}|{name}\n{s}\n")
    for line in open(f):
        if line.startswith('>'):
            flush()
            seq = []
            m = re.search(r'name=([^|]+)', line)
            name = m.group(1) if m else None
        else:
            seq.append(line.strip())
    flush()
name = None
seq = []
recs = []
for line in open(phaffia):
    if line.startswith('>'):
        if name:
            recs.append((name, ''.join(seq)))
        m = re.search(r'gene=(HD\d).*protein_id=(\S+?)\]', line)
        name = f"PHAFFIA_{m.group(2)}|{m.group(1)}"
        seq = []
    else:
        seq.append(line.strip())
recs.append((name, ''.join(seq)))
for n, s in recs:
    out.append(f">{n}\n{s}\n")
open(outp, 'w').write(''.join(out))
print(len(out), 'queries')
