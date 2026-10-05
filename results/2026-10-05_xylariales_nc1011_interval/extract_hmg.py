"""Pick MATA_HMG-family domain regions (HMG_box PF00505 / HMG_box_2 PF09011) from proteomes.

Reads: proteins.faa, hmmsearch domtblout. Writes hmg_domains.faa: best HMG-box domain per protein,
env coordinates extended by 35 aa downstream (the paper's 'core' = HMG box + ~35 aa), min 60 aa.
MATalpha_HMGbox (PF04769) hits are listed separately (different family; not aligned here).
"""
import sys
faa, dom, out, listing = sys.argv[1:5]
seqs = {}; k = None
for l in open(faa):
    if l[0] == '>': k = l[1:].strip(); seqs[k.split()[0]] = [k, []]
    else: seqs[k.split()[0]][1].append(l.strip())
best = {}; alpha = set()
for l in open(dom):
    if l[0] == '#': continue
    f = l.split()
    q, acc, score, ienv, jenv = f[0], f[4], float(f[13]), int(f[19]), int(f[20])
    if acc.startswith('PF04769'): alpha.add(q); continue
    if acc.startswith(('PF00505', 'PF09011')):
        if q not in best or score > best[q][0]: best[q] = (score, ienv, jenv, acc)
n = 0
with open(out, 'w') as o, open(listing, 'w') as li:
    li.write('protein\thmm\tscore\tenv_start\tenv_end\theader\n')
    for q, (s, a, b, acc) in best.items():
        full = ''.join(seqs[q][1]); sub = full[a - 1:min(len(full), b + 35)]
        if len(sub) < 60: continue
        o.write('>%s\n%s\n' % (q, sub)); n += 1
        li.write('%s\t%s\t%.1f\t%d\t%d\t%s\n' % (q, acc, s, a, b, seqs[q][0][:120]))
print('hmg domains:', n, ' alpha-box proteins:', len(alpha))
