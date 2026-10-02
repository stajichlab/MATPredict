"""Adjacency of HD-class hits (from hits/) to candidate non-Agaricomycete
anchors (hits2/: HELI, U403, U400, IES1 from Ustilaginomycotina b-loci;
SF3B5 from Microbotryomycetes). Anchor = top tblastn HSP."""
import sys, os, collections
W = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
D = os.path.dirname(os.path.abspath(__file__))
exec(open(f'{D}/analyze_anchors.py').read().split('rows = []')[0].replace("W = int(sys.argv[1]) if len(sys.argv) > 1 else 20000", ""))
def load2(g):
    hs = []
    for l in open(f'{D}/hits2/{g}.tsv'):
        q, s, pid, ln, a, b, e, bits = l.split('\t')
        if q == 'none': continue
        hs.append((q.split('|')[0], q, s, float(pid), int(ln), min(int(a), int(b)), max(int(a), int(b)), float(e), float(bits)))
    return hs
A = ['HELI', 'U403', 'U400', 'IES1', 'SF3B5']
print('subphylum\torder\tspecies\t' + '\t'.join(f'{a}(pid:nearHD)' for a in A))
agg = collections.defaultdict(lambda: collections.Counter())
for asm, taxid, cls, order, fam, sp in pilot:
    hd = load(asm); h2 = load2(asm)
    cells = []
    s = sub.get(cls, '?'); agg[s]['n'] += 1
    for a in A:
        t = top(h2, a)
        n = near(t, hd, 'HD', W) if t else []
        if n: agg[s][a] += 1
        cells.append(f'{t[3]:.0f}:{len(n)}' if t else '-')
    print('\t'.join([s, order, sp] + cells))
print('\nsubphylum\tgenomes\t' + '\t'.join(A))
for s, c in agg.items(): print('\t'.join([s, str(c['n'])] + [str(c[a]) for a in A]))
