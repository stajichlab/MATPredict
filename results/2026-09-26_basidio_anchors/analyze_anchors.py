"""Adjacency of HD-class hits to MIP1 / beta-fg, and receptor hits to ICMT,
per pilot genome. Anchor = top-bitscore tblastn HSP of that query.
Usage: python analyze_anchors.py [window_bp]"""
import sys, collections, os
W = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
D = os.path.dirname(os.path.abspath(__file__))
pilot = [l.rstrip('\n').split('\t') for l in open(f'{D}/pilot.tsv')]
sub = {'Agaricomycetes':'Agaricomycotina','Dacrymycetes':'Agaricomycotina','Tremellomycetes':'Agaricomycotina',
       'Ustilaginomycetes':'Ustilaginomycotina','Malasseziomycetes':'Ustilaginomycotina','Exobasidiomycetes':'Ustilaginomycotina',
       'Microbotryomycetes':'Pucciniomycotina','Pucciniomycetes':'Pucciniomycotina','Mixiomycetes':'Pucciniomycotina',
       'Wallemiomycetes':'Wallemiomycotina'}
def load(g):
    hs = []
    for l in open(f'{D}/hits/{g}.tsv'):
        q, s, pid, ln, a, b, e, bits = l.split('\t')
        if q == 'none': continue
        hs.append((q.split('|')[0], q, s, float(pid), int(ln), min(int(a), int(b)), max(int(a), int(b)), float(e), float(bits)))
    return hs
def top(hs, cls):
    c = [h for h in hs if h[0] == cls]
    return max(c, key=lambda h: h[8]) if c else None
def near(anchor, hs, cls, w):
    if not anchor: return []
    return [h for h in hs if h[0] == cls and h[2] == anchor[2] and h[5] - w <= anchor[6] and h[6] + w >= anchor[5]]
def clusters(hs, cls, gap=10000):
    c = sorted([h for h in hs if h[0] == cls], key=lambda h: (h[2], h[5]))
    out = []
    for h in c:
        if out and out[-1][0] == h[2] and h[5] <= out[-1][2] + gap:
            out[-1][2] = max(out[-1][2], h[6]); out[-1][3].append(h)
        else: out.append([h[2], h[5], h[6], [h]])
    return out
rows = []
print(f'window {W} bp; anchor = top tblastn HSP; HD/receptor hit = any HSP e<=1e-5')
print('subphylum\torder\tspecies\tMIP(pid,contig)\tHD_near_MIP(best pid)\tBFG(pid)\tHD_near_BFG\tMIP-BFG_same_contig_dist\tHD_clusters\tbg_HD_frac\tICMT_pid\tPR_near_ICMT\tPR_clusters')
for asm, taxid, cls, order, fam, sp in pilot:
    hs = load(asm)
    glen = sum(int(l.split('\t')[1]) for l in open(f'{D}/hits/{asm}.len'))
    mip, bfg, icmt = top(hs, 'MIP1'), top(hs, 'BFG'), top(hs, 'ICMT')
    hm, hb = near(mip, hs, 'HD', W), near(bfg, hs, 'HD', W)
    pr = near(icmt, hs, 'PR', W)
    hdc, prc = clusters(hs, 'HD'), clusters(hs, 'PR')
    # background: fraction of genome within W of any HD cluster (random-point chance of adjacency)
    bg = min(1.0, sum((c[2] - c[1]) + 2 * W for c in hdc) / glen) if glen else 0
    d = ''
    if mip and bfg and mip[2] == bfg[2]:
        d = str(max(0, max(mip[5], bfg[5]) - min(mip[6], bfg[6])))
    fmt = lambda h: f'{h[3]:.0f}' if h else '-'
    rows.append((sub.get(cls, '?'), order, sp, bool(hm), bool(hb), bool(pr), bg))
    print('\t'.join([sub.get(cls, '?'), order, sp, f'{fmt(mip)}', f'{len(hm)}({max([h[3] for h in hm]) if hm else 0:.0f})',
                     fmt(bfg), str(len(hb)), d or '-', str(len(hdc)), f'{bg:.3f}', fmt(icmt), str(len(pr)), str(len(prc))]))
print()
agg = collections.defaultdict(lambda: [0, 0, 0, 0, 0.0])
for s, o, sp, a, b, p, bg in rows:
    x = agg[s]; x[0] += 1; x[1] += a; x[2] += b; x[3] += p; x[4] += bg
print('subphylum\tgenomes\tHD_near_MIP\tHD_near_BFG\tPR_near_ICMT\tmean_bg_HD_frac')
for s, x in agg.items():
    print(f'{s}\t{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]/x[0]:.3f}')
