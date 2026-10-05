"""Campaign overview tables and SVG figures from the committed campaign tables (stdlib only).

Inputs (all in results/): 2026-10-03_ascomycota_v060/{by_class,loci}.tsv, 2026-10-03_basidiomycota_v060/{by_order,loci}.tsv,
2026-10-03_mucoromycotina_mat/calls.tsv.  Outputs next to this script: campaign_summary.tsv, three SVG figures.
Colours: reference palette slots 1-4 (blue, orange, aqua, yellow), validated adjacent-pair CVD and normal-vision floors; the aqua and yellow
contrast warning is covered by direct labels plus the tables in the analysis note.
"""
import csv, collections, os, html
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
OUT = os.path.dirname(os.path.abspath(__file__))
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']

def rd(p): return list(csv.DictReader(open(os.path.join(R, p)), delimiter='\t'))

def esc(s): return html.escape(str(s))

def svg_open(w, h, title, desc):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-labelledby="t d" '
            f'font-family="system-ui, -apple-system, Segoe UI, Helvetica, Arial, sans-serif">\n<title id="t">{esc(title)}</title><desc id="d">{esc(desc)}</desc>\n'
            f'<rect width="{w}" height="{h}" fill="{SURF}"/>\n')

def hbars(path, title, subtitle, rows, legend, xlabel='% of genomes with at least one MAT call', tick_name=None):
    """rows: (label, n, pct, color_index[, tick_pct]); tick_pct draws a marker (e.g. the rate without PR-only calls)"""
    W, L, Rm, top, bh, gap = (1060 if tick_name else 900), 250, (230 if tick_name else 70), 80, 14, 9
    H = top + len(rows) * (bh + gap) + 52
    pw = W - L - Rm
    s = svg_open(W, H, title, subtitle)
    s += f'<text x="20" y="28" font-size="17" font-weight="600" fill="{INK}">{esc(title)}</text>\n'
    s += f'<text x="20" y="48" font-size="12.5" fill="{INK2}">{esc(subtitle)}</text>\n'
    x = L
    for i, name in legend:
        s += f'<rect x="{x}" y="58" width="12" height="12" rx="3" fill="{C[i]}"/><text x="{x+18}" y="68" font-size="12.5" fill="{INK2}">{esc(name)}</text>\n'
        x += 18 + 7 * len(name) + 24
    for k in (0, 25, 50, 75, 100):
        gx = L + pw * k / 100
        s += f'<line x1="{gx:.1f}" y1="{top-8}" x2="{gx:.1f}" y2="{H-40}" stroke="{GRID}" stroke-width="1"/>\n'
        s += f'<text x="{gx:.1f}" y="{H-24}" font-size="11.5" fill="{INK2}" text-anchor="middle">{k}%</text>\n'
    if tick_name:
        s += f'<line x1="{x+6}" y1="57" x2="{x+6}" y2="71" stroke="{INK}" stroke-width="2.5"/><text x="{x+14}" y="68" font-size="12.5" fill="{INK2}">{esc(tick_name)}</text>\n'
    for i, row in enumerate(rows):
        lab, n, pct, ci = row[:4]; tick = row[4] if len(row) > 4 else None
        y = top + i * (bh + gap)
        s += f'<text x="{L-10}" y="{y+bh-3}" font-size="12.5" fill="{INK}" text-anchor="end">{esc(lab)} <tspan fill="{INK2}">n={n:,}</tspan></text>\n'
        s += f'<rect x="{L}" y="{y}" width="{max(pw*pct/100,1):.1f}" height="{bh}" rx="3" fill="{C[ci]}"><title>{esc(lab)}: {pct:.1f}% of {n:,} genomes</title></rect>\n'
        vx = L + pw * pct / 100 + 6
        if tick is not None:
            tx = L + pw * tick / 100
            s += f'<line x1="{tx:.1f}" y1="{y-3}" x2="{tx:.1f}" y2="{y+bh+3}" stroke="{INK}" stroke-width="2.5"><title>{esc(lab)}: {tick:.1f}% without PR-only calls</title></line>\n'
        s += f'<text x="{vx:.1f}" y="{y+bh-3}" font-size="12" fill="{INK}">{pct:.1f}%' + (f' <tspan fill="{INK2}">({tick:.1f}% without PR-only)</tspan>' if tick is not None else '') + '</text>\n'
    s += f'<text x="{L+pw/2}" y="{H-6}" font-size="11.5" fill="{INK2}" text-anchor="middle">{esc(xlabel)}</text>\n</svg>\n'
    open(path, 'w').write(s)

def stacked(path, title, subtitle, bars, cats, note):
    """bars: (label, {cat: count}); cats: [(key, name, color_index)]"""
    W, L, Rm, top, bh, gap = 900, 250, 40, 100, 30, 24
    H = top + len(bars) * (bh + gap) + 56
    pw = W - L - Rm
    s = svg_open(W, H, title, subtitle)
    s += f'<text x="20" y="28" font-size="17" font-weight="600" fill="{INK}">{esc(title)}</text>\n'
    s += f'<text x="20" y="48" font-size="12.5" fill="{INK2}">{esc(subtitle)}</text>\n'
    x = L
    for key, name, ci in cats:
        s += f'<rect x="{x}" y="60" width="12" height="12" rx="3" fill="{C[ci]}"/><text x="{x+18}" y="70" font-size="12.5" fill="{INK2}">{esc(name)}</text>\n'
        x += 18 + 7 * len(name) + 24
    for i, (lab, d) in enumerate(bars):
        y = top + i * (bh + gap); tot = sum(d.get(k, 0) for k, _, _ in cats)
        s += f'<text x="{L-10}" y="{y+bh/2+4}" font-size="12.5" fill="{INK}" text-anchor="end">{esc(lab)} <tspan fill="{INK2}">{tot:,} loci</tspan></text>\n'
        cx = L
        for key, name, ci in cats:
            v = d.get(key, 0)
            if not v: continue
            w = pw * v / tot
            s += f'<rect x="{cx:.1f}" y="{y}" width="{max(w-2,1):.1f}" height="{bh}" rx="3" fill="{C[ci]}"><title>{esc(lab)} {esc(name)}: {v:,} ({100*v/tot:.1f}%)</title></rect>\n'
            if w > 60: s += f'<text x="{cx+(w-2)/2:.1f}" y="{y+bh/2+4}" font-size="12" font-weight="600" fill="{INK}" text-anchor="middle">{100*v/tot:.0f}%</text>\n'
            cx += w
    s += f'<text x="20" y="{H-8}" font-size="11.5" fill="{INK2}">{esc(note)}</text>\n</svg>\n'
    open(path, 'w').write(s)

# --- Ascomycota by class
asc = rd('2026-10-03_ascomycota_v060/by_class.tsv')
rows = []
for r in asc:
    n = int(r['genomes'])
    if n < 20: continue
    rt = dict(kv.split(':') for kv in r['routing'].split('; ') if ':' in kv)
    fb = int(rt.get('phylum_fallback', 0)) / n
    rows.append((r['class'] if r['class'] != 'NA' else 'class not assigned', n, float(r['pct_called']), 1 if fb > 0.5 else 0, fb))
rows.sort(key=lambda r: -r[2])
hbars(os.path.join(OUT, 'fig1_ascomycota_called_by_class.svg'), 'Ascomycota v0.6.0: genomes with a MAT call, by class',
      f'{sum(r[1] for r in rows):,} genomes in classes with 20 or more; "called" is detection, not ground truth',
      [(a, b, c, d) for a, b, c, d, _ in rows], [(0, 'routed by lineage'), (1, 'mostly phylum fallback (no record for the clade)')])

# --- Basidiomycota by order
bas = rd('2026-10-03_basidiomycota_v060/by_order.tsv')
npr = {r['order']: int(r['called_without_PR_only']) for r in rd('2026-10-03_basidiomycota_v060/by_order_nonPR.tsv')}
rows2 = []
for r in bas:
    n = int(r['n'])
    if n < 40: continue
    fb = int(r['phylum_fallback']) / n
    rows2.append((r['order'], n, float(r['pct']), 1 if fb > 0.5 else 0, 100 * npr[r['order']] / n))
rows2.sort(key=lambda r: -r[2])
hbars(os.path.join(OUT, 'fig2_basidiomycota_called_by_order.svg'), 'Basidiomycota v0.6.0: genomes with a MAT call, by order',
      f'{sum(r[1] for r in rows2):,} genomes in orders with 40 or more; bars include PR-only calls (1,360 are CAAX-only, unverified); the marker excludes them',
      rows2, [(0, 'routed by lineage'), (1, 'mostly phylum fallback (no record for the clade)')], tick_name='without PR-only calls')

# --- locus class composition
def classes(path, flt=None):
    c = collections.Counter()
    for r in rd(path):
        if flt and not flt(r): continue
        c[r['locus_class']] += 1
    return c
a_c = classes('2026-10-03_ascomycota_v060/loci.tsv'); b_c = classes('2026-10-03_basidiomycota_v060/loci.tsv')
m_c = classes('2026-10-03_mucoromycotina_mat/calls.tsv', lambda r: r['source'] == 'BFD' and r['status'] == 'called')
cats = [('mat_locus', 'full locus', 0), ('partial_locus', 'partial locus', 1), ('idiomorph_gene_only', 'gene only', 2), ('homothallic_candidate', 'homothallic cand.', 3)]
stacked(os.path.join(OUT, 'fig3_locus_class_by_campaign.svg'), 'What the calls contain: locus class by campaign',
        'Share of loci: full locus = MAT genes plus flank structure; gene only = the idiomorph gene without the locus',
        [('Ascomycota v0.6.0', a_c), ('Basidiomycota v0.6.0', b_c), ('Mucoromycotina BFD (calls)', m_c)], cats,
        'Segments under 5% are unlabelled; exact counts are in the analysis note table.')

# --- summary table
conf = lambda p, flt=None: collections.Counter(r['confidence'] for r in rd(p) if (not flt or flt(r)))
ca, cb = conf('2026-10-03_ascomycota_v060/loci.tsv'), conf('2026-10-03_basidiomycota_v060/loci.tsv')
cm = conf('2026-10-03_mucoromycotina_mat/calls.tsv', lambda r: r['source'] == 'BFD' and r['status'] == 'called')
with open(os.path.join(OUT, 'campaign_summary.tsv'), 'w') as o:
    o.write('campaign\tmetric\tvalue\n')
    for name, cc, lc in (('Ascomycota v0.6.0', ca, a_c), ('Basidiomycota v0.6.0', cb, b_c), ('Mucoromycotina BFD', cm, m_c)):
        for k in ('high', 'medium', 'low'): o.write(f'{name}\tconfidence_{k}\t{cc.get(k,0)}\n')
        for k, _, _ in cats: o.write(f'{name}\tlocus_class_{k}\t{lc.get(k,0)}\n')
print('Ascomycota classes', len(rows), 'Basidiomycota orders', len(rows2)); print(dict(a_c), dict(b_c), dict(m_c)); print(dict(ca), dict(cb), dict(cm))
