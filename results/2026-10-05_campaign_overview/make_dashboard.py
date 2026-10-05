"""Build dashboard.html (HTML fragment for the Artifact tool and for the repo) from the committed campaign tables.

Everything is drawn as inline SVG with CSS classes, so light and dark themes follow the page tokens.
Data: results/2026-10-03_{ascomycota,basidiomycota}_v060, 2026-10-03_mucoromycotina_mat, 2026-10-05_xylariales_nc1011_interval.
"""
import csv, collections, os, re, html, statistics as st
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, '..')
def rd(p): return list(csv.DictReader(open(os.path.join(R, p)), delimiter='\t'))
E = html.escape
def pct(a, b): return 100 * a / b

# ------------------------------------------------------------------ data
ASC = rd('2026-10-03_ascomycota_v060/genomes.tsv'); ASC_L = rd('2026-10-03_ascomycota_v060/loci.tsv'); ASC_C = rd('2026-10-03_ascomycota_v060/by_class.tsv')
BAS = rd('2026-10-03_basidiomycota_v060/genomes.tsv'); BAS_L = rd('2026-10-03_basidiomycota_v060/loci.tsv'); BAS_O = rd('2026-10-03_basidiomycota_v060/by_order.tsv')
NPR = {r['order']: int(r['called_without_PR_only']) for r in rd('2026-10-03_basidiomycota_v060/by_order_nonPR.tsv')}
MUC = [r for r in rd('2026-10-03_mucoromycotina_mat/calls.tsv') if r['source'] == 'BFD']
MUC_C = [r for r in MUC if r['status'] == 'called']
XG = {r['asm']: r for r in rd('2026-10-05_xylariales_nc1011_interval/genomes_miniprot.tsv')}
XGENES = collections.defaultdict(dict)
for r in rd('2026-10-05_xylariales_nc1011_interval/genes_miniprot.tsv'): XGENES[r['asm']][r['gene']] = (int(r['start']), int(r['end']), r['strand_norm'])

n_asc, n_bas, n_muc = len(ASC), len(BAS), len({r['genome'] for r in MUC})
called_asc = sum(r['status'] == 'called' for r in ASC); called_bas = sum(r['status'] == 'called' for r in BAS); called_muc = len({r['genome'] for r in MUC_C})
loci_asc, loci_bas, loci_muc = len(ASC_L), len(BAS_L), len(MUC_C)

# ------------------------------------------------------------------ svg helpers
def hbars(rows, tick=False, lab_w=170, w=640, rowh=22):
    """rows: (label, n, pct, series_index[, tick_pct])"""
    pw = w - lab_w - 150 if tick else w - lab_w - 56
    h = len(rows) * rowh + 34
    o = f'<svg viewBox="0 0 {w} {h}" class="chart" role="img">'
    for k in (0, 25, 50, 75, 100):
        x = lab_w + pw * k / 100
        o += f'<line x1="{x:.1f}" x2="{x:.1f}" y1="2" y2="{h-26}" class="grid"/><text x="{x:.1f}" y="{h-10}" class="axis" text-anchor="middle">{k}%</text>'
    for i, r in enumerate(rows):
        lab, n, p, si = r[:4]; t = r[4] if len(r) > 4 else None
        y = 6 + i * rowh
        o += f'<text x="{lab_w-8}" y="{y+11}" class="lab" text-anchor="end">{E(lab)} <tspan class="dim">{n:,}</tspan></text>'
        o += f'<rect x="{lab_w}" y="{y}" width="{max(pw*p/100,1.5):.1f}" height="12" rx="3" class="s{si}"><title>{E(lab)}: {p:.1f}% of {n:,} genomes</title></rect>'
        if t is not None:
            tx = lab_w + pw * t / 100
            o += f'<line x1="{tx:.1f}" x2="{tx:.1f}" y1="{y-3}" y2="{y+15}" class="tick"><title>{E(lab)}: {t:.1f}% without PR-only calls</title></line>'
        o += f'<text x="{lab_w+pw*p/100+6:.1f}" y="{y+11}" class="val">{p:.1f}%' + (f'<tspan class="dim"> ({t:.0f}% w/o PR-only)</tspan>' if t is not None and abs(t - p) >= 3 else '') + '</text>'
    return o + '</svg>'

def stacked(rows, cats, w=640, lab_w=170, bh=24, gap=14):
    """rows: (label, {key: count}); cats: [(key, name, series)]"""
    pw = w - lab_w - 16; h = len(rows) * (bh + gap) + 8
    o = f'<svg viewBox="0 0 {w} {h}" class="chart" role="img">'
    for i, (lab, d) in enumerate(rows):
        y = 4 + i * (bh + gap); tot = sum(d.get(k, 0) for k, _, _ in cats) or 1
        o += f'<text x="{lab_w-8}" y="{y+bh/2+4}" class="lab" text-anchor="end">{E(lab)} <tspan class="dim">{tot:,}</tspan></text>'
        cx = lab_w
        for k, name, si in cats:
            v = d.get(k, 0)
            if not v: continue
            sw = pw * v / tot
            o += f'<rect x="{cx:.1f}" y="{y}" width="{max(sw-2,1.5):.1f}" height="{bh}" rx="3" class="s{si}"><title>{E(lab)} · {E(name)}: {v:,} ({100*v/tot:.1f}%)</title></rect>'
            if sw > 46: o += f'<text x="{cx+(sw-2)/2:.1f}" y="{y+bh/2+4}" class="{"segd" if si in (3, 4) else "seg"}" text-anchor="middle">{100*v/tot:.0f}%</text>'
            cx += sw
    return o + '</svg>'

def legend(items):
    return '<div class="legend">' + ''.join(f'<span><i class="sw s{si}"></i>{E(n)}</span>' for si, n in items) + '</div>'

# ------------------------------------------------------------------ tab 1: summary
asc_rows = []
for r in ASC_C:
    n = int(r['genomes'])
    if n < 20: continue
    rt = dict(kv.split(':') for kv in r['routing'].split('; ') if ':' in kv)
    asc_rows.append((r['class'] if r['class'] != 'NA' else 'class not assigned', n, float(r['pct_called']), 2 if int(rt.get('phylum_fallback', 0)) / n > .5 else 1))
asc_rows.sort(key=lambda r: -r[2])
bas_rows = []
for r in BAS_O:
    n = int(r['n'])
    if n < 40: continue
    bas_rows.append((r['order'], n, float(r['pct']), 2 if int(r['phylum_fallback']) / n > .5 else 1, pct(NPR[r['order']], n)))
bas_rows.sort(key=lambda r: -r[2])
loc_cats = [('mat_locus', 'full locus', 1), ('partial_locus', 'partial locus', 2), ('idiomorph_gene_only', 'gene only', 3), ('homothallic_candidate', 'homothallic candidate', 4)]
cnt = lambda rows, key='locus_class': collections.Counter(r[key] for r in rows)
loc_rows = [('Ascomycota', cnt(ASC_L)), ('Basidiomycota', cnt(BAS_L)), ('Mucoromycotina', cnt(MUC_C))]
conf = lambda rows: collections.Counter(r['confidence'] for r in rows)
ca, cb, cm = conf(ASC_L), conf(BAS_L), conf(MUC_C)

def tile(k, v, sub): return f'<div class="tile"><div class="k">{E(k)}</div><div class="v">{v}</div><div class="s">{E(sub)}</div></div>'
tab_summary = f'''
<div class="tiles">
{tile('Genomes searched', f'{n_asc+n_bas+n_muc:,}', f'Ascomycota {n_asc:,} · Basidiomycota {n_bas:,} · Mucoromycotina {n_muc}')}
{tile('MAT loci reported', f'{loci_asc+loci_bas+loci_muc:,}', f'{loci_asc:,} · {loci_bas:,} · {loci_muc}')}
{tile('Genomes with a call', f'{pct(called_asc,n_asc):.0f}% · {pct(called_bas,n_bas):.0f}% · {pct(called_muc,n_muc):.0f}%', 'Ascomycota · Basidiomycota · Mucoromycotina; detection, not ground truth')}
{tile('High-confidence loci', f"{pct(ca['high'],loci_asc):.0f}% · {pct(cb['high'],loci_bas):.0f}% · {pct(cm['high'],loci_muc):.0f}%", 'same order of campaigns')}
</div>
<section><h2>Call rate by Ascomycota class</h2>
<p class="note">Classes with 20 or more genomes ({sum(r[1] for r in asc_rows):,} of {n_asc:,}). Orange classes had no curated record for their clade and were searched against the whole phylum, so a low rate there is a reference gap, not a measured absence.</p>
{legend([(1,'routed by lineage'),(2,'mostly phylum fallback')])}{hbars(asc_rows)}</section>
<section><h2>Call rate by Basidiomycota order</h2>
<p class="note">Orders with 40 or more genomes. The marker is the rate after removing PR-only calls; 1,360 PR calls come only from a CAAX-motif scan and are labelled unverified. Cantharellales falls from 91% to 17% and Trichosporonales from 47% to 8%.</p>
{legend([(1,'routed by lineage'),(2,'mostly phylum fallback')])}{hbars(bas_rows, tick=True)}</section>
<section><h2>What the calls contain</h2>
<p class="note">Share of loci by class. Half of the Basidiomycota loci are the idiomorph gene alone, without a recognised locus around it.</p>
{legend([(c[2], c[1]) for c in loc_cats])}{stacked(loc_rows, loc_cats)}
<details><summary>Table view</summary><table><thead><tr><th>Campaign</th>{''.join(f'<th>{E(c[1])}</th>' for c in loc_cats)}<th>High</th><th>Medium</th><th>Low</th></tr></thead><tbody>
{''.join(f"<tr><td>{n}</td>" + ''.join(f"<td>{d.get(c[0],0):,}</td>" for c in loc_cats) + f"<td>{cf['high']:,}</td><td>{cf['medium']:,}</td><td>{cf['low']:,}</td></tr>" for (n, d), cf in zip(loc_rows, (ca, cb, cm)))}
</tbody></table></details></section>'''

# ------------------------------------------------------------------ tab 2: exemplars and outliers
def card(kind, title, body, src): return f'<article class="card {kind}"><span class="tag">{"Exemplar" if kind=="ex" else "Outlier"}</span><h3>{E(title)}</h3><p>{body}</p><p class="src">{E(src)}</p></article>'
dothi = [r for r in ASC_L if r['family_called'] == 'Ascomycota:MAT' and r['locus_class'] == 'mat_locus' and r['class_'] == 'Dothideomycetes']
dboth = pct(sum('APN2' in r['genes_found'] and 'SLA2' in r['genes_found'] for r in dothi), len(dothi))
cards = [
 card('ex', 'Wallemiales: 0 to 51 of 51', 'No MAT record existed. One record built from an undeposited locus calls every one of the 51 genomes (33 high, 18 medium) and separates two versions of the locus; the four W. ichthyophaga genomes in the minority are the published inverted strains.', 'docs/notes/publication-highlights.md'),
 card('ex', 'Rhodotorula P/R: 61 of 62 held-out strains', 'The called P/R allele matches an independent study for 61 of 62 held-out strains; the miss is a reported hybrid. HD calls rose from 0 to 217 genomes.', 'docs/publication-notable-findings/002'),
 card('ex', 'Lineage-routed orders at 95 to 100%', 'Ustilaginales 100%, Wallemiales 100%, Sporidiobolales 99.6%, Polyporales 99.5%, Agaricales 98.8%: orders with their own curated record.', 'results/2026-10-03_basidiomycota_v060/by_order.tsv'),
 card('ex', 'Zygomycete benchmark 23 of 23', 'All 23 Zygo genomes are typed on both scaffold and contig input, and the regression panel shows no changed locus across five panels.', 'results/2026-10-03_regression_v060'),
 card('out', 'Reference gaps: Orbiliomycetes 9%, Dipodascomycetes 29%', 'Both classes are searched through the phylum fallback. Low rates mark missing curation.', 'results/2026-10-03_ascomycota_v060/by_class.tsv'),
 card('out', 'PR-only calls inflate fallback orders', 'Cantharellales 91% drops to 17% and Trichosporonales 47% to 8% without PR-only calls. The CAAX-only calls are under a chance-level test in a draft PR (#33): about 13% of the unverified calls sit at chance.', 'analysis/2026-10-05_caax-receptor-test.md (draft)'),
 card('out', f'Dothideomycetes: SLA2 detached ({dboth:.0f}% with both flanks)', f'Only {dboth:.0f}% of full Dothideomycete loci carry both APN2 and SLA2, against 82 to 91% in the other filamentous classes. A genome-wide search finds SLA2 in every sampled genome but inside the locus in only 14% (80% in Sordariomycetes); Cladosporiales are the exception. Mostly organisation, not missed detection.', 'analysis/2026-10-05_dothideomycetes-sla2.md'),
 card('out', 'Early-diverging fungi: few calls', 'Mortierellomycota 6 of 100, Kickxellomycota 7 of 190, and Lichtheimiaceae 1 of 30 in the 2026-09-26 scan. No Mucorales-type synteny support in Mortierella or Kickxella.', 'results/2026-09-26_early_diverging/summary.txt'),
 card('out', 'Xylaria NC1011: flanks, no MAT gene', 'SLA2, COX13 and APN2 are intact and expressed, with 475 and 118 bp between them and no HMG-box or alpha-box domain in the 47-kb region. The one call is an HMG gene nearest NCU03481 (14% support).', 'analysis/2026-10-05_xylariales-synteny.md'),
 card('out', 'Mycotypha: sexM 153 kb from sexP', 'A homothallic species with its two idiomorph genes about 150 kb apart on one scaffold. v0.6.0 did not call the sexM gene and the classifier carried it among its negatives; the classifier was rebuilt on 2026-10-04 and the new outcome has not been re-checked here.', 'analysis/2026-10-04_unclassified-hmg-in-gene-trees.md'),
]
tab_ex = '<div class="cards">' + ''.join(cards) + '</div>'

# ------------------------------------------------------------------ tab 3: Xylariales synteny
def parse(s): return [(m.group(1), m.group(2)) for m in re.finditer(r'(\w+)\(([+-])\)', s)]
def renorm(g):
    d = dict(g); return [(n, '+' if s == '-' else '-') for n, s in reversed(g)] if d.get('SLA2') == '-' else g
def state(r):
    g = renorm(parse(r['order_string'])); names = [n for n, s in g]; d = dict(g)
    if not all(k in d for k in ('SLA2', 'COX13', 'APN2')): return 'incomplete'
    t = [n for n in names if n in ('SLA2', 'COX13', 'APN2')]
    if t == ['SLA2', 'APN2', 'COX13'] and d['APN2'] == '-' and d['COX13'] == '+': return 'SAC'
    if t == ['SLA2', 'COX13', 'APN2'] and d['COX13'] == '-' and d['APN2'] == '+':
        i, j = names.index('SLA2'), max(names.index('COX13'), names.index('APN2'))
        return 'SCA+' if [n for n in names[i+1:j] if n not in ('COX13', 'APN2')] else 'SCA'
    return 'other'
X = [r for r in XG.values() if r['order'] == 'Xylariales']
for r in X: r['state'] = state(r)
fam_rows = []
fams = collections.defaultdict(collections.Counter)
for r in X: fams[r['family'] or 'family not assigned'][r['state']] += 1
for f, c in sorted(fams.items(), key=lambda x: -sum(x[1].values())):
    fam_rows.append((f, c))
st_cats = [('SAC', 'outgroup-like order (MAT interval kept)', 1), ('SCA', 'COX13 beside SLA2', 2), ('SCA+', 'COX13 beside SLA2, plus H609 inside', 3), ('other', 'other order', 4)]
sp_state = collections.defaultdict(collections.Counter)
for r in X: sp_state[r['species']][r['state']] += 1
sp_maj = collections.Counter(c.most_common(1)[0][0].replace('SCA+', 'SCA') for c in sp_state.values())
n_sp = len(sp_state); eutypa = sum(1 for r in X if r['species'] == 'Eutypa lata')

# locus maps
EX = [('GCA_000213195.1_v1.0', 'Neurospora tetrasperma', 'outgroup, Sordariales'), ('GCA_000259975.2_FO_MN25_V1', 'Fusarium oxysporum', 'outgroup, Hypocreales'),
      ('GCF_022984875.1_Hypfra2', 'Hypoxylon fragiforme', 'SAC, 12.5-kb interval'), ('GCA_000966885.1_ASM96688v1', 'Xylaria sp. JS573', 'SAC, 7.5-kb interval'),
      ('GCA_022478725.1_Dales1', 'Daldinia eschscholtzii', 'SAC, interval shrunk to 2.9 kb'), ('GCA_022453505.1_Xylcub1', 'Xylaria flabelliformis NC1011', 'SCA, 2.2 kb'),
      ('GCF_000349385.1_UCREL1V03', 'Eutypa lata UCREL1', 'SCA with H609 inside, 22.5 kb')]
GCLS = {'SLA2': 'g1', 'COX13': 'g2', 'APN2': 'g3', 'APC5': 'g5', 'CIA30': 'g5', 'H609': 'g4'}
def frame(asm):
    r = XG[asm]; g = XGENES[asm]; order = [n for n, s in parse(r['order_string'])]
    by_start = [n for n, v in sorted(g.items(), key=lambda kv: kv[1][0])]
    flipped = len(by_start) > 1 and by_start == order[::-1]
    out = {}
    for n, (a, b, s) in g.items():
        out[n] = (-b, -a, s) if flipped else (a, b, s)
    if out['SLA2'][2] == '-':
        out = {n: (-b, -a, '+' if s == '-' else '-') for n, (a, b, s) in out.items()}
    x0 = out['SLA2'][0]
    return {n: (a - x0, b - x0, s) for n, (a, b, s) in out.items()}
def locus_maps():
    lo, hi = -22000, 52000; w = 760; lab_w = 210; pw = w - lab_w - 10; rowh = 46
    sx = lambda x: lab_w + pw * (min(max(x, lo), hi) - lo) / (hi - lo)
    h = len(EX) * rowh + 40
    o = f'<svg viewBox="0 0 {w} {h}" class="chart wide" role="img">'
    for k in range(-20, 51, 10):
        x = sx(k * 1000); o += f'<line x1="{x:.1f}" x2="{x:.1f}" y1="4" y2="{h-26}" class="grid"/><text x="{x:.1f}" y="{h-10}" class="axis" text-anchor="middle">{k:+d} kb</text>'
    for i, (asm, name, note) in enumerate(EX):
        y = 8 + i * rowh; g = frame(asm); mid = y + 14
        o += f'<text x="{lab_w-10}" y="{y+12}" class="lab" text-anchor="end">{E(name)}</text><text x="{lab_w-10}" y="{y+26}" class="dim" text-anchor="end" font-size="11">{E(note)}</text>'
        o += f'<line x1="{lab_w}" x2="{w-10}" y1="{mid}" y2="{mid}" class="axisline"/>'
        sla, cox, apn = g['SLA2'], g.get('COX13'), g.get('APN2')
        if apn and sla and (apn[0] - sla[1]) > 3500 and apn[0] > sla[1] and (not cox or cox[0] > apn[0] or True):
            inner_start = sla[1]; inner_end = min(apn[0], cox[0]) if cox and cox[0] > sla[1] else apn[0]
            if inner_end - inner_start > 3500:
                o += f'<rect x="{sx(inner_start):.1f}" y="{mid-12}" width="{sx(inner_end)-sx(inner_start):.1f}" height="24" rx="3" class="gap"><title>{name}: {(inner_end-inner_start)/1000:.1f} kb between SLA2 and the next flank gene</title></rect>'
        for n, (a, b, s) in sorted(g.items(), key=lambda kv: kv[1][0]):
            if b < lo or a > hi: continue
            x1, x2 = sx(a), sx(b); cls = GCLS.get(n, 'g0')
            tip = max(min(7, (x2 - x1) * .5), 3)
            pts = f'{x1:.1f},{mid-8} {x2-tip if s=="+" else x2:.1f},{mid-8} {x2:.1f},{mid} {x2-tip if s=="+" else x2:.1f},{mid+8} {x1:.1f},{mid+8}' if s == '+' else f'{x2:.1f},{mid-8} {x1+tip:.1f},{mid-8} {x1:.1f},{mid} {x1+tip:.1f},{mid+8} {x2:.1f},{mid+8}'
            o += f'<polygon points="{pts}" class="{cls}"><title>{n} ({s}) {a/1000:.1f} to {b/1000:.1f} kb</title></polygon>'
            if n == 'SLA2' and (x2 - x1) > 26: o += f'<text x="{(x1+x2)/2:.1f}" y="{mid+3.5}" class="gl" text-anchor="middle">SLA2</text>'
    return o + '</svg>'

# adjacency
ADJ = [r for r in rd('2026-10-05_xylariales_nc1011_interval/adjacency_miniprot.tsv')]
def grp(a):
    r = XG[a['asm']]
    if r['order'] == 'Xylariales': return 'Xylariales ' + ('SAC' if r['state'] == 'SAC' else 'SCA' if r['state'].startswith('SCA') else 'other')
    if r['group'] == 'sister': return 'Amphisphaeriales'
    if r['order'] in ('Sordariales', 'Hypocreales', 'Glomerellales', 'Diaporthales', 'Coniochaetales'): return 'outgroups'
    return None
for r in X: pass
AT = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
hdr = ['ARM-FRE', 'FRE-APC5', 'APC5-CIA30', 'CIA30-SLA2', 'SLA2-COX13', 'COX13-APN2']
for a in ADJ:
    g = grp(a)
    if not g: continue
    for h in hdr:
        v = a[h]
        if v == 'NA': continue
        AT[g][h][1] += 1; AT[g][h][0] += v == '1'
gs = ['outgroups', 'Amphisphaeriales', 'Xylariales SAC', 'Xylariales SCA']
def heat(v):
    if v is None: return '<td class="na">n/a</td>'
    return f'<td class="hm" style="--p:{v:.0f}">{v:.0f}%</td>'
adj_tbl = '<table class="heat"><thead><tr><th>Neighbours in NC1011</th>' + ''.join(f'<th>{g}</th>' for g in gs) + '</tr></thead><tbody>' + ''.join(
    f'<tr class="{"key" if h in ("CIA30-SLA2","SLA2-COX13") else ""}"><th>{h}</th>' + ''.join(heat(pct(*AT[g][h]) if AT[g][h][1] else None) for g in gs) + '</tr>' for h in hdr) + '</tbody></table>'
# gap histogram
bins = [(0, 1000, '<1'), (1000, 3000, '1-3'), (3000, 6000, '3-6'), (6000, 10000, '6-10'), (10000, 15000, '10-15'), (15000, 25000, '15-25'), (25000, 10**9, '25+')]
def hist(sel):
    c = [0] * len(bins); n = 0
    for r in X:
        if sel(r) and r['sla2_apn2_gap_bp']:
            gp = int(r['sla2_apn2_gap_bp']); n += 1
            for i, (a, b, _) in enumerate(bins):
                if a <= gp < b: c[i] += 1; break
    return c, n
hs, ns = hist(lambda r: r['state'] == 'SAC'); hc, nc = hist(lambda r: r['state'].startswith('SCA') and r['species'] != 'Eutypa lata')
he, ne = hist(lambda r: r['species'] == 'Eutypa lata')
def grouped_hist():
    w, h, L, B = 640, 220, 44, 34; pw = w - L - 10; ph = h - B - 14; m = -(-max(max(hs), max(hc), max(he)) // 30) * 30 or 30
    o = f'<svg viewBox="0 0 {w} {h}" class="chart" role="img">'
    for k in range(0, 5):
        v = m * k / 4; y = 8 + ph - ph * k / 4
        o += f'<line x1="{L}" x2="{w-8}" y1="{y:.1f}" y2="{y:.1f}" class="grid"/><text x="{L-6}" y="{y+4:.1f}" class="axis" text-anchor="end">{v:.0f}</text>'
    slot = pw / len(bins)
    for i, (_, _, lab) in enumerate(bins):
        x0 = L + slot * i
        for j, (vals, si, nm) in enumerate(((hs, 1, 'SAC'), (hc, 2, 'SCA (excluding E. lata)'), (he, 3, 'Eutypa lata'))):
            bw = slot / 4; bx = x0 + slot * .12 + j * (bw + 2); bhh = ph * vals[i] / m
            o += f'<rect x="{bx:.1f}" y="{8+ph-bhh:.1f}" width="{bw:.1f}" height="{max(bhh,0):.1f}" rx="2" class="s{si}"><title>{nm}: {vals[i]} genomes with {lab} kb</title></rect>'
        o += f'<text x="{x0+slot/2:.1f}" y="{h-18}" class="axis" text-anchor="middle">{lab}</text>'
    o += f'<text x="{L+pw/2}" y="{h-3}" class="axis" text-anchor="middle">SLA2 to APN2 distance (kb)</text>'
    return o + '</svg>'

tab_syn = f'''
<section><h2>Gene order near the MAT region in Xylariales</h2>
<p class="note">Seven genomes drawn on the same frame, SLA2 at zero and pointing right. Gold boxes mark the interval between SLA2 and the next flank gene where a MAT locus sits in outgroups. In outgroups and the SAC genomes the interval is open; in NC1011 it holds only COX13.</p>
<div class="legend"><span><i class="sw s1"></i>SLA2</span><span><i class="sw s2"></i>COX13</span><span><i class="sw s3"></i>APN2</span><span><i class="sw s4"></i>H609</span><span><i class="sw sg5"></i>APC5, CIA30</span><span><i class="sw sg0"></i>other neighbours</span><span><i class="sw sgap"></i>interval ≥ 3.5 kb</span></div>
{locus_maps()}</section>
<section><h2>Sampling: 257 genomes, {n_sp} species</h2>
<p class="note">{eutypa} of the 257 genomes are strains of <em>Eutypa lata</em>, all in the SCA state with H609 inside the interval, so genome counts overstate how often SCA arose. By species, SAC is {sp_maj['SAC']} and SCA is {sp_maj['SCA']} (majority state per species; {sp_maj['other']} other, {sp_maj['incomplete']} incomplete).</p>
{legend([(c[2], c[1]) for c in st_cats])}{stacked(fam_rows, st_cats)}</section>
<section><h2>Which neighbours stay together</h2>
<p class="note">Share of genomes in which two NC1011 neighbours are still adjacent with the same relative orientation. Two junctions change: CIA30 next to SLA2 (absent in outgroups, shared with Amphisphaeriales), and SLA2 next to COX13 (only in SCA). COX13 and APN2 never separate.</p>
{adj_tbl}</section>
<section><h2>How much room is there between SLA2 and APN2</h2>
<p class="note">SAC keeps a gap of several kb, enough for a MAT locus; SCA genomes mostly have only COX13 between them. The deletions in a few SAC genomes (about 3 kb) are a second route to losing the region.</p>
{legend([(1,f'SAC, n={ns}'),(2,f'SCA excluding E. lata, n={nc}'),(3,f'Eutypa lata, n={ne}')])}{grouped_hist()}</section>
<p class="note">Sources: results/2026-10-05_xylariales_nc1011_interval, analysis/2026-10-05_xylariales-synteny.md. Gene-level only; breakpoints are bounded by gene ends.</p>'''

# ------------------------------------------------------------------ SLA2 position test (Dothideomycetes)
SLA2T = [r for r in rd('2026-10-05_dothideo_sla2_test/sla2_distance.tsv') if r['gene'] == 'SLA2']
def sla2_cat(r):
    if r['where'] == 'within_locus': return 'in'
    if r['where'] == 'same_contig':
        d = int(r['dist_best_bp'] or 0)
        return 'in' if d == 0 else 'near' if d <= 100000 else 'far'
    return 'other'
sla_rows = []
def add_row(label, sel):
    c = collections.Counter(sla2_cat(r) for r in sel)
    if sum(c.values()): sla_rows.append((label, dict(c)))
add_row('Sordariomycetes (control)', [r for r in SLA2T if r['class'] == 'Sordariomycetes'])
dothi = [r for r in SLA2T if r['class'] == 'Dothideomycetes']
add_row('Dothideomycetes, all sampled', dothi)
bo = collections.defaultdict(list)
for r in dothi: bo[r['order'] or 'order not assigned'].append(r)
for o, v in sorted(bo.items(), key=lambda x: -len(x[1])):
    if len(v) >= 20: add_row(o, v)
sla_cats = [('in', 'inside the called locus', 1), ('near', '20 to 100 kb away, same contig', 3), ('far', 'over 100 kb away, same contig', 4), ('other', 'on another contig', 2)]
n_dothi_sla = len(dothi); pos_med = sorted(float(r['positives']) for r in dothi if r['where'] != 'no_hit')[len([1 for r in dothi if r['where'] != 'no_hit']) // 2]
n_hit = sum(r['where'] != 'no_hit' for r in dothi)

# ------------------------------------------------------------------ tab 4: size and content
muc = collections.defaultdict(list)
for r in MUC_C:
    try: sz = int(r['locus_end']) - int(r['locus_start'])
    except ValueError: continue
    muc[r['name'].split()[0]].append((sz, set(r['genes_found'].split(',')), r['locus_class']))
muc_rows = sorted(((g, v) for g, v in muc.items() if len(v) >= 5), key=lambda x: -len(x[1]))
def range_chart(rows, w=640, lab_w=170, xmax=None, unit='kb'):
    """rows: (label, n, min, q1, med, q3, max)"""
    xmax = xmax or max(r[6] for r in rows) * 1.05
    pw = w - lab_w - 20; rowh = 22; h = len(rows) * rowh + 48
    sx = lambda v: lab_w + pw * min(v, xmax) / xmax
    o = f'<svg viewBox="0 0 {w} {h}" class="chart" role="img">'
    step = 25 if xmax > 100 else 10 if xmax > 40 else 5
    for k in range(0, int(xmax) + 1, step):
        x = sx(k); o += f'<line x1="{x:.1f}" x2="{x:.1f}" y1="2" y2="{h-38}" class="grid"/><text x="{x:.1f}" y="{h-24}" class="axis" text-anchor="middle">{k}</text>'
    o += f'<text x="{lab_w+pw/2}" y="{h-6}" class="axis" text-anchor="middle">locus length ({unit}); line shows the full range</text>'
    for i, (lab, n, mn, q1, md, q3, mx) in enumerate(rows):
        y = 6 + i * rowh; m = y + 6
        o += f'<text x="{lab_w-8}" y="{y+10}" class="lab" text-anchor="end">{E(lab)} <tspan class="dim">{n}</tspan></text>'
        o += f'<line x1="{sx(mn):.1f}" x2="{sx(mx):.1f}" y1="{m}" y2="{m}" class="whisk"/>' + (f'<text x="{sx(xmax)+2:.1f}" y="{m+4}" class="axis">›<title>{E(lab)}: maximum {mx:.0f} kb, beyond the axis</title></text>' if mx > xmax else '')
        o += f'<rect x="{sx(q1):.1f}" y="{m-5}" width="{max(sx(q3)-sx(q1),2):.1f}" height="10" rx="2" class="s1"><title>{E(lab)}: median {md:.1f}, IQR {q1:.1f} to {q3:.1f}, range {mn:.1f} to {mx:.1f} {unit} (n={n})</title></rect>'
        o += f'<circle cx="{sx(md):.1f}" cy="{m}" r="3.6" class="med"/>'
    return o + '</svg>'
def quart(v):
    v = sorted(v); n = len(v); return (v[0], v[n//4], st.median(v), v[(3*n)//4], v[-1])
muc_chart_rows = [(g, len(v)) + tuple(x / 1000 for x in quart([a for a, _, _ in v])) for g, v in muc_rows]
muc_chart_rows = [(a, b, c, d, e, f, g_) for a, b, c, d, e, f, g_ in muc_chart_rows]
genes = ['sexP', 'sexM', 'tptA', 'rnhA', 'algA', 'glrA', 'btbA']
mt = '<table class="heat"><thead><tr><th>Genus (loci)</th>' + ''.join(f'<th>{g}</th>' for g in genes) + '</tr></thead><tbody>' + ''.join(
    f'<tr><th>{E(g)} <span class="dim">{len(v)}</span></th>' + ''.join(heat(pct(sum(gn in s for _, s, _ in v), len(v))) for gn in genes) + '</tr>' for g, v in muc_rows) + '</tbody></table>'
fil = collections.defaultdict(list)
for r in ASC_L:
    if r['family_called'] == 'Ascomycota:MAT' and r['locus_class'] == 'mat_locus':
        fil[r['class_']].append((int(r['end']) - int(r['start']), 'APN2' in r['genes_found'] and 'SLA2' in r['genes_found'], len(r['genes_found'].split('|'))))
fcl = ['Sordariomycetes', 'Eurotiomycetes', 'Dothideomycetes', 'Leotiomycetes', 'Lecanoromycetes']
fil_rows = [(c, len(fil[c])) + tuple(x / 1000 for x in quart([a for a, _, _ in fil[c]])) for c in fcl]
fil_flank = [(c, len(fil[c]), pct(sum(b for _, b, _ in fil[c]), len(fil[c])), 1 if c != 'Dothideomycetes' else 2) for c in fcl]
tab_size = f'''
<section><h2>Mucoromycotina: locus length by genus</h2>
<p class="note">Genera with five or more called loci in the BFD set. Bar is the interquartile range, dot the median, line the full range. Median across all {len(MUC_C)} loci is 14.5 kb. <em>Umbelopsis</em> has the longest loci (median 41.5 kb); <em>Syncephalastrum</em> and <em>Apophysomyces</em> have the shortest (5.3 and 7.8 kb).</p>
{range_chart(muc_chart_rows, xmax=120)}</section>
<section><h2>Mucoromycotina: which genes the locus carries</h2>
<p class="note">Share of each genus's loci in which the gene was found. Content differs by genus: btbA is found only in <em>Rhizopus</em> (77% of its loci); glrA is missing from <em>Cunninghamella</em>, <em>Syncephalastrum</em>, <em>Apophysomyces</em> and <em>Phycomyces</em>; rnhA is missing from <em>Backusella</em> and <em>Apophysomyces</em>; tptA is missing from <em>Syncephalastrum</em>. Percentages for sexP and sexM can add to more than 100 because a locus can list both.</p>{mt}</section>
<section><h2>Filamentous Ascomycota: locus length by class</h2>
<p class="note">Full loci of the generic Ascomycota MAT family. Sordariomycetes loci are twice as long as Dothideomycete loci at the median (17.1 against 8.9 kb), and Dothideomycetes have a long upper tail.</p>
{range_chart(fil_rows, xmax=100)}</section>
<section><h2>Filamentous Ascomycota: APN2 and SLA2 in the same locus</h2>
<p class="note">Share of full loci that carry both flank genes. Dothideomycetes stand out.</p>
{hbars([(c, n, p, s) for c, n, p, s in fil_flank])}
<p class="note">Sources: results/2026-10-03_mucoromycotina_mat/calls.tsv (BFD, called), results/2026-10-03_ascomycota_v060/loci.tsv.</p></section>
<section><h2>Dothideomycetes: SLA2 is detached from the MAT locus, not missed</h2>
<p class="note">A genome-wide search for SLA2 in {n_dothi_sla} sampled Dothideomycete genomes finds it in {n_hit} (best hit {pos_med*100:.0f}% positives against the NC1011 protein), so the gene is present and recognisable. Only a minority place it inside the called locus, against 80% in the Sordariomycete control. Cladosporiales keep it there; in most other orders SLA2 sits tens of kb away or on another contig. "Another contig" includes assembly breaks and misses (17% in the control). Sample: up to 25 genomes per order, one locus each.</p>
{legend([(c[2], c[1]) for c in sla_cats])}{stacked(sla_rows, sla_cats, lab_w=235)}
<p class="note">Source: results/2026-10-05_dothideo_sla2_test/sla2_distance.tsv (miniprot, SLA2 query from Xylaria NC1011), analysis/2026-10-05_dothideomycetes-sla2.md.</p></section>'''

# ------------------------------------------------------------------ page
CSS = '''
:root{--bg:#f4f5f4;--panel:#fcfcfb;--ink:#101312;--ink2:#505451;--dim:#7b807d;--grid:#e3e5e2;--line:#cfd2cf;--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a;--s4:#eda100;--g5:#8a98b0;--g0:#c4c9c6;--gap:#f3d98a;--gapline:#c98500;--accent:#2a78d6;--ok:#1a7f55;--bad:#b4442a;--font:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;--mono:"IBM Plex Mono",ui-monospace,Menlo,monospace}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#121413;--panel:#1a1c1b;--ink:#f1f2f0;--ink2:#c3c6c2;--dim:#979c98;--grid:#2c2f2d;--line:#3a3e3b;--s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--g5:#6f7d96;--g0:#4a504d;--gap:#4a3e12;--gapline:#c98500;--accent:#3987e5;--ok:#4cc79a;--bad:#ee8a70;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#121413;--panel:#1a1c1b;--ink:#f1f2f0;--ink2:#c3c6c2;--dim:#979c98;--grid:#2c2f2d;--line:#3a3e3b;--s1:#3987e5;--s2:#d95926;--s3:#199e70;--s4:#c98500;--g5:#6f7d96;--g0:#4a504d;--gap:#4a3e12;--gapline:#c98500;--accent:#3987e5;--ok:#4cc79a;--bad:#ee8a70;color-scheme:dark}
*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);font:14px/1.5 var(--font);margin:0}
.wrap{max-width:1040px;margin:0 auto;padding:20px 16px 48px}
header h1{font-size:22px;line-height:1.2;margin:0 0 4px;font-weight:600;letter-spacing:-.01em;text-wrap:balance}
header p{margin:0;color:var(--ink2);max-width:70ch}
nav{display:flex;gap:4px;flex-wrap:wrap;margin:18px 0 18px;border-bottom:1px solid var(--line)}
nav button{font:inherit;background:none;border:0;border-bottom:2px solid transparent;color:var(--ink2);padding:8px 12px;cursor:pointer;margin-bottom:-1px}
nav button:hover{color:var(--ink)}
nav button[aria-selected=true]{color:var(--ink);border-bottom-color:var(--accent);font-weight:600}
nav button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.tab{display:none}.tab.on{display:block}
section{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:16px 16px 14px;margin:0 0 14px;min-width:0}
h2{font-size:15px;margin:0 0 4px;font-weight:600}
.note{color:var(--ink2);margin:0 0 10px;max-width:78ch;font-size:13px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:10px;margin:0 0 14px}
.tile{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px 14px;min-width:0}
.tile .k{font-size:12px;color:var(--ink2);text-transform:uppercase;letter-spacing:.06em}
.tile .v{font:600 22px/1.25 var(--mono);margin:4px 0 2px;font-variant-numeric:tabular-nums}
.tile .s{font-size:12px;color:var(--dim)}
.chart{width:100%;height:auto;display:block}
.chart.wide{min-width:560px}
section>svg,section>div:has(svg){overflow-x:auto}
.lab{fill:var(--ink);font-size:12px}.dim{fill:var(--dim);color:var(--dim)}.val{fill:var(--ink);font-size:11.5px;font-variant-numeric:tabular-nums}
.axis{fill:var(--ink2);font-size:11px}.grid{stroke:var(--grid);stroke-width:1}.axisline{stroke:var(--line);stroke-width:1}
.tick{stroke:var(--ink);stroke-width:2.5}.whisk{stroke:var(--ink2);stroke-width:1.5}.med{fill:var(--ink)}
.seg{fill:#fff;font-size:11.5px;font-weight:600}
.s1{fill:var(--s1);background:var(--s1)}.s2{fill:var(--s2);background:var(--s2)}.s3{fill:var(--s3);background:var(--s3)}.s4{fill:var(--s4);background:var(--s4)}
.segd{fill:#0b0b0b;font-size:11.5px;font-weight:600}
.g1{fill:var(--s1)}.g2{fill:var(--s2)}.g3{fill:var(--s3)}.g4{fill:var(--s4)}.g5{fill:var(--g5)}.g0{fill:var(--g0)}.gap{fill:var(--gap);stroke:var(--gapline);stroke-width:1;stroke-dasharray:3 2}
.gl{fill:#fff;font-size:10px;font-weight:600}.gl2{fill:var(--ink2);font-size:9.5px}
.legend{display:flex;gap:6px 16px;flex-wrap:wrap;font-size:12px;color:var(--ink2);margin:0 0 8px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.sw{width:11px;height:11px;border-radius:3px;display:inline-block}
.sg5{background:var(--g5)}.sg0{background:var(--g0)}.sgap{background:var(--gap);outline:1px dashed var(--gapline)}
details{margin-top:8px;font-size:13px}summary{cursor:pointer;color:var(--ink2)}
table{border-collapse:collapse;width:100%;font-size:12.5px;margin-top:8px}
th,td{padding:5px 8px;text-align:right;border-bottom:1px solid var(--grid);font-variant-numeric:tabular-nums}
th:first-child,td:first-child{text-align:left}thead th{color:var(--ink2);font-weight:500}
.table-scroll,section:has(table){overflow-x:auto}
td.hm{background:color-mix(in srgb,var(--s1) calc(var(--p)*0.55%),transparent)}
td.na{color:var(--dim);text-align:center}
tr.key th,tr.key td{font-weight:600}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:12px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px 14px;min-width:0}
.card h3{font-size:14px;margin:4px 0 6px;font-weight:600}.card p{margin:0 0 6px;font-size:13px;color:var(--ink2)}
.card .src{font:11.5px var(--mono);color:var(--dim);overflow-wrap:anywhere}
.tag{display:inline-block;font-size:11px;letter-spacing:.06em;text-transform:uppercase;font-weight:600;padding:1px 8px;border-radius:99px;border:1px solid}
.card.ex .tag{color:var(--ok);border-color:var(--ok)}.card.out .tag{color:var(--bad);border-color:var(--bad)}
.foot{color:var(--dim);font-size:12px;margin-top:10px}
@media (max-width:520px){.tile .v{font-size:19px}}
'''
JS = '''
(function(){var tabs=document.querySelectorAll('nav button'),pages=document.querySelectorAll('.tab');
function show(id){var f=false;tabs.forEach(function(b){var on=b.dataset.t===id;b.setAttribute('aria-selected',on);if(on)f=true});if(!f)id=tabs[0].dataset.t;
pages.forEach(function(p){p.classList.toggle('on',p.id==='tab-'+id)});tabs.forEach(function(b){b.setAttribute('aria-selected',b.dataset.t===id)})}
tabs.forEach(function(b){b.addEventListener('click',function(){show(b.dataset.t);try{history.replaceState(null,'','#'+b.dataset.t)}catch(e){}})});
show((location.hash||'').slice(1))})();
'''
page = f'''<title>MATPredict Campaign Dashboard</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header><h1>MATPredict campaign dashboard</h1><p>Whole-clade runs on the BFD genomes with code v0.6.0 (2026-10-03 to 05): how many genomes were called, what the calls contain, and where the locus is arranged differently.</p></header>
<nav role="tablist"><button role="tab" data-t="summary">Summary</button><button role="tab" data-t="exemplars">Exemplars and outliers</button><button role="tab" data-t="synteny">Xylariales synteny</button><button role="tab" data-t="size">Locus size and content</button></nav>
<div id="tab-summary" class="tab">{tab_summary}</div>
<div id="tab-exemplars" class="tab">{tab_ex}</div>
<div id="tab-synteny" class="tab">{tab_syn}</div>
<div id="tab-size" class="tab">{tab_size}</div>
<p class="foot">Built by results/2026-10-05_campaign_overview/make_dashboard.py from the committed tables. "Called" means at least one MAT call; most clades have no curated reference, so uncalled does not mean absent.</p>
</div>
<script>{JS}</script>
'''
open(os.path.join(HERE, 'dashboard.html'), 'w').write(page)
print('dashboard.html', len(page), 'bytes; species', n_sp, 'eutypa', eutypa, 'maj', dict(sp_maj))

# ------------------------------------------------------------------ print version (all tabs, with a reading guide) for the PDF
GUIDE = '''
<section class="guide"><h2>How to read this report</h2>
<p>MATPredict finds mating-type (MAT) loci in fungal genomes. This report summarises three whole-clade runs of code version 0.6.0 over the BFD genome collection
(Ascomycota, Basidiomycota and the Mucoromycotina part of Mucoromycota), then looks at two kinds of change: gene order near the MAT region in Xylariales, and locus
length and gene content across groups. Everything is drawn from tables committed in the repository; nothing here is new computation.</p>
<dl>
<dt>Called</dt><dd>A genome with at least one MAT call. This measures detection, not truth. Most clades have no curated reference locus, so an uncalled genome is a
reference gap until shown otherwise.</dd>
<dt>Routed by lineage / phylum fallback</dt><dd>Genomes whose order has its own curated record are searched with it (blue). Others are searched against the whole phylum
(orange) and call less reliably.</dd>
<dt>Locus class</dt><dd>Full locus: the MAT gene(s) plus the expected flank structure. Partial locus: some of it. Gene only: the idiomorph gene without a recognised locus.
Homothallic candidate: both idiomorphs at one locus.</dd>
<dt>Confidence</dt><dd>High, medium or low, from how complete the evidence is.</dd>
<dt>PR-only calls</dt><dd>Basidiomycota genomes whose only call is the pheromone-receptor family. 1,360 of these come from a motif scan alone and are labelled unverified.</dd>
<dt>SAC and SCA (Xylariales)</dt><dd>Orders of the flank genes SLA2, COX13 and APN2. SAC keeps the outgroup order with an open interval between SLA2 and APN2 where a MAT locus
sits. SCA has COX13 between SLA2 and APN2, so the interval is closed.</dd>
</dl>
<h3>Panel by panel</h3>
<ul>
<li><b>Summary.</b> Headline counts; call rate per class or order (a marker shows the rate without PR-only calls); and the composition of the calls. Takeaway: rates are high
where a clade has its own record and low where it does not.</li>
<li><b>Exemplars and outliers.</b> Short cards with a source file each. Exemplars show what the approach recovers; outliers show where it falls short or where biology differs.</li>
<li><b>Xylariales synteny.</b> A gene-order map of seven genomes on one frame; the sampling bar chart (genomes versus species); a table of which neighbours remain adjacent; and
a histogram of the SLA2 to APN2 distance. Takeaway: the NC1011 order is shared by many species and its MAT interval is closed.</li>
<li><b>Locus size and content.</b> Locus length per Mucoromycotina genus and per filamentous Ascomycota class (bar = interquartile range, dot = median, line = range), a table
of which genes each Mucoromycotina genus carries, and the share of Ascomycota loci carrying both flank genes.</li>
</ul>
<p class="note">Limits: counts are descriptive; several runs used earlier code; genome counts reflect uneven strain sampling; synteny is gene-level only. The markdown overview is
<code>analysis/2026-10-05_campaign-overview.md</code> and the guide is <code>analysis/2026-10-05_campaign-dashboard-guide.md</code>.</p></section>'''
PRINT_CSS = '''
@page{size:Letter;margin:13mm 12mm}
body{background:#fff;color:#101312;font-size:11.5px}
:root{--bg:#fff;--panel:#fff;--grid:#e3e5e2;--line:#cfd2cf}
.wrap{max-width:none;padding:0}
nav{display:none}
.tab{display:block !important;break-before:page}
.tab:first-of-type{break-before:auto}
.tab>h1{font-size:18px;margin:0 0 8px}
section,.card,.tile{break-inside:avoid}
section{border-color:#cfd2cf}
.chart.wide{min-width:0}
details{display:block}details>summary{display:none}
.guide dl{display:grid;grid-template-columns:170px 1fr;gap:4px 12px;margin:8px 0}.guide dt{font-weight:600}.guide dd{margin:0;color:var(--ink2)}
.guide h3{font-size:13px;margin:12px 0 4px}.guide li{margin:0 0 4px}
.cards{grid-template-columns:repeat(2,1fr)}
'''
def tabbed(title, body): return f'<div class="tab"><h1>{E(title)}</h1>{body}</div>'
print_page = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><title>MATPredict campaign report</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>{CSS}{PRINT_CSS}</style></head><body><div class="wrap">
<header><h1>MATPredict campaign report</h1><p>Whole-clade runs on the BFD genomes with code v0.6.0 (2026-10-03 to 05). Generated {__import__("datetime").date.today().isoformat()} from committed tables by results/2026-10-05_campaign_overview/make_dashboard.py.</p></header>
{GUIDE}
{tabbed('1. Summary', tab_summary)}{tabbed('2. Exemplars and outliers', tab_ex)}{tabbed('3. Xylariales synteny', tab_syn)}{tabbed('4. Locus size and content', tab_size)}
</div></body></html>'''
open(os.path.join(HERE, 'campaign_report_print.html'), 'w').write(print_page)
print('print version written')
