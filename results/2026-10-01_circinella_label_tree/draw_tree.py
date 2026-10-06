"""Named PDF of the Circinella-label tree.

Tips: species + strain + RSA/LCG label. Colour = label (Plus orange, Minus blue,
none grey); curated/Zygo references bold, coloured by type; outgroups purple.
Right strip = clade placement (sexP / sexM) from clades.py on the same tree.
Rooted on Fusarium MAT1-2-1. Node dots: filled UFBoot >= 95 (or FastTree >= 0.95),
open 80-94.
Usage: draw_tree.py TREE CLADES_TSV OUT.pdf "title"
"""
import csv, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo
tf, ctsv, outpdf, title = sys.argv[1:5]
tips = {r['tip']: r for r in csv.DictReader(open('tips.tsv'), delimiter='\t')}
GN = {r['gid']: r['notes'] for r in csv.DictReader(open('genomes.tsv'), delimiter='\t')}
for r in tips.values():
    if 'misidentified' in GN.get(r['gid'], ''):
        r['note'] = 'misidentified: ' + GN[r['gid']].split('prior finding: ')[1]
place = {r['tip']: r['placement'] for r in csv.DictReader((l for l in open(ctsv) if not l.startswith('#')), delimiter='\t')}
t = Phylo.read(tf, 'newick')
for c in t.get_terminals():
    c.name = c.name.replace('__', '|')
def support(cl):
    """UFBoot (IQ-TREE 'SH/UFB'), RAxML/FastTree confidence; returns 0-100 or None."""
    v = None
    if cl.name and '/' in cl.name:
        v = float(cl.name.split('/')[-1])
    elif cl.confidence is not None:
        v = cl.confidence * 100 if cl.confidence <= 1 else cl.confidence
    elif cl.name:
        try: v = float(cl.name); v = v * 100 if v <= 1 else v
        except ValueError: pass
    return v
root = [c for c in t.get_terminals() if c.name.startswith('OUT|MAT1-2-1|5518')][0]
t.root_with_outgroup(root)
t.ladderize(reverse=True)
def label(r):
    g = r['gid']
    if r['group'] == 'outgroup':
        return f"{r['species']}  [outgroup]"
    if r['group'] == 'curated_ref':
        return f"REF {g}  [{r['label']}]"
    if r['group'] == 'zygo_truth':
        return f"{g.replace('_', ' ')}  [Zygo truth {r['label']}]"
    sp = (r['species'] or g).split(' (')[0]
    strain = g.replace(sp.replace(' ', '_'), '').strip('_').replace('_', ' ')
    lab = f"  [label {r['label']}]" if r['label'] else ''
    note = ''
    if 'misidentified' in r['note']:
        note = '  (misidentified' + ('; = ' + r['note'].split(': ')[1] if ': ' in r['note'] else '') + ')'
    elif r['source'] == 'fallback': note = f"  ({r['note']})"
    return f"{sp} {strain}{lab}{note}".replace('  ', ' ', 0)
COL = {'Plus': '#e6550d', 'Minus': '#2b8cbe', '': '#555555'}
REFCOL = {'sexP': '#c41e1e', 'Plus': '#c41e1e', 'sexM': '#1f4fbf', 'Minus': '#1f4fbf'}
terms = t.get_terminals(); n = len(terms)
ypos = {c: i for i, c in enumerate(terms)}
depth = t.depths()
def y_of(cl):
    if cl in ypos: return ypos[cl]
    ys = [y_of(c) for c in cl.clades]; ypos[cl] = (min(ys) + max(ys)) / 2; return ypos[cl]
y_of(t.root)
maxx = max(depth.values())
fig, ax = plt.subplots(figsize=(13, max(8, n * 0.14 + 2)))
for cl in t.find_clades():
    x, y = depth[cl], ypos[cl]
    if cl.clades:
        ys = [ypos[c] for c in cl.clades]
        ax.plot([x, x], [min(ys), max(ys)], color='black', lw=0.5)
        for c in cl.clades:
            ax.plot([x, depth[c]], [ypos[c], ypos[c]], color='black', lw=0.5)
        s = support(cl)
        if s is not None and cl is not t.root:
            if s >= 95: ax.plot(x, y, 'o', ms=2.4, color='black')
            elif s >= 80: ax.plot(x, y, 'o', ms=2.4, mfc='white', mec='black', mew=0.4)
            if s >= 50 and len(cl.get_terminals()) >= 4:
                ax.text(x - maxx * 0.004, y - 0.35, f'{s:.0f}', fontsize=3.5, ha='right', va='bottom', color='#444444')
for c in terms:
    r = tips[c.name]
    ref = r['group'] in ('curated_ref', 'zygo_truth')
    col = '#6a3d9a' if r['group'] == 'outgroup' else REFCOL.get(r['label']) if ref else COL.get(r['label'], '#555555')
    ax.text(depth[c] + maxx * 0.006, ypos[c], label(r), va='center', ha='left', fontsize=5.3, color=col,
            fontweight='bold' if ref else 'normal')
sx = maxx * 1.55
for c in terms:
    r = tips[c.name]
    p = place.get(c.name) or ({'sexP': 'sexP', 'Plus': 'sexP', 'sexM': 'sexM', 'Minus': 'sexM'}.get(r['label']) if r['group'] in ('curated_ref', 'zygo_truth') else '')
    if p in ('sexP', 'sexM'):
        ax.plot([sx, sx], [ypos[c] - 0.5, ypos[c] + 0.5], color=REFCOL[p], lw=5, solid_capstyle='butt')
ax.text(sx, -1.5, 'clade', ha='center', fontsize=6)
ax.set_ylim(n, -1); ax.set_xlim(0, maxx * 1.58); ax.axis('off')
ax.plot([0, 0.5], [n - 0.2, n - 0.2], color='black', lw=0.8)
ax.text(0.25, n + 0.9, '0.5 substitutions/site', ha='center', fontsize=6)
h = [Line2D([], [], color=c, lw=0, marker='s', ms=6, label=l) for c, l in [
    ('#c41e1e', 'curated sexP reference / Zygo truth Plus (bold)'), ('#1f4fbf', 'curated sexM reference / Zygo truth Minus (bold)'),
    ('#e6550d', 'strain label Plus'), ('#2b8cbe', 'strain label Minus'), ('#555555', 'no label'),
    ('#6a3d9a', 'outgroup (root: Fusarium MAT1-2-1)')]]
h += [Line2D([], [], color='#c41e1e', lw=5, label='strip: placed with sexP references only'),
      Line2D([], [], color='#1f4fbf', lw=5, label='strip: placed with sexM references only'),
      Line2D([], [], color='black', lw=0, marker='o', ms=4, label='support >= 95'),
      Line2D([], [], mfc='white', mec='black', lw=0, marker='o', ms=4, label='support 80-94')]
ax.legend(handles=h, loc='upper left', fontsize=7, frameon=False)
ax.set_title(title, fontsize=9, loc='left')
fig.savefig(outpdf, bbox_inches='tight')
print('wrote', outpdf, n, 'tips')
