"""Root on Fusarium MAT1-2-1; report placement of each tip relative to the curated sexP / sexM references.

For every tip: the smallest clade containing it plus >=1 curated/Zygo reference;
the reference types in that clade; its support. Also monophyly of reference sets.
Usage: clades.py tree.newick > out.tsv
"""
import sys, csv
from Bio import Phylo
t = Phylo.read(sys.argv[1], 'newick')
for c in t.get_terminals():
    c.name = c.name.replace('__', '|')
root = [c for c in t.get_terminals() if c.name.startswith('OUT|MAT1-2-1|5518')][0]
t.root_with_outgroup(root)
tips = {r['tip']: r for r in csv.DictReader(open('tips.tsv'), delimiter='\t')}
def rtype(n):
    r = tips.get(n)
    if not r: return None
    if r['group'] in ('curated_ref', 'zygo_truth'):
        return {'sexM': 'M', 'Minus': 'M', 'sexP': 'P', 'Plus': 'P'}.get(r['label'])
    return None
def sup(c):
    if c.confidence is not None: return c.confidence
    if c.name and '/' in str(c.name): return c.name
    return ''
P = [n.name for n in t.get_terminals() if rtype(n.name) == 'P']
M = [n.name for n in t.get_terminals() if rtype(n.name) == 'M']
half = {}
for lab, s in (('sexP_refs', P), ('sexM_refs', M)):
    ca = t.common_ancestor(s)
    half[lab] = (set(x.name for x in ca.get_terminals()), sup(ca))
    inside = [x.name for x in ca.get_terminals()]
    others = [x for x in inside if x not in s]
    print(f'# {lab}: {len(s)} refs; MRCA holds {len(inside)} tips; non-ref tips inside: {len(others)}; support {sup(ca)}')
    print('#   other ref-type inside:', [x for x in others if rtype(x)])
print('tip\tlabel\tclade_refs_P\tclade_refs_M\tclade_size\tsupport\tplacement\tin_sexP_ref_mrca\tin_sexM_ref_mrca')
for n in t.get_terminals():
    if rtype(n.name) or n.name.startswith('OUT'):
        continue
    path = t.get_path(n)
    for anc in reversed([t.root] + path[:-1]):
        refs = [rtype(x.name) for x in anc.get_terminals() if rtype(x.name)]
        if refs:
            p, m = refs.count('P'), refs.count('M')
            pl = 'sexP' if p and not m else 'sexM' if m and not p else 'mixed'
            print(f"{n.name}\t{tips[n.name]['label'] or '.'}\t{p}\t{m}\t{len(anc.get_terminals())}\t{sup(anc)}\t{pl}\t{half['sexP_refs'][1] if n.name in half['sexP_refs'][0] else 'no'}\t{half['sexM_refs'][1] if n.name in half['sexM_refs'][0] else 'no'}")
            break
