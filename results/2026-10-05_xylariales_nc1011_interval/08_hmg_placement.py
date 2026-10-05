"""Place Xylariales HMG-box proteins in the RAxML-NG tree relative to MAT1-2-1/MAT1-1-3 references, NCU03481 and fmf-1.

Class of a labelled reference: MAT_ref (database MAT1-2-1/MAT1-1-3-type proteins), NCU03481, fmf1.
Writes tree/hmg_placement.tsv: nearest labelled reference (patristic distance), the three nearest, and whether the
smallest clade (midpoint-rooted) containing the protein and at least one labelled reference is single-class.
"""
import re, sys, collections
from Bio import Phylo
D = sys.argv[1]; T = f'{D}/tree'
hdr = {}
for l in list(open(f'{T}/hmg_domains.tsv'))[1:]:
    f = l.rstrip('\n').split('\t'); hdr[f[0]] = f[5]
def cls(name):
    h = hdr.get(name, '')
    if name.startswith('DBREF|') and re.search(r'(MAT1-2-1|mat a-1|mt a-1|MAT1-1-3|matA-3|mat A-3)', name, re.I): return 'MAT_ref'
    if 'NCU03481' in h or 'NCU03481' in name: return 'NCU03481'
    if 'female and male fertility' in h or 'NCU09387' in h: return 'fmf1'
    return None
tree = Phylo.read(f'{T}/rx.raxml.support', 'newick')
tree.root_at_midpoint()
tips = {t.name: t for t in tree.get_terminals()}
lab = {n: cls(n) for n in tips if cls(n)}
print('labelled references:', collections.Counter(lab.values()))
queries = [n for n in tips if n.startswith(('XYL_ours|', 'XYL_paper|'))]
rows = []
for q in queries:
    d = sorted((tree.distance(tips[q], tips[r]), r) for r in lab)
    near = d[0]
    top3 = ';'.join(f'{lab[r]}:{x:.2f}' for x, r in d[:3])
    # smallest clade containing q and >= 1 labelled reference
    path = tree.get_path(tips[q])
    clade = None
    for node in reversed([tree.root] + path[:-1]):
        names = {t.name for t in node.get_terminals()}
        if names & set(lab): clade = node; break
    cl = collections.Counter(lab[n] for n in {t.name for t in clade.get_terminals()} if n in lab)
    ntip = len(clade.get_terminals())
    rows.append((q, hdr.get(q, '')[:60], f'{near[0]:.3f}', lab[near[1]], top3, ';'.join(f'{k}={v}' for k, v in cl.items()), ntip, clade.confidence if clade.confidence is not None else ''))
with open(f'{T}/hmg_placement.tsv', 'w') as o:
    o.write('protein\theader\tnearest_ref_dist\tnearest_ref_class\tthree_nearest\tref_classes_in_smallest_clade\tclade_tips\tclade_support\n')
    for r in rows: o.write('\t'.join(str(x) for x in r) + '\n')
c = collections.Counter(r[3] for r in rows); print('nearest labelled class among Xylariales HMG proteins:', dict(c))
# is the MAT_ref set monophyletic (midpoint-rooted)?
mat = [n for n in lab if lab[n] == 'MAT_ref']
if len(mat) > 1:
    mrca = tree.common_ancestor([tips[n] for n in mat]); inside = {t.name for t in mrca.get_terminals()}
    print('MAT_ref MRCA clade: tips', len(inside), 'support', mrca.confidence, 'other labelled classes inside:', collections.Counter(lab[n] for n in inside if n in lab and lab[n] != 'MAT_ref'), 'Xylariales inside:', sum(n in set(queries) for n in inside))
