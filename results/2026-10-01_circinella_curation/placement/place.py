"""Placement of target tips: for each tip, walk up toward the root and report
the first clade containing known-type tips (curated refs, Zygo truth, labelled
references), with its composition and support. Trees rooted on the Fusarium
MAT1-2-1 outgroup."""
import csv, sys
from Bio import Phylo
T = '/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_circinella_label_tree'
lab = {}
for r in csv.DictReader(open(f'{T}/tips.tsv'), delimiter='\t'):
    lab[r['tip'].replace('|', '__')] = (r['group'], r['label'])
def ktype(name):
    if name.startswith('REF__'):
        return 'Plus' if name.endswith('sexP') else 'Minus'
    g, l = lab.get(name, ('', ''))
    if name.startswith(('ZYGO__', 'LAB__')) and l in ('Plus', 'Minus'):
        return l
    return None
def support(cl):
    s = cl.name if cl.name else (str(cl.confidence) if cl.confidence is not None else '')
    return s
def place(treefile, targets, extra_known=None):
    t = Phylo.read(treefile, 'newick')
    for x in t.get_terminals(): x.name = x.name.replace('|', '__')
    out = [x for x in t.get_terminals() if x.name.startswith('OUT__MAT1-2-1__5518')]
    if out: t.root_with_outgroup(out[0])
    res = {}
    for tg in targets:
        leaf = [x for x in t.get_terminals() if tg in x.name]
        if not leaf: res[tg] = 'absent'; continue
        path = t.get_path(leaf[0])
        for anc in reversed(path[:-1]):
            kinds = [ktype(x.name) for x in anc.get_terminals()]
            k = [x for x in kinds if x]
            if k:
                p, m = k.count('Plus'), k.count('Minus')
                names=[x.name.split("__")[1][:28] for x in anc.get_terminals() if ktype(x.name)]
                res[tg] = f"first known-type clade: {len(anc.get_terminals())} tips, Plus {p} / Minus {m}, support {support(anc)}; known: {names}"
                break
        else:
            res[tg] = 'no known-type clade'
    return res
if __name__ == '__main__':
    tree = sys.argv[1]; targets = sys.argv[2:]
    for k, v in place(tree, targets).items(): print(f'{k}\t{v}')
