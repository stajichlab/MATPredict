"""Size and support of the smallest clade holding a named set of tips, and of
the Phascolomyces + Absidia NRRL 3163 pair. Rooted on Fusarium MAT1-2-1."""
import sys
from Bio import Phylo
def run(tf):
    t = Phylo.read(tf, 'newick')
    for x in t.get_terminals(): x.name = x.name.replace('|', '__')
    t.root_with_outgroup([x for x in t.get_terminals() if x.name.startswith('OUT__MAT1-2-1__5518')][0])
    find = lambda s: [x for x in t.get_terminals() if s in x.name]
    for label, keys in (('Phaart1+Abs3163', ['Phaart1', 'Absidia_sp']),
                        ('Phaart1+Abs3163+Cirumb1', ['Phaart1', 'Absidia_sp', 'Cirumb1']),
                        ('Circinella sexP-type tips (CG__Circinella_* in sexP half)', ['Cirumb1', 'Circinella_muscae_NRRL_1357', 'Circinella_angarensis_NRRL_1594'])):
        tips = [y for k in keys for y in find(k)]
        m = t.common_ancestor(tips)
        print(f"  {label}: MRCA {len(m.get_terminals())} tips, support {m.name or m.confidence}")
for f in sys.argv[1:]:
    print(f); run(f)
