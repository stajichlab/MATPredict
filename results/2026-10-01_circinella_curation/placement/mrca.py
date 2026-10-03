"""Is a target inside the MRCA clade of all known Plus (resp. Minus) tips, and does
that clade exclude the opposite type? Report MRCA size, opposite-type tips inside,
support. Also the same using ONLY tier-1/established curated refs (excluding the
new tier-2 genome-derived records 13706, 41833, 44442 and labelled strains)."""
import sys
from Bio import Phylo
sys.path.insert(0, '.')
from place import ktype
TIER2 = ('13706_', '41833_', '44442_')
def run(tf, targets):
    t = Phylo.read(tf, 'newick')
    for x in t.get_terminals(): x.name = x.name.replace('|', '__')
    out = [x for x in t.get_terminals() if x.name.startswith('OUT__MAT1-2-1__5518')]
    t.root_with_outgroup(out[0])
    for mode in ('all_known', 'established_refs_only'):
        for T, opp in (('Plus', 'Minus'), ('Minus', 'Plus')):
            def ok(n):
                k = ktype(n)
                if mode == 'established_refs_only':
                    return k == T and n.startswith('REF__') and not any(z in n for z in TIER2)
                return k == T
            tips = [x for x in t.get_terminals() if ok(x.name)]
            m = t.common_ancestor(tips)
            names = {x.name for x in m.get_terminals()}
            nopp = sum(1 for n in names if ktype(n) == opp)
            sup = m.name or (m.confidence if m.confidence is not None else '')
            inside = {tg: any(tg in n for n in names) for tg in targets}
            print(f"  {mode:22s} {T:5s} MRCA {len(names):3d} tips, opposite-type inside {nopp:2d}, support {sup}; " + ', '.join(f"{k} inside={v}" for k, v in inside.items()))
if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2:])
