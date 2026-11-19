"""Top phmmer hits of chosen S. schoenii ORFs against refs.faa. Usage: orf_hits.py data/longest_orfs.pep refs.faa ORF..."""
import sys, pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet, TextSequence
pep, refs, want = sys.argv[1], sys.argv[2], set(sys.argv[3:])
alpha = Alphabet.amino(); seqs = []; name = None; cur = {}
for l in open(pep):
    if l.startswith(">"): name = l[1:].split()[0]; cur[name] = ""
    elif name: cur[name] += l.strip()
q = [TextSequence(name=n.encode(), sequence=cur[n].rstrip("*")).digitize(alpha) for n in want]
with SequenceFile(refs, digital=True, alphabet=alpha) as sf: db = sf.read_block()
for hits in pyhmmer.phmmer(q, db, E=1):
    print(hits.query.name, len(cur[hits.query.name.decode() if isinstance(hits.query.name, bytes) else hits.query.name]) if False else "")
    for h in list(hits)[:7]:
        n = h.name.decode() if isinstance(h.name, bytes) else h.name
        print("   %.0f  %s" % (h.score, n[:90]))
