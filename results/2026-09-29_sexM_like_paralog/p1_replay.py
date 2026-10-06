"""Replay a third classifier class for the P1 sexM-like paralog: HMM built from the
single BFD (non-held-out) P1 copy GCA_000697295.1 (Mucor indicus B7402); a called core
protein is flagged when its P1 score exceeds both its sexM and sexP scores."""
import csv, collections, pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
from Bio import SeqIO
S = lambda x: x.decode() if isinstance(x, bytes) else x
abc = Alphabet.amino()
train = [(r.id, str(r.seq)) for r in SeqIO.parse('P1_train_bfd.faa', 'fasta') if r.id.startswith('GCA_000697295')]
seq = TextSequence(name=train[0][0].encode(), sequence=train[0][1]).digitize(abc)
p1, _, _ = pyhmmer.plan7.Builder(abc).build(seq, pyhmmer.plan7.Background(abc))
p1.name = b'P1_paralog'
with pyhmmer.plan7.HMMFile('clf_sexM.hmm') as f: hM = f.read()
with pyhmmer.plan7.HMMFile('clf_sexP.hmm') as f: hP = f.read()
hmms = [hM, hP, p1]
names = [S(h.name) for h in hmms]
m = list(csv.DictReader(open('models.tsv'), delimiter='\t'))
by = collections.defaultdict(list)
for x in m:
    if x['gene'] in ('sexM', 'sexP'): by[(x['org'], x['call'])].append(x)
core = {}
for k, v in by.items():
    want = {'Minus': 'sexM', 'Plus': 'sexP'}.get(v[0]['idiomorph'])
    c = [x for x in v if x['gene'] == want]
    if c:
        core[f"LCG|{k[0]}|{k[1]}|{v[0]['idiomorph']}|{v[0]['confidence']}|{v[0]['locus_class']}"] = max(c, key=lambda x: float(x['bitscore'] or 0))['protein']
for r in SeqIO.parse('jena_weak.faa', 'fasta'): core[f"JENA|{r.id}|0|Minus|weak|"] = str(r.seq)
block = DigitalSequenceBlock(abc, [TextSequence(name=k.encode(), sequence=v).digitize(abc) for k, v in core.items() if v])
res = collections.defaultdict(dict)
for hits in pyhmmer.hmmsearch(hmms, block, E=1e3):
    q = S(hits.query.name)
    for h in hits:
        res[S(h.name)][q] = max(res[S(h.name)].get(q, -1e9), h.score)
rows = []
for k in core:
    s = [res[k].get(n, 0.0) for n in names]
    rows.append((k, *s, s[2] > max(s[0], s[1])))
with open('p1_replay.tsv', 'w') as fo:
    fo.write('call\tsexM\tsexP\tP1\tflagged\n')
    for r in rows: fo.write(f"{r[0]}\t{r[1]:.1f}\t{r[2]:.1f}\t{r[3]:.1f}\t{r[4]}\n")
fl = [r for r in rows if r[4]]
print('scored', len(rows), 'flagged', len(fl))
for r in fl: print(' ', r[0][:72].ljust(73), 'M', round(r[1], 1), 'P', round(r[2], 1), 'P1', round(r[3], 1))
