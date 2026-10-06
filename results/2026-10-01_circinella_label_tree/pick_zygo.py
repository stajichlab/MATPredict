"""Pick the best-typed protein per Zygo truth locus (score in the truth direction)."""
import sys
sys.argv=[sys.argv[0]]
exec(open('extract.py').read().split('rows = list(')[0])  # reuse score()
seqs, nm = {}, None
for l in open(str(R / '2026-09-26_sexMP_hmm/zygo_locus_proteins.faa')):
    if l.startswith('>'): nm = l[1:].strip(); seqs[nm] = ''
    else: seqs[nm] += l.strip()
best = {}
for h, s in seqs.items():
    _, org, lab, pid = h.split('|')
    sc = score(s)
    k = 'sexP' if lab == 'Plus' else 'sexM'
    v = sc[k] - max(sc['sexM' if k == 'sexP' else 'sexP'], sc['P1'])
    if sc[k] >= 60 and (org not in best or v > best[org][0]):
        best[org] = (v, lab, pid, sc, s.replace('*', ''))
with open('zygo_best.faa', 'w') as fa, open('zygo_best.tsv', 'w') as t:
    t.write('org\tlabel\tprotein\tsexM\tsexP\tP1\tmargin\tlen\n')
    for org, (v, lab, pid, sc, s) in sorted(best.items()):
        fa.write(f'>ZYGO|{org}|{lab}|{pid}\n{s}\n')
        t.write(f"{org}\t{lab}\t{pid}\t{sc['sexM']}\t{sc['sexP']}\t{sc['P1']}\t{v:.1f}\t{len(s)}\n")
print(len(best))
