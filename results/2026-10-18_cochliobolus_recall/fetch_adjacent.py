"""Fetch the three Bipolaris maydis proteins adjacent to the MAT idiomorph (C4 deposit AF027687.1: GAP1 AAB84002.1, ORF1 AAB84003.1, BGL1 AAB84005.1)."""
import urllib.request
names = {'AAB84002.1': 'GAP1', 'AAB84003.1': 'ORF1', 'AAB84005.1': 'BGL1'}
url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=' + ','.join(names) + '&rettype=fasta&retmode=text'
txt = urllib.request.urlopen(url).read().decode()
out = []
for rec in txt.split('>')[1:]:
    head, *seq = rec.strip().split('\n')
    acc = head.split()[0]
    out.append(f'>{acc.split(".")[0]}|gene0|{names[acc]}\n' + ''.join(seq))
open('results/2026-10-18_cochliobolus_recall/mat_adjacent_refs.faa', 'w').write('\n'.join(out) + '\n')
print(len(out), 'proteins')
