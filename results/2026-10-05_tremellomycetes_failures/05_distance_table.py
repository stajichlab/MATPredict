import csv, statistics as st
R = list(csv.DictReader(open('joined.tsv'), delimiter='\t'))
def q(o, c, good='1'):
    return [r for r in R if r['order'] == o and r['cls'] == c and r['good_asm'] == good]
for o in ['Tremellales', 'Trichosporonales', 'Filobasidiales', 'Cystofilobasidiales']:
    print(o, {c: len(q(o, c)) for c in ['called_other', 'PR_only', 'uncalled']}, 'poor', {c: len(q(o, c, '0')) for c in ['called_other', 'PR_only', 'uncalled']})
    for c in ['called_other', 'PR_only', 'uncalled']:
        x = q(o, c)
        d = [int(r['d_kn_ste3']) for r in x if r['d_kn_ste3'] != '']
        if d:
            print('  ', c, 'KN-STE3 same contig', len(d), '/', len(x), 'median', st.median(d), '<=150kb', sum(v <= 150000 for v in d))
