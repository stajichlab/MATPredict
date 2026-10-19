"""Compare a run's calls (summarize_runs.py TSV) with the independent six-frame scan truth (clade_gates.tsv: alpha-box ORF -> MAT1-1,
HMG ORF -> MAT1-2). usage: eval_calls.py runs.tsv label [out.tsv]"""
import csv, re, sys
from collections import Counter, defaultdict
runs = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
label = sys.argv[2]
gates = list(csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_gates.tsv'), delimiter='\t'))
exp = defaultdict(dict)
meta = {}
for g in gates:
    idm = 'MAT1-1' if g['gene'] == 'MAT1-1-1' else 'MAT1-2'
    m = re.match(r'(.+):(\d+)-(\d+)', g['scan_where'])
    exp[g['genome']][idm] = (m.group(1), int(m.group(2)), int(m.group(3)))
    meta[g['genome']] = (g['genus'], g['species'], g['clade'])
calls = defaultdict(list)
for r in runs:
    if r['n_detected'] != '0':
        calls[r['genome']].append(r)
allg = {r['genome'] for r in runs}
res = Counter(); rows = []
for g in sorted(allg):
    e = exp.get(g, {})
    c = calls.get(g, [])
    got = {}
    on_orf = {}
    for r in c:
        got[r['idiomorph']] = r
        if r['idiomorph'] in e:
            ctg, s, en = e[r['idiomorph']]
            on_orf[r['idiomorph']] = r['contig'] == ctg and int(r['end']) >= s - 500 and int(r['start']) <= en + 500
    gen = meta.get(g, ('?', '', ''))[0]
    if not c:
        out = 'no_call'
    else:
        right = [i for i in got if i in e and on_orf.get(i)]
        wrong = [i for i in got if i not in right]
        if wrong: out = 'WRONG_or_unsupported:' + ','.join(wrong)
        elif set(right) == set(e): out = 'all_expected_idiomorphs'
        else: out = 'subset_of_expected'
    res[(gen, out)] += 1
    rows.append(dict(genome=g, genus=gen, expected=','.join(sorted(e)), called=','.join(sorted(got)), outcome=out,
                     conf=','.join(r['confidence'] for r in c), cls=','.join(r['locus_class'] for r in c), genes=';'.join(r['genes_found'] for r in c)))
for k in sorted(res): print(label, k, res[k])
print(label, 'total genomes', len(allg), 'called', len(calls), 'expected-idiomorph-correct (any)', sum(1 for r in rows if r['outcome'] in ('all_expected_idiomorphs', 'subset_of_expected')))
if len(sys.argv) > 3:
    with open(sys.argv[3], 'w') as f:
        w = csv.DictWriter(f, list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
