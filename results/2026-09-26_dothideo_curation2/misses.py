"""Explain every genome with no Ascomycota:MAT call in a run: routing, what was
withheld and why, and assembly contiguity (from the BFD library genome)."""
import csv, gzip, sys
from pathlib import Path
import yaml

HERE = Path(__file__).parent
TAG = sys.argv[1]
PREV = Path('/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_dothideo_curation/dca6ccc/runs')
LIB = Path('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes')
S = {r['ASMID']: r for r in csv.DictReader(open('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv'))}

def called(p):
    if not p.exists():
        return None
    r = yaml.safe_load(p.read_text())
    return r, [x for x in (r.get('detected') or []) if x['family'] == 'Ascomycota:MAT']

def asm_stats(asmid):
    lens, n, cur = [], 0, 0
    with gzip.open(LIB / f'{asmid}.fa.gz', 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if cur: lens.append(cur)
                cur = 0
            else:
                l = line.strip(); cur += len(l); n += l.upper().count('N')
    if cur: lens.append(cur)
    lens.sort(reverse=True); tot = sum(lens); acc = 0; n50 = 0
    for L in lens:
        acc += L
        if acc >= tot / 2: n50 = L; break
    return len(lens), tot, n50, n

prev_uncalled = sorted(d.name for d in PREV.iterdir() if not (called(d / 'detection_report.yaml') or (None, []))[1])
print(f'previous pilot (dca6ccc) uncalled: {len(prev_uncalled)}')
for g in prev_uncalled:
    res = called(HERE / TAG / 'runs' / g / 'detection_report.yaml')
    s = S.get(g, {})
    ctg, tot, n50, nn = asm_stats(g)
    head = f"{g}  {s.get('SPECIES','?')}  {s.get('ORDER','?')}/{s.get('FAMILY','?')}  contigs={ctg} size={tot/1e6:.1f}Mb N50={n50/1e3:.1f}kb N={nn}"
    if res is None:
        print(head, '| NO REPORT (suppressed or failed)'); continue
    r, calls = res
    print(head)
    print(f"   now: {'CALLED ' + ','.join(x['idiomorph'] + '/' + x['confidence'] for x in calls) if calls else 'uncalled'}  routing={r.get('routing_mode')}")
    for nd in r.get('not_detected') or []:
        if nd.get('family') == 'Ascomycota:MAT':
            print(f"   not_detected: {nd.get('reason')} | found={nd.get('genes_found')} best_fraction={nd.get('best_fraction_found')}")
    for sl in r.get('suppressed_loci') or []:
        if sl.get('family') == 'Ascomycota:MAT':
            print(f"   withheld: {sl['contig']}:{sl['start']}-{sl['end']} idio={sl.get('idiomorph')} polished={sl.get('polished_genes')} genes={sl.get('genes_found')} reason={sl.get('reason','')}")
