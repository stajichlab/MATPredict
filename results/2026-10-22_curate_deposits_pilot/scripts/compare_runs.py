"""Compare base vs pilot detection reports genome by genome; print every change."""
import sys, yaml, csv, json, collections
from pathlib import Path
P = Path(sys.argv[1]); out = Path(sys.argv[2])
meta = {r['asm']: r for r in csv.DictReader(open(P / 'replay_genomes.tsv'), delimiter='\t')}
def load(sub, asm):
    f = P / 'runs' / sub / asm / 'detection_report.yaml'
    return yaml.safe_load(f.read_text()) if f.exists() else None
def core(l):
    return {g['gene']: g for g in l.get('gene_evidence') or [] if g.get('role') == 'core_MAT'}
def key(l): return (l['family'], l['contig'])
def summ(l):
    cg = core(l)
    return dict(family=l['family'], contig=l['contig'], start=l['start'], end=l['end'], idiomorph=l.get('idiomorph'),
                confidence=l.get('confidence'), locus_class=l.get('locus_class'), detection_pass=l.get('detection_pass'),
                genes=sorted(l.get('genes_found') or []), polished=l.get('polished_genes'),
                core={g: dict(ref=e.get('reference_record'), ident=e.get('identity'), cov=e.get('coverage'),
                              span=(e.get('start'), e.get('end')), nexon=len(e.get('exons') or []), status=e.get('status'),
                              bits=e.get('bitscore')) for g, e in cg.items()},
                partial_core=l.get('partial_core_note') or l.get('partial_core'))
def overlap(a, b): return a['contig'] == b['contig'] and a['start'] <= b['end'] and b['start'] <= a['end'] and a['family'] == b['family']
rows = []; stats = collections.Counter(); refdiff = collections.Counter()
for asm, m in meta.items():
    b, p = load('base', asm), load('pilot', asm)
    if b is None or p is None:
        stats['missing_' + ('base' if b is None else 'pilot')] += 1; continue
    stats['compared'] += 1
    g = m['group']
    for k in ('routing_mode', 'not_searched_reason', 'genetic_code'):
        if b.get(k) != p.get(k): rows.append((g, asm, m['species'], 'RUN_FIELD', k, b.get(k), p.get(k)))
    fb = sorted(f if isinstance(f, str) else json.dumps(f, sort_keys=True) for f in (b.get('families_attempted') or []))
    fp = sorted(f if isinstance(f, str) else json.dumps(f, sort_keys=True) for f in (p.get('families_attempted') or []))
    if fb != fp: rows.append((g, asm, m['species'], 'FAMILIES', '', fb, fp))
    rb = {l.split('|')[0][1:] for l in (P/'runs/base'/asm/'_reference.faa').read_text().splitlines() if l.startswith('>')}
    rp = {l.split('|')[0][1:] for l in (P/'runs/pilot'/asm/'_reference.faa').read_text().splitlines() if l.startswith('>')}
    nb = sum(1 for l in (P/'runs/base'/asm/'_reference.faa').read_text().splitlines() if l.startswith('>'))
    npl = sum(1 for l in (P/'runs/pilot'/asm/'_reference.faa').read_text().splitlines() if l.startswith('>'))
    refdiff[(nb, npl, tuple(sorted(rp - rb)), tuple(sorted(rb - rp)))] += 1
    for sect in ('detected', 'suppressed_loci'):
        LB = [summ(l) for l in b.get(sect) or []]; LP = [summ(l) for l in p.get(sect) or []]
        used = set()
        for lb in LB:
            j = next((i for i, lp in enumerate(LP) if i not in used and overlap(lb, lp)), None)
            if j is None:
                rows.append((g, asm, m['species'], sect.upper() + '_LOST', lb['family'], json.dumps(lb, default=str), '')); continue
            used.add(j); lp = LP[j]
            diffs = {f: (lb[f], lp[f]) for f in lb if lb[f] != lp[f]}
            if diffs:
                rows.append((g, asm, m['species'], sect.upper() + '_CHANGED', lb['family'] + ' ' + lb['contig'],
                             json.dumps({f: v[0] for f, v in diffs.items()}, default=str),
                             json.dumps({f: v[1] for f, v in diffs.items()}, default=str)))
        for i, lp in enumerate(LP):
            if i not in used:
                rows.append((g, asm, m['species'], sect.upper() + '_NEW', lp['family'], '', json.dumps(lp, default=str)))
    stats['det_base'] += len(b.get('detected') or []); stats['det_pilot'] += len(p.get('detected') or [])
    stats['sup_base'] += len(b.get('suppressed_loci') or []); stats['sup_pilot'] += len(p.get('suppressed_loci') or [])
with open(out, 'w') as fh:
    w = csv.writer(fh, delimiter='\t'); w.writerow(['group', 'asm', 'species', 'change', 'what', 'base', 'pilot'])
    w.writerows(rows)
print(dict(stats))
print('changes by type:', collections.Counter(r[3] for r in rows))
print('changes by group:', collections.Counter((r[0], r[3]) for r in rows))
print('reference fasta (n_base, n_pilot, added, removed) -> genomes:')
for k, v in refdiff.items(): print('  ', k, v)
