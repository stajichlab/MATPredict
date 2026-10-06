"""Per-strain table: label, detection call, locus gene, classifier scores of the tree protein,
clade placement per tree, verdict.

placement = reference types in the smallest clade holding the tip and >=1 curated/Zygo
reference; support = that clade's support. half = support of the MRCA of all sexP
(or sexM) references when the tip sits inside it ('no' otherwise).
Verdict (labelled strains): matches / inverted when every tree agrees; else unresolved.
Usage: strain_table.py name=clades.tsv ...  > strain_table.tsv
"""
import csv, sys
trees = [a.split('=') for a in sys.argv[1:]]
sys.argv = [sys.argv[0]]
exec(open('extract.py').read().split('rows = list(')[0])
seqs = read_fasta('tree_set.faa')
G = {r['gid']: r for r in csv.DictReader(open('genomes.tsv'), delimiter='\t')}
calls = {r['genome_id']: r for r in csv.DictReader(open('../2026-10-01_lichtheimiaceae/all.tsv'), delimiter='\t')}
rm = {r['gid']: r for r in csv.DictReader(open('rnhA_models.tsv'), delimiter='\t')}
tips = {r['tip']: r for r in csv.DictReader(open('tips.tsv'), delimiter='\t')}
loci = {r['gid']: r for r in csv.DictReader(open('loci.tsv'), delimiter='\t')}
pl = {}
for name, f in trees:
    for r in csv.DictReader((l for l in open(f) if not l.startswith('#')), delimiter='\t'):
        half = r['in_sexP_ref_mrca'] if r['in_sexP_ref_mrca'] != 'no' else r['in_sexM_ref_mrca']
        pl.setdefault(r['tip'], {})[name] = (r['placement'], r['support'] or 'NA', half or 'NA')
cols = ['strain', 'group', 'label', 'label_disputed', 'detection_call', 'rnhA_gene', 'dist_rnhA_kb', 'len', 'sexM', 'sexP', 'P1', 'sexP_minus_sexM']
for n, _ in trees:
    cols += [f'{n}_placement', f'{n}_support', f'{n}_refhalf_support']
cols += ['verdict', 'note']
print('\t'.join(cols))
for tip, r in tips.items():
    if r['group'] not in ('circinella_group', 'reference_labelled', 'reference_mucorales'):
        continue
    gid = r['gid']; g = G[gid]
    m = rm.get(gid, {})
    if r['source'].startswith('rnhA_model'):
        gene = 'sexP-type' if m['query'].startswith('sexP') else 'sexM-type'; d = m['dist_rnhA_kb']
    elif r['source'] == 'fallback':
        gene, d = 'other', ''
    else:
        l = loci.get(gid, {})
        gene = 'standard ' + ('sexP' if l.get('mp_query', '').endswith('sexP') else 'sexM'); d = l.get('dist_to_rnhA_kb', '')
    sc = score(seqs[tip])
    pt = pl.get(tip, {})
    vals = []
    for n, _ in trees:
        vals += list(pt.get(n, ('NA', 'NA', 'NA')))
    places = {pt[n][0] for n, _ in trees if n in pt}
    lab = g['label']; v = ''
    if lab in ('Plus', 'Minus'):
        want = 'sexP' if lab == 'Plus' else 'sexM'
        v = ('label matches phylogeny' if places == {want} else 'label inverted') if len(places) == 1 and places <= {'sexP', 'sexM'} else 'unresolved'
    note = '; '.join(x for x in (r['note'] if r['source'] == 'fallback' else '', g.get('notes', '')) if x)
    c = calls.get(gid, {})
    print('\t'.join(map(str, [gid, g['group'], lab, c.get('label_disputed', ''), c.get('new_call', ''), gene, d, len(seqs[tip]),
                              sc['sexM'], sc['sexP'], sc['P1'], f"{sc['sexP'] - sc['sexM']:.1f}"] + vals + [v, note])))
