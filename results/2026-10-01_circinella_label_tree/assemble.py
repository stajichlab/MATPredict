"""Assemble the tree protein set and tip metadata (tree_set.faa, tips.tsv)."""
import csv, re
R = '/bigdata/stajichlab/jstajich/projects/MATPredict/results'
P1 = '/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db/Mucoromycota/classifiers/MAT/paralogs/P1.faa'
def rf(p):
    s, n = {}, None
    for l in open(p):
        if l.startswith('>'): n = l[1:].strip(); s[n] = ''
        else: s[n] += l.strip()
    return s
G = {r['gid']: r for r in csv.DictReader(open('genomes.tsv'), delimiter='\t')}
cur = {r['genome_id'] if 'genome_id' in r else r.get('gid'): r for r in csv.DictReader(open(f'{R}/2026-09-28_lcg_holdout/curator_table.tsv'), delimiter='\t')} if False else {}
seqs, tips = {}, []
def add(tip, seq, **kw):
    seqs[tip] = seq.replace('*', '')
    tips.append(dict(tip=tip, len=len(seqs[tip]), **kw))
# Circinella group: rnhA-adjacent model (miniprot, full length)
rm = rf('rnhA_models.faa')
rmt = {r['gid']: r for r in csv.DictReader(open('rnhA_models.tsv'), delimiter='\t')}
loci = rf('loci.faa'); mpq = rf('mpq_best.faa'); extra = rf('extra.faa')
fallback = {  # genomes with no gene modelled at rnhA: source, note
    'Circinella_muscae_NRRL_1360': (loci, 'Circinella_muscae_NRRL_1360', 'standard sexM 0.9 kb from rnhA (annot)'),
    'Circinella_muscae_NRRL_1362': (mpq, 'Circinella_muscae_NRRL_1362|sexPtype', 'sexP-type on a scaffold without rnhA'),
    'Circinella_naumovii_NRRL_5846': (loci, 'Circinella_naumovii_NRRL_5846', 'standard sexM 0.1 kb from rnhA'),
    'Circinella_simplex_van_Tieghem_NRRL_2407': (loci, 'Circinella_simplex_van_Tieghem_NRRL_2407', 'standard sexP 0.3 kb from rnhA'),
    'Circinella_simplex_CBS_142.35': (loci, 'Circinella_simplex_CBS_142.35', 'partial 103 aa, not at rnhA'),
    'Circinella_tenella_NRRL_A-23557': (extra, 'Circinella_tenella_NRRL_A-23557|sexPlike', 'exonerate model, not at rnhA'),
    'Thamnostylum_piriforme_NRRL_A-21589': (mpq, 'Thamnostylum_piriforme_NRRL_A-21589|sexMtype', 'sexM-type on a scaffold without rnhA'),
    'Thamnostylum_repens_Tieghem_Upadhyay_NRRL_6240': (loci, 'Thamnostylum_repens_Tieghem_Upadhyay_NRRL_6240', 'misidentified: Circinomucor circinelloides'),
    'GCA_982397305.1_T17-F': (mpq, 'GCA_982397305.1_T17-F|sexMtype', 'partial 143 aa'),
}
for gid, g in G.items():
    if g['group'] != 'circinella_group':
        continue
    if gid in rm:
        add('CG|' + gid, rm[gid], group='circinella_group', gid=gid, species=g['species'], label=g['label'],
            source='rnhA_model:' + rmt[gid]['query'], note=f"{rmt[gid]['dist_rnhA_kb']} kb from rnhA")
    elif gid in fallback:
        src, k, note = fallback[gid]
        add('CG|' + gid, src[k], group='circinella_group', gid=gid, species=g['species'], label=g['label'], source='fallback', note=note)
    else:
        print('MISSING', gid)
# labelled / unlabelled reference strains (skip those already in Zygo)
zy = rf('zygo_best.faa')
zorgs = {k.split('|')[1] for k in zy}
for gid, g in G.items():
    if g['group'] == 'circinella_group' or gid in zorgs:
        continue
    k = gid if gid in loci else None
    if k is None and f'{gid}|sexMlike' in extra:
        k = None
    if gid in loci:
        add('LAB|' + gid, loci[gid], group=g['group'], gid=gid, species=g['species'], label=g['label'], source='loci.faa', note='')
    else:
        print('no seq', gid)
for k, s in zy.items():
    _, org, lab, pid = k.split('|')
    add(f'ZYGO|{org}|{pid}', s, group='zygo_truth', gid=org, species=org, label=lab, source='zygo_best', note='')
for k, s in rf(f'{R}/2026-10-01_lichtheimiaceae/sexMP_refs.faa').items():
    add('REF|' + k.split()[0], s, group='curated_ref', gid=k.split('|')[0], species='', label=k.split('|')[2], source='sexMP_refs', note='')
for k, s in rf(P1).items():
    add('OUT|P1|' + k.split()[0], s, group='outgroup', gid=k, species='Mucor indicus P1 paralog', label='P1', source='P1.faa', note='HMG box only')
for k, s in rf(f'{R}/2026-09-26_sexMP_phylogeny/refs_MAT121.faa').items():
    t = '|'.join(k.split('|')[:3])
    add(t, s, group='outgroup', gid=k.split('|')[2], species='Fusarium MAT1-2-1' if '5518' in k else 'Tuber MAT1-2-1',
        label='MAT1-2-1', source='refs_MAT121', note='root' if '5518' in k else '')
with open('tree_set.faa', 'w') as f:
    for k, s in seqs.items():
        f.write(f'>{k}\n{s}\n')
with open('tips.tsv', 'w') as f:
    cols = ['tip', 'group', 'gid', 'species', 'label', 'source', 'note', 'len']
    f.write('\t'.join(cols) + '\n')
    for t in tips:
        f.write('\t'.join(str(t.get(c, '')) for c in cols) + '\n')
print(len(seqs))
