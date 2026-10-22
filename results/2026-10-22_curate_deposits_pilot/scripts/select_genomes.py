"""Fixed-seed, stratified genome selection for the before/after replay."""
import csv, random, collections
avail = {l.strip()[:-6] for l in open('bfd_genomes.txt') if l.strip().endswith('.fa.gz') and not l.strip().endswith('masked.fasta.gz')}
rows = [r for r in csv.DictReader(open('ncbi_accessions_taxonomy.csv')) if r['ASM_ACCESSION'] in avail and r['NCBI_TAXID'].isdigit()]
rng = random.Random(20261022)
chosen = {}
def take(group, pool, n, by='GENUS'):
    pool = [r for r in pool if r['ASM_ACCESSION'] not in chosen]
    strata = collections.defaultdict(list)
    for r in sorted(pool, key=lambda r: r['ASM_ACCESSION']):
        strata[r[by] or r['SPECIES']].append(r)
    for v in strata.values(): rng.shuffle(v)
    keys = sorted(strata); rng.shuffle(keys)
    got = 0
    while got < n and any(strata[k] for k in keys):
        for k in keys:
            if strata[k] and got < n:
                r = strata[k].pop(); chosen[r['ASM_ACCESSION']] = (group, r); got += 1
    return got
pilot_species = ["Beauveria bassiana","Ophiocordyceps sinensis","Ustilaginoidea virens","Calonectria naviculata",
  "Cordyceps militaris","Diplocarpon coronariae","Monilinia laxa","Hymenoscyphus fraxineus","Ascochyta rabiei",
  "Diplodia sapinea","Lasiodiplodia theobromae","Berkeleyomyces basicola","Cryphonectria parasitica","Didymella pinodes"]
for asm in ['GCF_000280675.1_ASM28067v1','GCA_008080495.1_ASM808049v1','GCA_009829825.1_ASM982982v1','GCA_003031705.1_ASM303170v1']:
    r = next((r for r in rows if r['ASM_ACCESSION']==asm), None)
    if r: chosen[asm] = ('pilot_species', r)
    else: print('deposit-strain genome not available', asm)
for sp in pilot_species:
    pool = [r for r in rows if r['SPECIES'] == sp]
    # prefer the deposit's own strain if present
    take('pilot_species', pool, 3, by='STRAIN')
take('Hypocreales', [r for r in rows if r['ORDER']=='Hypocreales'], 60)
take('Helotiales', [r for r in rows if r['ORDER']=='Helotiales'], 60)
take('other_pilot_orders', [r for r in rows if r['ORDER'] in ('Pleosporales','Botryosphaeriales','Diaporthales','Microascales')], 16, by='ORDER')
pilot_orders = {'Hypocreales','Helotiales','Pleosporales','Botryosphaeriales','Diaporthales','Microascales'}
asc = [r for r in rows if r['PHYLUM']=='Ascomycota' and r['ORDER'] not in pilot_orders and r['ORDER']]
yeast_classes = {'Saccharomycetes','Pichiomycetes','Dipodascomycetes','Alloascoideomycetes','Lipomycetes','Trigonopsidomycetes','Sporopachydermiomycetes','Schizosaccharomycetes','Taphrinomycetes','Pneumocystomycetes','Neolectomycetes','Archaeorhizomycetes'}
pez = [r for r in asc if r['CLASS'] not in yeast_classes]
big = [o for o,_ in collections.Counter(r['ORDER'] for r in pez).most_common(8)]
for o in ['Saccharomycetales','Serinales']:
    take('control_major_orders', [r for r in asc if r['ORDER']==o], 3)
print('major orders', big)
for o in big:
    take('control_major_orders', [r for r in asc if r['ORDER']==o], 5)
take('control_ascomycota_other', asc, 40, by='ORDER')
take('control_basidiomycota', [r for r in rows if r['PHYLUM']=='Basidiomycota'], 12, by='ORDER')
take('control_mucoromycota', [r for r in rows if r['PHYLUM']=='Mucoromycota'], 8, by='ORDER')
with open('replay_genomes.tsv','w') as fh:
    fh.write('asm\ttaxid\tgroup\tspecies\torder\n')
    for a,(g,r) in chosen.items():
        fh.write(f"{a}\t{r['NCBI_TAXID']}\t{g}\t{r['SPECIES_IN']}\t{r['ORDER']}\n")
c = collections.Counter(g for g,_ in chosen.values()); print(c, sum(c.values()))
print(collections.Counter(r['ORDER'] for g,r in chosen.values() if g.startswith('control_')))
print([ (r['SPECIES_IN']) for g,r in chosen.values() if g=='pilot_species'])
