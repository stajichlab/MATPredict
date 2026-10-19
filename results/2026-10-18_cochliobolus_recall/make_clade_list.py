import csv
rows = list(csv.DictReader(open('results/2026-10-18_dothideomycetes_audit/genome_audit_final.tsv'), delimiter='\t'))
lst = {r['ASM_ACCESSION']: r for r in csv.DictReader(open('results/2026-10-18_dothideomycetes_audit/genome_list.tsv'), delimiter='\t')}
GEN = ('Bipolaris', 'Curvularia', 'Exserohilum', 'Ophiosphaerella', 'Cochliobolus', 'Setosphaeria')
# the whole Pleosporaceae/Phaeosphaeriaceae that are no-call, as comparison
with open('results/2026-10-18_cochliobolus_recall/clade_genomes.tsv', 'w') as f:
    f.write('genome\ttaxid\tfamily\tgenus\tspecies\tstrain\tstatus\tclade\n')
    for r in rows:
        clade = r['genus'] in GEN
        other_nc = r['family'] in ('Pleosporaceae', 'Phaeosphaeriaceae') and r['status'] != 'called' and not clade
        if clade or other_nc:
            l = lst.get(r['genome'], {})
            f.write('\t'.join([r['genome'], r['taxid'], r['family'], r['genus'], r['species'], l.get('STRAIN', ''), r['status'], 'clade' if clade else 'other_nocall']) + '\n')
