"""Build candidate metadata.yaml records for the Sporidiobolales redPR and
Pucciniales rustHD families from GenBank records. Coordinates, exons,
codon_start, transl_table and protein accessions are read from the deposit's
own CDS features, never typed by hand. Usage: build_records.py specs.yaml"""
import re, sys, time, urllib.request
import yaml
from Bio import SeqIO

D = '/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_puccinio_curation'
EMAIL = 'jason.stajich@ucr.edu'


def lineage(taxid):
    x = urllib.request.urlopen(f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxid}&retmode=xml&email={EMAIL}').read().decode()
    time.sleep(0.5)
    ranks = dict(re.findall(r'<ScientificName>([^<]+)</ScientificName>\s*<Rank>([^<]+)</Rank>', x))
    inv = {v: k for k, v in ranks.items()}
    name = re.search(r'<ScientificName>([^<]+)</ScientificName>', x).group(1)
    parts = [f'k__{inv.get("kingdom", "Fungi")}', f'p__{inv["phylum"]}']
    if 'subphylum' in inv:
        parts.append(f'sc__{inv["subphylum"]}')
    for r, p in (('class', 'c'), ('order', 'o'), ('family', 'f'), ('genus', 'g')):
        if r in inv:
            parts.append(f'{p}__{inv[r]}')
    parts.append('s__' + name.replace(' ', '_'))
    return ';'.join(parts), name


def gene_entry(i, f, name, role, offset):
    loc = f.location
    strand = '+' if loc.strand == 1 else '-'
    exons = [{'start': int(p.start) + 1, 'end': int(p.end)} for p in loc.parts]
    exons.sort(key=lambda e: -e['start'] if strand == '-' else e['start'])
    q = f.qualifiers
    return {'gene_index': i, 'name': name, 'protein_accession': 'ncbi_protein:' + q['protein_id'][0],
            'role': role, 'present': True, 'locus_tag': q.get('locus_tag', [None])[0], 'segment_index': 0,
            'start': int(loc.start) + 1, 'end': int(loc.end), 'strand': strand, 'order_in_locus': i,
            'exons': exons, 'codon_start': int(q.get('codon_start', ['1'])[0]),
            'transl_table': int(q.get('transl_table', ['1'])[0])}


def build(spec):
    gb = SeqIO.read(f'{D}/ncbi/{spec["accession"]}.gbk', 'genbank')
    src = [x for x in gb.features if x.type == 'source'][0].qualifiers
    taxid = int([x for x in src['db_xref'] if x.startswith('taxon:')][0].split(':')[1])
    cds = {f.qualifiers['protein_id'][0]: f for f in gb.features if f.type == 'CDS' and 'protein_id' in f.qualifiers}
    genes = [gene_entry(i, cds[pid], name, role, 0) for i, (pid, name, role) in enumerate(spec['genes'])]
    genes.sort(key=lambda g: g['start'])
    for i, g in enumerate(genes):
        g['order_in_locus'] = i
    lin, sci = lineage(taxid)
    start = min(g['start'] for g in genes)
    end = max(g['end'] for g in genes)
    tier = spec['tier']
    ev = lambda m: {'tier': tier, 'experimental_method': m, 'citations': spec['citations']}
    strain = spec['strain']
    rec = {
        'record_id': spec['record_id'], 'record_version': 1,
        'taxonomy': {'taxid': taxid, 'lineage': lin, 'lineage_resolved_date': '2026-09-26'},
        'organism': {'species': sci, 'strain': {'name': strain, 'known': True,
                     'culture_collection_ids': [strain], 'differs_from_sequenced': False}},
        'mating_type': {'locus_name': spec['locus'], 'idiomorphs': [spec['idiomorph']],
                        'system': 'heterothallic'},
        'locus': {'coordinate_provenance': spec.get('provenance', 'published_explicit'),
                  'excluded_from_coordinate_benchmark': True,
                  'core': {'completeness': spec['completeness'], 'reference_orientation': spec['orientation'],
                           'definition_note': spec['note'],
                           'segments': [{'segment_index': 0, 'sequence_source': {'type': 'insdc_nucleotide',
                                         'accession': gb.id, 'seq_region': gb.id}, 'start': start, 'end': end,
                                         'contig_edge_distance': None, 'sequence_checksum': None}]},
                  'extended_flank': []},
        'genes': genes,
        'evidence': {'locus_existence': ev(spec['methods'][0]), 'boundaries': ev(spec['methods'][1]),
                     'idiomorph_assignment': ev(spec['methods'][2])},
        'validation': {'status': 'needs_review', 'rejection_reason': None, 'accession_resolved': None,
                       'accession_resolved_date': None, 'accession_resolved_version': None,
                       'sequence_match': None, 'taxonomy_current': None},
        'curation': {'proposed_by': 'claude-literature-search', 'proposal_dedupe_key': spec['dedupe'],
                     'reviewed_by': None, 'reviewed_date': None,
                     'notes': 'Proposed 2026-09-26 on curator request (Sporidiobolales, Pucciniales). '
                              'PENDING CURATOR SIGN-OFF: accepted on the unmerged branch curation-puccinio '
                              'only so detect can be measured with it.'},
        'model_provenance': None,
    }
    out = f'{D}/records/{spec["record_id"]}.yaml'
    open(out, 'w').write(yaml.safe_dump(rec, sort_keys=False, width=100))
    print('wrote', spec['record_id'], taxid, sci, len(genes), 'genes', start, end)


if __name__ == '__main__':
    for s in yaml.safe_load(open(sys.argv[1])):
        build(s)
