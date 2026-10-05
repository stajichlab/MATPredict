"""Build candidate metadata.yaml files for the round-2 Dothideomycetes records.
Coordinates and protein ids are copied from the GenBank deposits (fetched
2026-09-26). Evidence text is limited to what the cited abstract states."""
import sys, yaml
from pathlib import Path

OUT = Path(sys.argv[1])
TODAY = '2026-09-26'
LIN = {
 5022: 'k__Fungi;p__Ascomycota;subphylum__Pezizomycotina;c__Dothideomycetes;o__Pleosporales;f__Leptosphaeriaceae;g__Plenodomus;s__Plenodomus_lingam',
 13684: 'k__Fungi;p__Ascomycota;subphylum__Pezizomycotina;c__Dothideomycetes;o__Pleosporales;f__Phaeosphaeriaceae;g__Parastagonospora;s__Parastagonospora_nodorum',
 5499: 'k__Fungi;p__Ascomycota;subphylum__Pezizomycotina;c__Dothideomycetes;o__Mycosphaerellales;f__Mycosphaerellaceae;g__Fulvia;s__Fulvia_fulva',
 1873960: 'k__Fungi;p__Ascomycota;subphylum__Pezizomycotina;c__Dothideomycetes;o__Mycosphaerellales;f__Mycosphaerellaceae;g__Pseudocercospora;s__Pseudocercospora_fijiensis',
}
SPECIES = {5022: 'Plenodomus lingam', 13684: 'Parastagonospora nodorum', 5499: 'Fulvia fulva', 1873960: 'Pseudocercospora fijiensis'}
CIT = {
 'lm': {'pmid': '12679880', 'doi': '10.1007/s00294-003-0391-6'},
 'pn': {'pmid': '12948511', 'doi': '10.1016/s1087-1845(03)00062-8'},
 'ff': {'pmid': '17178244', 'doi': '10.1016/j.fgb.2006.11.004'},
 'pf': {'pmid': '20507483', 'doi': '10.1111/j.1364-3703.2006.00376.x'},
}
EV = {
 'lm': ('Targeted sequencing of the MAT1-1 and MAT1-2 regions from isolates of each mating '
        'type, with transcript analysis of the MAT and flanking genes and a three-primer '
        'mating-type PCR assay (Cozijnsen & Howlett 2003, abstract). The abstract does not '
        'describe crosses; the full text is closed access and was not read.'),
 'pn': ('Idiomorphs identified "using DNA from a pair of isolates from Poland and Georgia, '
        'USA that are known to mate"; MAT-specific primers then typed field isolates of both '
        'mating types (Bennett et al. 2003, abstract).'),
 'ff': ('Cloning and sequencing of MAT1-1-1 and MAT1-2-1 from the "presumed asexual" '
        'Cladosporium fulvum; mating types typed in 86 strains (Stergiopoulos et al. 2007, '
        'abstract). No sexual state is known, so there is no cross; the evidence is targeted '
        'locus sequencing, like the accepted C. albicans MTL records.'),
 'pf': ('Targeted isolation of both idiomorphs: the MAT1-2 HMG box by degenerate PCR and DNA '
        'walking from the DNA lyase gene, MAT1-1 by long-range PCR from the flanks '
        '(Conde-Ferraez et al. 2007, abstract). The abstract does not describe crosses; the '
        'full text could not be retrieved.'),
}
R = [
 # key, taxid, strain, known, idiomorph, accession, genes[(name, prot, role, exons)], orientation, note
 ('lm', 5022, 'unknown-1', False, 'MAT1-1', 'AY174048.1',
  [('MAT1-1-1', 'AAO37757.1', 'core_MAT', [(4362, 4614), (4660, 5729)])],
  'GAP1->ORF1->MAT1-1-1->APN2',
  'Leptosphaeria maculans (NCBI Plenodomus lingam, taxid 5022) MAT1-1 region, deposit AY174048.1 '
  '(9,823 bp). The deposit names no strain. The paper calls the gene MAT1-1 (alpha box, 441 aa); '
  'it is MAT1-1-1. Curated: MAT1-1-1 (AAO37757.1). Not curated: GAP1 (AAO37755.1) and an unknown '
  'ORF (AAO37756.1), not roster genes; the DNA lyase (APN2, AAO37758.1) is partial at the deposit '
  'end with codon_start=2 and is left out.'),
 ('lm', 5022, 'unknown-2', False, 'MAT1-2', 'AY174049.1',
  [('MAT1-2-1', 'AAO37761.1', 'core_MAT', [(4780, 5245), (5301, 6034)])],
  'GAP1->ORF1->MAT1-2-1->APN2',
  'Leptosphaeria maculans (NCBI Plenodomus lingam, taxid 5022) MAT1-2 region, deposit AY174049.1 '
  '(10,124 bp). The deposit names no strain. MAT1-2 (HMG box, 397 aa) is MAT1-2-1. Curated: '
  'MAT1-2-1 (AAO37761.1). Not curated: GAP1 (AAO37759.1), an unknown ORF (AAO37760.1), and the '
  'partial DNA lyase (AAO37762.1).'),
 ('pn', 13684, 'sn435pl98', True, 'MAT1-1', 'AY212018.1',
  [('MAT1-1-1', 'AAO31740.1', 'core_MAT', [(1941, 2160), (2216, 3030)])],
  'ORF1->MAT1-1-1',
  'Parastagonospora nodorum (as Phaeosphaeria nodorum) isolate SN435PL98 (Poland) MAT1-1, deposit '
  'AY212018.1 (6,193 bp). Curated: MAT1-1-1 (AAO31740.1). Not curated: ORF1 (AAO31739.1), not a '
  'roster gene; per the paper the idiomorph begins inside ORF1.'),
 ('pn', 13684, 'sn436ga98', True, 'MAT1-2', 'AY212019.1',
  [('MAT1-2-1', 'AAO31742.1', 'core_MAT', [(2272, 2761), (2815, 3362)])],
  'ORF1->MAT1-2-1',
  'Parastagonospora nodorum isolate SN436GA98 (Georgia, USA) MAT1-2, deposit AY212019.1 (6,401 bp). '
  'Curated: MAT1-2-1 (AAO31742.1). Not curated: ORF1 (AAO31741.1).'),
 ('ff', 5499, 'alenya-b', True, 'MAT1-1', 'DQ659350.2',
  [('MAT1-1-1', 'ABG45907.1', 'core_MAT', [(2156, 2404), (2453, 2544), (2593, 2761), (2815, 3474)])],
  'ORF1-1-2->MAT1-1-1',
  'Fulvia fulva (syn. Cladosporium fulvum, Passalora fulva) strain Alenya B MAT1-1 idiomorph, deposit '
  'DQ659350.2 (5,433 bp; idiomorph 662..4823). Curated: MAT1-1-1 (ABG45907.1). Not curated: '
  'ORF1-1-2 (ABK41478.1), not a roster gene. NCBI places Fulvia in Mycosphaerellaceae.'),
 ('ff', 5499, 'imi-day9-054980', True, 'MAT1-2', 'DQ659351.2',
  [('MAT1-2-1', 'ABG49507.1', 'core_MAT', [(2140, 2312), (2364, 2605), (2658, 2728), (2781, 3449)])],
  'ORF1-2-2->MAT1-2-1',
  'Fulvia fulva strain IMI Day9 054980 MAT1-2 idiomorph, deposit DQ659351.2 (6,344 bp; idiomorph '
  '881..4434). Curated: MAT1-2-1 (ABG49507.1). Not curated: ORF1-2-2 (ABK41479.1).'),
 ('pf', 1873960, 'unknown-1', False, 'MAT1-1', 'DQ787015.1',
  [('MAT1-1-1', 'ABH04239.1', 'core_MAT', [(2524, 2773), (2825, 2918), (2974, 3142), (3192, 3845)])],
  'MAT1-1-1',
  'Pseudocercospora fijiensis (as Mycosphaerella fijiensis) MAT1-1, deposit DQ787015.1 (5,243 bp; '
  'idiomorph 1004..4877). The deposit names no strain. Curated: MAT1-1-1 (ABH04239.1).'),
 ('pf', 1873960, 'unknown-2', False, 'MAT1-2', 'DQ787016.1',
  [('APN2', 'ABH04240.1', 'flanking_conserved', [(895, 2763)]),
   ('MAT1-2-1', 'ABH04241.1', 'core_MAT', [(6550, 6734), (6768, 7129), (7180, 7955)])],
  'APC5->APN2->MAT1-2-1',
  'Pseudocercospora fijiensis MAT1-2 with flanks, deposit DQ787016.1 (9,799 bp; idiomorph 4848..9254). '
  'The deposit names no strain. Curated: APN2 (annotated "DNA lyase-like protein", ABH04240.1) and '
  'MAT1-2-1 (ABH04241.1). Not curated: the partial anaphase-promoting complex protein (ABH04242.1), '
  'not a roster gene.'),
]
for key, taxid, strain, known, idio, acc, genes, orient, note in R:
    rid = f'{taxid}_{strain}_MAT_{idio}'
    gl = []
    for i, (name, prot, role, exons) in enumerate(genes):
        gl.append({'gene_index': i, 'name': name, 'protein_accession': f'ncbi_protein:{prot}',
                   'role': role, 'present': True, 'locus_tag': None, 'segment_index': 0,
                   'start': exons[0][0], 'end': exons[-1][1], 'strand': '+', 'order_in_locus': i,
                   'exons': [{'start': s, 'end': e} for s, e in exons]})
    cit = CIT[key]
    rec = {
     'record_id': rid, 'record_version': 1,
     'taxonomy': {'taxid': taxid, 'lineage': LIN[taxid], 'lineage_resolved_date': TODAY},
     'organism': {'species': SPECIES[taxid], 'strain': {'name': strain if known else 'unknown', 'known': known,
                  'culture_collection_ids': [], 'differs_from_sequenced': False}},
     'mating_type': {'locus_name': 'MAT', 'idiomorphs': [idio], 'system': 'heterothallic'},
     'locus': {'coordinate_provenance': 'published_explicit', 'excluded_from_coordinate_benchmark': False,
               'core': {'completeness': 'complete', 'reference_orientation': orient, 'definition_note': note,
                        'segments': [{'segment_index': 0,
                                      'sequence_source': {'type': 'insdc_nucleotide', 'accession': acc, 'seq_region': acc},
                                      'start': gl[0]['start'], 'end': gl[-1]['end'], 'contig_edge_distance': None}]},
               'extended_flank': []},
     'genes': gl,
     'evidence': {
       'locus_existence': {'tier': 1, 'experimental_method': EV[key], 'citations': [cit]},
       'boundaries': {'tier': 1, 'experimental_method': "Boundaries are the GenBank deposit's own; the record segment spans the curated genes on it.", 'citations': [cit]},
       'idiomorph_assignment': {'tier': 1, 'experimental_method': 'Idiomorph from MAT gene content (alpha-box MAT1-1-1 vs HMG-box MAT1-2-1) of the targeted locus sequence.', 'citations': [cit]}},
     'validation': {'status': 'needs_review', 'rejection_reason': None, 'accession_resolved': False,
                    'accession_resolved_date': None, 'accession_resolved_version': None,
                    'sequence_match': None, 'taxonomy_current': True},
     'curation': {'proposed_by': 'claude-literature-search',
                  'proposal_dedupe_key': f"{cit['doi']}|{taxid}|MAT|{idio}",
                  'reviewed_by': None, 'reviewed_date': None,
                  'notes': ('Proposed 2026-09-26 on curator request (round 2: Leptosphaeria, Parastagonospora, '
                            'Pseudocercospora, Fulvia) from docs/notes/2026-09-24_mat-reference-gap-literature.md. '
                            'PENDING CURATOR SIGN-OFF: accepted on the unmerged branch curation-mucor-dothideo only so '
                            'detect can be measured with it.')},
     'model_provenance': None,
    }
    p = OUT / f'{rid}.yaml'
    p.write_text(yaml.safe_dump(rec, sort_keys=False))
    print(p)
