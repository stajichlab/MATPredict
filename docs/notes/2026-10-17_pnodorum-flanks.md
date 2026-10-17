# P. nodorum flank records (2026-10-17)

Branch `pnodorum-flanks` (from origin/main 5b9a9e0). One new record, no code change, no change to the split-locus limits.
Runs: HPCC, pixi env mat-prof, worktree synced to /bigdata/stajichlab/jstajich/projects/pnod, jobs `pnod-*`, at most 9-10 tasks at once.
Raw outputs: /bigdata/stajichlab/jstajich/projects/pnod-runs/{base,new,sla2,minhits1_*,flank}/ (reports, evidence_diagnostics.jsonl).

## 1. Short answer

- Added `13684_sn15_MAT_MAT1-1` (P. nodorum SN15, RefSeq GCF_000146915.1, scaffold NW_001884549.1:865888-873147): MAT1-1-1, APN2, COX13. validate: accession resolved, 3 of 3 genes 100% identity / 100% coverage, taxonomy_current true. The two AY212018/AY212019 records are untouched (no `record_version` bump).
- **SLA2 is not in the record.** In SN15 SLA2 (SNOG_00689, XP_001791366.1) is 690 kb from MAT1-1-1 on the same scaffold, and also in the chromosome-level SN15 assembly GCA_016801405.1 (CP069023.1: MAT1-1-1 at 1.879 Mb, APN2 1.884 Mb, COX13 1.886 Mb, SLA2 2.571 Mb; my tblastn of the SN15 proteins, not an annotation read). So SLA2 is not "beside the MAT locus" in P. nodorum. The brief assumed it was.
- Effect: flank identity of the Parastagonospora flank clusters rises to 98.8-100% (from 86-89%), but **none of the four test genomes gets a call**, and no call anywhere in the 44-genome panel is gained, lost or changed. Flank records alone cannot fire the split-locus rule (section 5).

## 2. Assemblies (NCBI Assembly esummary, 2026-10-17)

| need | candidate | in the 186 / test set? | use |
|---|---|---|---|
| SN15 | GCF_000146915.1 (Broad ASM14691v2, RefSeq, scaffolds, annotated SNOG_*) = GCA_000146915.1 (status "replaced") | no | **used** |
| SN15 | GCA_000146915.2 | does not exist in the esummary (only .1) | n/a, flagged |
| SN15 | GCA_008452785.1 (ANU chromosome), GCA_051013865.1 (ANU complete) | **yes, among the 186** | not used |
| SN15 | GCA_016801405.1 (Curtin chromosome, annotated) | no | read only to check SLA2 position (tblastn) |
| Sn4 | GCA_002267005.1 (LDN03-Sn4) | **yes** | not used |
| Sn79 | GCA_002267025.1 (Sn79-1087) | **yes** | not used |
| Sn79 | GCA_002216185.1 (Curtin Sn79_1.0, N50 17 kb) | no | no gene annotation (0 genes). tblastn: MAT1-2-1 (98.9-99.4%) and APN2 (98.3%) on one contig, COX13 and SLA2 elsewhere. Not curated |

So no MAT1-2 same-strain record with flanks was possible without a test genome or hand annotation. APN2 and COX13 lie outside the idiomorph (and the proteins are 98-100% identical between SN15 and the test genomes), so the MAT1-1 record's flanks serve both idiomorphs. SN15 is not among the 186, so no `--exclude-records` was needed.

## 3. Decisions and conventions

- New same-strain record rather than editing the AY records: those are locus-specific deposits of other strains (SN435PL98, SN436GA98); SN15 flank coordinates cannot be placed on them. Pattern follows `27337_vdls-17_MAT_MAT1-2` (RefSeq-derived, `coordinate_provenance: curator_derived`, `assembly_accession_basis: sequence_link`). B10 (>=99% over >=90% of each core segment) is met trivially: the record is cut from that assembly, validate gives 100/100 per gene.
- Core gene MAT1-1-1 (XP_001791062.1, 344 aa) is 100% identical to AAO31740.1 (the SN435PL98 protein). SN15 = MAT1-1 per Syme et al. 2013 (cited from the brief; not re-read).
- APN2 and COX13 found by blastp against the Af293 records (50.6%, 48.2%) and by order (MAT1-1-1, two small genes, APN2, COX13). SNOG_00377 (APN2) is annotated 5'-partial.
- The record is `accepted` through the normal propose/validate/accept path; `reviewed_by` is null (no curator sign-off yet). validate first ran without taxonkit on PATH (taxonomy_current null, explicit error); re-run with the pixi env bin on PATH gave taxonomy_current true. I did not hand-set the flag. A stale `taxonomy_error` key from the first run was deleted from the metadata.
- Roster: `db/Ascomycota/order.yml`, locus MAT (scope 716546) declares SLA2, APN2 and COX13 as `flanking_conserved`. Covered; no change.

## 4. Tests

tests/db + tests/detect on the HPCC (job 30034518, test env): 1116 passed, 4 skipped.

## 5. Detect before/after

Base = origin/main db; new = main + the record. Same code. `--taxid` from the campaign list, MATPREDICT_HTML=0. 44 genomes x 2 = 88 runs, all rc 0: the 4 test genomes, 20 other Parastagonospora (seeded sample, seed 20261017, from the 186), 20 other Pleosporales (one per genus, seeded).

**Test genomes: calls and flank identity**

| genome | species | call base -> new | flank cluster best identity base -> new | suppressed loci |
|---|---|---|---|---|
| GCA_003497985.1 | P. nodorum | none -> none | 88.6 -> 100.0 (APN2/COX13 contig QVWR01000146.1) | 38 -> 38 |
| GCA_003503065.1 | P. nodorum | none -> none | 86.4 -> 98.8 (QVWJ01000028.1) | 13 -> 13 |
| GCA_022559985.1 | P. avenae | none -> none | 88.6 -> 91.5 (JAFJTR010000071.1) | 39 -> 40 |
| GCA_003503165.1 | P. avenae | none -> none | 88.6 -> 91.5 (QVWU01000052.1) | 7 -> 7 |

Per-gene best tblastn HSP identity (my own run, all flank proteins of the db, best bitscore): APN2 48.7% -> 98.8% (nodorum), 91.5% (avenae); COX13 35.3% -> 74.1%; SLA2 74.4% -> 74.4% (unchanged, no SLA2 record). COX13 stays at 74.1% because the best HSP spans introns (170 aligned columns for a 141 aa query; the SN15 protein against the SN15 chromosome also gives 74.1%). detect's own cluster/polished identity for COX13 is 100%, but the split-locus rule uses the HSP value.

**Why no call.** (a) The lone MAT1-2-1 cluster (99.4%) has one distinct non-superseded gene, so it fails the evidence floor (min_hits 2) and is never polished; the split-locus rule needs a modelled core gene. (b) Split-locus needs 2 flanks within 10 points of the core: with APN2 98.8%, COX13 74.1% and SLA2 74.4% only one qualifies. (c) GCA_003497985.1 is 2.2 kb from the contig end (limit 500 bp), unchanged by design.

**Side-effect panel.** Parastagonospora (20): every call, idiomorph, confidence, locus class and suppressed count identical; APN2/COX13 identity in the called locus 64.0/51.2% -> 98.7-100/100%. Other Pleosporales (20): 0 calls gained or lost; 8 genomes have extra withheld loci (+1 to +4; Lojkania +2, Setomelanomma +1, Beverwykella +2, Clohesyomyces +2, Melanomma +1, Foliophoma +1, Lophiotrema +2, Westerdykella +4; Pseudomonodictys -1); locus coordinates move by 3-6 bp in 7; flank identities rise in the called loci (COX13 to 62-92%, APN2 to 53-85%), idiomorph unchanged. Same pattern as the Verticillium note (new references add weak withheld clusters); not traced to specific hits.

## 6. SLA2 variant (not committed)

Variant db (HPCC /bigdata/stajichlab/jstajich/projects/pnod-sla2/db): the same record plus SLA2 as a second segment (NW_001884549.1:1557598-1563291, XP_001791366.1, 1609 aa, annotated partial at both ends), `order_in_locus` 3. It builds and passes the gene writer.
- Default settings, 44 genomes: no call changes; withheld loci +1 to +2 in most Parastagonospora genomes (counts not tabulated for the other Pleosporales); one called Clathrospora locus gains a 37% SLA2 gene.
- With `--min-hits 1` (a diagnostic, not proposed): GCA_003503065.1 gets MAT1-2 `partial_locus`, confidence low, `split_locus` (core 100%, 317 bp from the edge; APN2 98.8%, SLA2 90.4%, i.e. just inside the 90.0 limit). GCA_003497985.1 stays uncalled. Without the SLA2 record, `--min-hits 1` also gave no call.
So recovering even one genome needs both the SLA2 reference and a relaxed evidence floor, and the SLA2 reference rests on a gene 690 kb from MAT.

## 7. Not verified

- Syme et al. 2013 idiomorph assignments (SN15 MAT1-1, Sn4/Sn79 MAT1-2): taken from the brief.
- SLA2 distance: tblastn positions in two SN15 assemblies; no annotation of GCA_016801405.1 was read. Whether SLA2 is syntenic in other P. nodorum strains or P. avenae is not known.
- The 74.1% COX13 HSP artefact: inferred from alignment length, not traced inside detect.
- Withheld-locus increases in other Pleosporales were not traced to hits.
- P. avenae (2 of the 4 test genomes) has only the P. nodorum flanks at 91.5%; no P. avenae reference exists.
- GCA_003503165.1 (N50 5 kb) is an assembly gap, unchanged.

## 8. Questions for the curator

1. SLA2 is 690 kb from MAT in SN15 (two assemblies). Keep it out (as committed), or add it as a second segment of the record?
2. Accept that flank records alone give no new call, or authorise a separate change (evidence floor for a lone >=95% core with strong flanks; the earlier trace proposed `strong_single_identity`)?
3. COX13 never exceeds 74% at HSP level because of introns. Should the split-locus rule use the polished identity for flanks? (Code change, not made.)
4. Sign off `13684_sn15_MAT_MAT1-1` (reviewed_by is null)?
5. Is GCA_000146915.2 a different accession you meant (NCBI shows only .1 and RefSeq GCF_000146915.1)?
6. A MAT1-2 same-strain flank record needs a non-test, annotated Sn4 or Sn79 assembly; none exists outside the 186 (Curtin Sn79 is unannotated). Skip, or annotate GCA_002216185.1?
