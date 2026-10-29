# Botrytis genomes labelled MAT1-2 in BioSample but called MAT1-1 by v0.6.0

Date 2026-10-17, branch `botrytis-mislabels`. No change to `src/` or `db/`.

## Verdict

(b) The software call is wrong; the BioSample label is right. All 17 genomes carry the B. cinerea MAT1-2 idiomorph (MAT1-2-1 plus MAT1-2-4, now MAT1-2-10) at 98-100% identity to GenBank KF944386.1 (strain SAS405, `type: MAT1-2`). The step that fails is the reference set feeding the vote: the curated DB has a B. cinerea MAT1-1 record (KF944385.1) but no B. cinerea MAT1-2 record, so the real MAT1-2-1 matches a distant Sclerotinia protein at 77% (222 bits) while a 5'-truncated MAT1-1-1 remnant matches its close B. cinerea reference at 91-95% (~415 bits) and wins the bitscore vote. Confidence: 95% for "label right, software wrong"; 90% for the minimal fix below.

## Evidence

Method: miniprot (`-u -N 3`) of KF944385.1 (MAT1-1-5 AHX22633.1, MAT1-1-1 AHX22632.1) and KF944386.1 (MAT1-2-1 AHX22634.1, MAT1-2-4 AHX22635.1) proteins, plus the curated Sclerotinia APN2/SLA2, against all 32 BioProject PRJNA624742 assemblies, T4 (GCA_000292645.1) and B05.10 (GCF_000143535.2). Table: `results/2026-10-17_botrytis_mislabels/idiomorph_gene_coverage.tsv` (reference coverage / identity per gene).

1. Locus structure. In all 17 the locus is on one contig, in the order SLA2 - MAT1-1-1 remnant - MAT1-2-4 - MAT1-2-1 - APN2 (COX13 beyond APN2). This is the T4 layout (JH795966.1: APN2 691034, MAT1-2-1 693735-694982, MAT1-2-4 695728-696439, MAT1-1-1 remnant 696591-697284, SLA2 698934). T4 is therefore MAT1-2 and B05.10 (NC_037310.1: APN2 810328, MAT1-1-5 813478-814778, MAT1-1-1 815137-816297, SLA2 817947) is MAT1-1; the other 15 BioProject isolates have the B05.10 layout.
2. The MAT1-2-1 gene is real. 3 exons, 376 aa, 98.4-100% identical to AHX22634.1; the HMG-box profile scores it 220-223 bits and the alpha profile 0. It is the same gene the campaign reported as "MAT1-2-1 77%" (same coordinates, e.g. Rf1 JACVFJ010000135.1:30178-31416). The "MAT1-2-10" at 65% is MAT1-2-4/MAT1-2-10 (99% to AHX22635.1, 204 aa). Gene name: MAT1-2-4 (Doughan and Rollins 2016, Fungal Biol 120:1105, doi:10.1016/j.funbio.2016.06.007) was renamed MAT1-2-10 per Wilken et al. 2017 (as recorded in `db/Ascomycota/order.yml`; I did not re-verify Wilken in PubMed).
3. The MAT1-1-1 "hit" is a remnant. It aligns reference residues 129/140-353 only (61% of the 353 aa reference; 33% in Rf2, which is split across contigs), 215 aa, 91-96% identity. It lacks the N-terminus (residues 1-139), including the alpha-box core (`YILFR`). Alpha-box profile 119-127 bits versus 311 bits for the intact B05.10 protein; HMG profile 0. The same fragment lies unannotated inside KF944386.1 itself (about nt 3,270-3,970), so SAS405 has it too. The mirror image exists in MAT1-1: B05.10 carries a 42-aa MAT1-2-1 relic (813023-813148, 83%, internal stop), as Wilson et al. 2023 (doi:10.1371/journal.pone.0292619) describe. These are relics of the homothallic ancestral arrangement (cf. S. sclerotiorum, where both idiomorphs are fused).
4. Profile scores for the campaign's own proteins agree with the biology: HMG 221 bits (MAT1-2-1) beats alpha 125 bits (MAT1-1-1 remnant) by about 96 bits, so the core-profile classifier would call MAT1-2.
5. 32/32 concordance. Every assembly has exactly one idiomorph: 17 have full-length MAT1-2-1 and MAT1-2-4 and no MAT1-1-5; 15 have full-length MAT1-1-5 and MAT1-1-1 and no MAT1-2-4. This matches the BioSample labels for all 32. No genome has both.
6. Not a mixed or heterokaryotic sample. Assemblies are haploid-looking with one idiomorph each. Read check: k-mer depth of MAT1-1-5 and the MAT1-1-1 N-terminus in the first 4M read-1 reads of Rf1 (SRR11527080) is 0.09 and 0.06 versus 2.8-3.6 for MAT1-2-1, MAT1-2-4 and an SLA2 control (about 2-3%, the background of 25-mer matching), so there is no MAT1-1 nucleus population. The same check on Pv1 and Vv1 stalled on the ENA FTP stream and was cancelled; only Rf1 was checked at read level.

## Provenance

BioProject PRJNA624742, INRAE BIOGER (A. Simon), "Genomic determinants of host specialization in Botrytis cinerea"; 32 isolates (Vv grapevine, Sl tomato, Rf bramble, Hm, Pv1 = B. fabae from bean, Tunisia). Publication: Mercier et al. 2021, Phytopathology 111:2355 (doi:10.1094/phyto-07-20-0302-fi; PMID 33829853), not open access, so I could not read its methods. BioSample attribute `mating_type` is `1-2` or `1-1` with no method stated. Because the 32 labels match the sequence in 32/32 cases, the submitter very likely derived them from the genome sequence or a MAT PCR. The 15 sibling isolates are labelled `1-1` and are called MAT1-1 correctly. The raw value `1-2` is the standard short form of MAT1-2; the label does not mean something else.

## What goes wrong in the software

- Search: fine. All genes are found, at the right place.
- Polish: fine; models are correct for the references given. Exonerate on Rf1 with the B. cinerea references: MAT1-2-1 raw score 1914 (99.7%, full length) versus MAT1-1-1 1085 (92.7%, residues 129-348).
- Vote (`assign_idiomorph` / `idiomorph_candidates`): wrong outcome, because the inputs are unequal. The MAT1-1 side has a 100% reference, the MAT1-2 side only a 77% cross-species one, and a truncated gene is allowed to vote at full strength.
- Labelling: it faithfully reports the vote. The confidence label `high` is unwarranted, since the campaign's own report shows the margin came from one gene with partial reference coverage and the HMG gene that "lost" was polished and full length.

## Proposed changes (not applied)

Minimal fix, database only: curate a B. cinerea MAT1-2 record (e.g. `40559_sas405_MAT_MAT1-2`) from KF944386.1 (SAS405, MAT1-2-1 AHX22634.1 and MAT1-2-4 AHX22635.1, to be named MAT1-2-10 in the DB). Its own metadata says "a companion MAT1-2 record can be built from it". Exonerate scores on Rf1 then become roughly 1914 (MAT1-2-1) against 1085 (MAT1-1-1), i.e. about 730 versus 416 bits in the report's units, and the existing vote calls MAT1-2 with no code change. Also add the KF944386 MAT1-1-1 remnant to the MAT1-2 record's note, so that it is understood.

Robustness fix, `src/MATPredict/detect/idiomorph.py`: let a hit vote with `score * min(1.0, reference_coverage)` (or drop votes of core genes below about 70% reference coverage when the other idiomorph has a full-length core gene). For polished hits this needs the reference length, which `gene_evidence` does not record for the exonerate path (`coverage: null`); record aligned/reference residues there. Sketch:

```python
# idiomorph_candidates, inside the hit loop
cov = hit.reference_coverage  # aligned aa / reference aa; None when unknown
if cov is not None:
    score *= min(1.0, cov)
```

With it, Rf1 MAT1-1-1 would score 416 * 0.61 = 254 versus MAT1-2-1 222 at 77% / near 1.0 coverage; still close, so the DB fix is the primary one. Alternative or addition: use the alpha/HMG profiles in `db/Ascomycota/classifiers/MAT_core/` (branch ascomycota-core-hmms) as the tiebreaker, which gives HMG 221 vs alpha 125 here.

Test case to add: Rf1 (GCA_015147965.1_ASM1514796v1) and T4 expect MAT1-2; B05.10 and Sl12 (GCA_015147975.1) expect MAT1-1.

## Files and scratch

Scratch on the HPCC: `/bigdata/stajichlab/jstajich/projects/botrytis-mislabels/` (miniprot gffs in `mp/`, locus fasta `loci.fa`, reference genomes `ref/`, k-mer script and output `kmer.py`, `kmer.out`). Jobs: botrytis-miniprot, botrytis-kmer.

## Caveats

- Rf2 (GCA_015148645.1) and Sl8 (GCA_015148665.1) have APN2 and SLA2 on different contigs, so synteny cannot be shown there; Rf2's gene content is still MAT1-2-1 plus MAT1-2-4 (MAT1-1-1 remnant 33% coverage).
- Read check covers one isolate (Rf1); SRA metadata shows each isolate is a single HiSeq 4000 run from one in vitro culture, with no sign of pooled samples.
