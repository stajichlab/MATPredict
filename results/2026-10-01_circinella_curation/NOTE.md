# Circinella-group MAT curation (tier 2)
Status: pending curator sign-off

## Question
Curator ruling 2026-10-01: curate the Circinella group (Circinella,
Thamnostylum, Fennellomyces, Zychaea) with one tier-2 Minus and one tier-2
Plus record. Let the MAT-gene gate accept rnhA-only flank support for this
group only. Does the change call the group correctly against
`label_treatment.tsv`, and what else changes?

## Data and code version
- Baseline: polish-scope-cuts a8863f3 (frozen worktree run-a8863f3).
- Candidate: curation-circinella f9ad092 (frozen worktree run-f9ad092).
- Label truth: results/2026-10-01_circinella_label_tree/label_treatment.tsv.
- RAxML-NG check: results/2026-10-01_circinella_label_tree/rx_aln_hmm.clades.tsv.
- Sources: NCBI datasets GCF_025766255.1 and GCA_025093555.1 (genomic.gff,
  protein.faa); contigs fetched by efetch.
- Regression: results/2026-10-01_regression_circinella/ (jobs 29328317-27).
  LCG Mucoromycotina, 621 genomes: `lcg/` (jobs 29328334-38). Both sides use
  `--phylum Mucoromycota --taxid <genus taxid>`.

## Method
1. Checked the RAxML-NG tree against label_treatment.tsv.
2. Translated each annotated CDS from the contig and compared it to the
   deposited protein. Found PF00505 with pyhmmer. Scored with the classifier
   (`verify_sources.py`).
3. Wrote the records (`make_records.py`). Built locus files with
   `curate-db build-gff`.
4. Added a taxon-scoped `flank_support_groups` entry to the gate.
5. Rebuilt the classifier with `scripts/build_idiomorph_hmms.py`. Ran it twice.
6. Ran the regression panel, the LCG set, label scoring (`score_labels.py`)
   and the record self-check.

## Results
**RAxML-NG.** All 10 labelled strains sit in the clade that
label_treatment.tsv gives. MRCA support: sexP references 42, sexM references
24. There is no contradiction.

**Records** (reviewed_by: null):
- `64656_rsa-1403_MAT_Minus`: *Zychaea mexicana* RSA 1403, GCF_025766255.1
  (RefSeq, because BFD uses it; GenBank copies are identical). Contains
  sexM XP_052979473.1 (191 aa, HMG box aa 46-113, E 5e-26) and rnhA
  XP_052979475.1.
- `101103_nrrl1351_MAT_Plus`: *Circinella umbellata* NRRL 1351,
  GCA_025093555.1. Contains sexP KAI7847721.1 (310 aa, HMG box aa 117-177,
  E 1.5e-15) and rnhA KAI7847723.1.
- All 10 proteins are 100% identical to their exon translations (table 1,
  codon_start 1).
- VMA1, pntA and gpmI are conserved in both genera (77-92% identity). They
  are in `extended_flank` only. A first draft added them to the roster as
  unsearched genes. That changed the scoring denominator for every Mucorales
  genome (test_partial_locus failed), so I reverted it.

**Scope decision.** The group stays in Mucoromycota:MAT, with a taxon-scoped
gate entry and no sub-family. It uses the same classifier and roster. Only
the gate's flank count differs. The rule applies only when the genome's
taxid lineage holds one of the four genera. Without a taxid it never applies.

**Classifier.** LOO 88/88 → 90/90. Worst correct margin 13.6 → 12.6 bits
(*U. vinacea* sexM). Gate 98.4 → 99.7 bits. The rerun gave byte-identical
HMMs. Circinella sexP-type proteins rise from 111-116 to 163-176 bits. These
are in-sample scores. LOO scores: *C. umbellata* 115.8, *Zychaea* 132.9.

**Tests.** 972 passed.

**Regression panel.** Ascomycota, Basidiomycota and record_sources: 0
changes. No call was lost and no confidence changed.
- Zygo23: 24 call-touching changes, all classifier shifts or model changes.
  No label changed.
- Mucoromycota: 28 call-touching changes. Three are gains:
  - *C. umbellata* Cirumb1: Plus.
  - *C. minor* GCA_016758965.1: Plus.
  - *Phascolomyces articulosus*: Plus at 116.9 bits. It is outside the
    group and is gained on score.
- Three *Umbelopsis* calls gain an unpolished rnhA hit at 30-33% identity,
  and their spans grow. Their labels do not change.
- The group rule kept 0 panel calls.

**LCG.** 418 call-touching changes:
- 0 calls lost.
- 19 gains:
  - 14 Circinella Plus calls, all on score.
  - *R. microsporus* NRRL A-17693 (C7): Plus, on score.
  - *Absidia* sp. NRRL 3163: Plus, on score.
  - Fennellomyces RSA 2442 and RSA 2446: undetermined, kept by the group
    rule.
  - Thamnostylum RSA 1405: Plus, kept by the group rule.
- 3 changes from undetermined to Minus: *M. hiemalis* ×2 and
  *Syncephalastrum* sp. NRRL 1507.
- The group rule kept only 3 calls, all in the group.
- New observations:
  - *C. minor* CBS 143.56 (LCG) now has both a Plus and a Minus call
    (C5 updated).
  - *C. muscae* NRRL 1355, 1363 and 2403 each have two sexP+rnhA Plus calls
    on separate scaffolds (164.9 and 169.0 bits).

**Label score** (8 non-provisional strains):
- Baseline: 6 correct, 2 uncalled.
- Candidate: 6 correct, 1 undetermined (RSA 2442), 1 uncalled (RSA 1415:
  sexP withheld by the modelled-gene bar).
- Both provisional Thamnostylum strains are now called Plus. They are not
  scored.

**Self-check.** Both new records call their own source genome: Minus/medium
and Plus/medium.

## What changed in detection
- Two new references.
- The rebuilt classifier.
- Gate `flank_support_groups`, with `mat_gene_gate_group` in the report.

## Limits
- The Plus assignment of the Circinella sexP clade rests on the classifier
  and two Fennellomyces labels. The HMG-only tree support is 0/19.
- Detect still models the Fennellomyces sexP badly (44-96 bits). No public
  Fennellomyces Plus genome exists.
- The group rule only makes 2 undetermined calls among the scored strains.

## Curator decisions
1. Sign off both records.
2. Keep or drop the group rule. It adds little beyond score.
3. The *Phascolomyces* and *Absidia* NRRL 3163 gains.
4. *C. minor* CBS 143.56 and *C. muscae* double calls.

## Files
`verify_sources.{py,tsv}`, `make_records.py`, `rescore.{py,tsv}`,
`manifest_before.yaml`, `build.log`, `hmm_md5_build1.txt`,
`score_labels.py`, `label_score.tsv`, `lcg/`.
