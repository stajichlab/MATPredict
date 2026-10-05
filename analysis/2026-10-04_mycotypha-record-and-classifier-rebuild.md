# Mycotypha africana combined record, cleaned negatives, classifier rebuild

Status: awaiting curator sign-off (record and rebuild).

## Question
Curate the Mycotypha africana sexM found in the gene trees
([non-MAT HMG report](2026-10-04_unclassified-hmg-in-gene-trees.md)), remove
probable MAT genes from the classifier negative set, rebuild the classifier,
and measure what changes.

## Ruling
J. Stajich, 2026-10-04: curate + clean + rebuild, then sign off; one combined
homothallic record; record the full sexM ORF but train on the HMG region;
remove the duplicate Mycotypha sexP from `training_extra.faa`.

## Record: `64632_nrrl-2978_MAT_combined` (Mucorales)
- Mycotypha africana NRRL 2978, GCF_025528875.1 (JGI Mycafr1), homothallic.
  Literature: Schulz et al. 2016, Endocytobiosis Cell Res 27(4):39-57, p. 46
  (sexP and sexM on one scaffold, about 150 kb apart; no DOI or PMID).
- Segment 0, NW_026515730.1:1,357,729-1,365,893: sexP (XP_052935775.1, 247 aa),
  rnhA (XP_052935776.1), glrA (XP_052935777.1). Validation: 100% identity and
  coverage for all three.
- Segment 1, NW_026515730.1:1,203,454-1,205,049 (-), 153 kb from sexP: sexM,
  a single-exon ORF, Met to stop, 531 aa, curator-derived. RefSeq has only an
  83-aa partial model here (XP_052935737.1). PF00505 at aa 33-98 (i-E 6.8e-16);
  two weak HMG-like repeats (aa 221-249, 367-395). Classifier (2026-10-01
  build) on the full ORF: sexM 171.6, sexP 70.4, P1 34.6 bits. No flank gene
  within 60 kb.
- New gene field `classifier_training_region` (schema, `classifier_build`,
  2 tests): the sexM trains on aa 26-218 only; the record and search keep the
  full ORF.

## Negative set
`db/Mucoromycota/classifiers/MAT/paralog_negatives_excluded.tsv` (read by
`scripts/make_paralog_negatives.py`; without it the old file is reproduced byte
for byte). Rule: inside a supported sexP or sexM clade of the 2026-10-03 gene
trees and classifier margin >= 25 bits. 189 -> 186:

| Protein | Scores (sexP / sexM) | Tree |
|---|---|---|
| Mycafr1 h1 (now the record's sexM) | 64.9 / 116.3 | sexM clade |
| Dicele1 h5 (Dichotomocladium elegans; genome uncalled) | 187.3 / 41.7 | Circinella-type sexP clade, UFBoot 96 |
| C. minor h6 (GCA_016758965.1; genome called Plus) | 157.2 / 31.6 | inside the sexP MRCA |

## Training set
`training_extra.faa`: the Mycotypha sexP entry (294 aa = the RefSeq 247 aa plus
a 47-aa N-terminal extension with no Met from the older MATPredict model) is
removed; the record's RefSeq protein replaces it.

## Rebuild
Locked environment (MAFFT 7.526, pyhmmer 0.12.3). Control: a rebuild of main's
database with this environment reproduced the committed sexP, sexM and P1 HMMs
byte for byte, so all changes below come from this work.

| | 2026-10-01 build | This build |
|---|---|---|
| Training proteins sexP / sexM | 79 / 11 | 79 / 12 |
| Leave-one-genus-out | 90/90, worst margin 12.6 | 91/91, worst margin 14.2 |
| Mycotypha margins (LOO) | - | sexP +151.4, sexM +50.6 |
| MAT-gene gate (95th pct of negatives) | 99.7 bits (189) | 97.8 bits (186) |

## Effect on calls
Regression panel + Zygo 23 vs v0.6.1 (`results/2026-10-04_regression_mycotypha/`):
- Ascomycota, Basidiomycota, record-source genomes: 0 changed loci.
- Zygo 23: 23/23 on scaffold and contig inputs; 6 loci with score shifts only.
- Mucoromycota panel: no idiomorph flip; 2 calls gained (Rhizomucor pusillus
  FCH_5_7 and CBS 183.67 Rhipu1, Plus, medium; withheld by the modelled-gene
  bar before); Mycotypha sexP now modelled full length (741 bp, 100%); the rest
  score shifts.
- Record self-calls: 34 called (Mycotypha record: Plus/high; its sexM is not
  called because a lone gene without flanks is not a locus), 5 missed as before.

All 288 BFD Mucoromycotina genomes, campaign mode, vs the 2026-10-03 campaign
(`results/2026-10-04_mucoro_bfd_rebuild/compare_vs_campaign.txt`):
- 286 same; 2 gained (the two R. pusillus above); 1 changed.
- Mucor griseocyanus CBS 116.08: Minus -> undetermined. Its sexM model is a
  ~112-aa HMG-box fragment (77.7% identity to the R7B sexM); the classifier
  scored it weakly before (Minus 35.2 vs Plus 4.1, margin 31.1, just above the
  25-bit minimum) and now 20.1 vs 3.8 (margin 16.3). A fragile call on a short
  model, not a reversal; identity still points to Minus.
  Cause (checked 2026-10-04 at the curator's request): the assembly
  (GCA_060309175.1, Oxford Nanopore MinION only, SRR32276708) has a 2-nt
  frameshift in an A5 homopolymer inside the sexM, after "...PKPSR". The curated
  R7B sexM aligns over aa 1-247 of 249 at 77.7% with that one frameshift; the
  model's translation stops there. The frame-corrected protein (253 aa) scores
  sexM 137.9 vs sexP 44.3 in the rebuild (156.6 vs 40.9 before). So this is a
  genome (assembly) problem, not a database or classifier problem. Notable
  finding 026.

## Curator comments (2026-10-04)
- Record: seems okay; the curator asked to see the multiple alignment that
  includes the Mycotypha record (`results/2026-10-04_mycotypha_alignment/`:
  `sexM_training.png` = the 12 training sexM as built, Mycotypha aa 26-218;
  `sexM_fulllength.png` = with the full 531-aa ORF) and a manuscript note on the
  long protein (notable finding 025).
- M. griseocyanus: find why the model is truncated (done: Nanopore frameshift;
  notable finding 026).
- Training region: agreed to use only the part of the Mycotypha sexM that
  matches the other sexM proteins up to their C-terminal end (aa 26-218 does
  this). The curator asked for further work on refining these alignments and
  testing for better match and placement (plan under Open).
- R. pusillus new Plus calls: acceptable; uncurated discoveries are part of the
  approach, not everything will be a curated set.

## Decision
Pending curator sign-off after the alignment review.

## Open
- Mycotypha sexM is still not called in its genome (no flanks); a report-only
  check for an unlinked second idiomorph would surface it.
- Dichotomocladium elegans h5: candidate missed sexP (review).
- Frameshift-aware classification (proposal): when exonerate reports a
  frameshift in a core MAT gene model, classify the frame-corrected translation
  and flag `frameshift_in_model`. Would restore M. griseocyanus (Minus) and help
  Nanopore-only assemblies. Not implemented; needs curator approval.
- Rhizomucor pusillus: no curated Rhizomucor record; the two new Plus calls
  are medium confidence.

## Plan: alignment refinement and placement tests (curator request 2026-10-04)
1. Alignment method: compare MAFFT L-INS-i (current), MAFFT E-INS-i, and
   hmmalign to the sexM/sexP HMMs, scored by LOO margins and column agreement
   in the HMG box and the C-terminal motif block.
2. Region boundary: rebuild with the Mycotypha sexM region ending at aa 200,
   218 (current) and 240, and starting at aa 1 vs 26; measure LOO margins, the
   gate and the regression panel. Keep the boundary that changes calls least
   and keeps the C-terminal motifs.
3. Placement: place the Mycotypha sexM region (and other long or partial
   models) on the full-length reference tree with EPA-ng and on a fresh
   IQ-TREE tree; check that it falls in sexM with support.
4. Apply the same region rule to any future record whose protein is much
   longer than its family's training set (flag by length outliers in the
   build manifest).

## Files
Record `db/Mucoromycota/Mucorales/64632_nrrl-2978_MAT_combined/`; classifier
`db/Mucoromycota/classifiers/MAT/`; `results/2026-10-04_regression_mycotypha/`;
`results/2026-10-04_mucoro_bfd_rebuild/` (reports archive, comparison).
