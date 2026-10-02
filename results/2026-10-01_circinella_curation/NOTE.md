# Circinella-group MAT curation (tier 2)
Status: pending curator sign-off (records); group rule dropped

## Question
Curator ruling 2026-10-01: curate the Circinella group (Circinella,
Thamnostylum, Fennellomyces, Zychaea) with one tier-2 Minus and one tier-2
Plus record. Follow-up rulings, same day:
1. Before sign-off, check that the Plus record's sexP sits in the sexP clade
   and the Minus record's sexM in the sexM clade, in both alignments and
   with IQ-TREE and RAxML-NG.
2. Drop the rnhA-only gate rule (`flank_support_groups`). Keep the records
   and the rebuilt classifier.
3. Examine the two out-of-group Plus gains (*Phascolomyces articulosus*,
   *Absidia* sp. NRRL 3163) and the 14 Circinella-group gains.

## Data and code version
- Baseline: polish-scope-cuts a8863f3 (frozen worktree run-a8863f3; equal to
  the remote branch head on 2026-10-01).
- Candidate v1 (records + classifier + group rule): curation-circinella
  f9ad092 (run-f9ad092).
- Candidate v2 (group rule removed): curation-circinella 94d9d37
  (run-94d9d37).
- Label truth: results/2026-10-01_circinella_label_tree/label_treatment.tsv.
- Regression v1: results/2026-10-01_regression_circinella/ (jobs
  29328317-27). Regression v2: results/2026-10-01_regression_circinella_v2/
  (candidate-only jobs 29335322-26, diffed against the v1 baseline outputs
  of the same a8863f3 tree).
- LCG Mucoromycotina, 621 genomes: `lcg/base`, `lcg/cand` (v1), `lcg/cand2`
  (v2, jobs 29335320-21). All runs use `--phylum Mucoromycota --taxid
  <genus taxid>`.
- Trees: `placement/ml/` (job 29335283), `placement/ft/`.

## Method
1. Records: translated each annotated CDS from the contig, compared it to
   the deposited protein, found PF00505 with pyhmmer, scored with the
   classifier (`verify_sources.py`). Wrote records (`make_records.py`), built
   locus files with `curate-db build-gff`.
2. Rebuilt the classifier with `scripts/build_idiomorph_hmms.py`, twice.
3. Placement (`placement/`): label-tree protein set (105 tips) plus detect's
   scored sexP model for the two out-of-group gains (`gains_extract.py`) =
   107 tips. hmmalign to PF00505 (69 match columns) and MAFFT L-INS-i full
   length, columns with >= 50% gaps removed (212 columns)
   (`build_aln.sh`). IQ-TREE 3.0.1 (ModelFinder, 1000 UFBoot, 1000 SH-aLRT)
   and RAxML-NG 2.0.2 (`--all`, 200 bootstraps, ModelFinder model) on both
   alignments. Rooted on *F. graminearum* MAT1-2-1. `mrca.py` tests whether
   a tip is inside the MRCA of all known Plus (or Minus) tips and whether
   that clade excludes the other type; `place.py` gives the first clade
   with a known-type tip; `group.py` gives the gains' clade with Circinella
   sexP.
4. Four gained sexP copies have no tip in the ML set (*C. minor* CBS 143.56
   Plus locus; the second sexP copy in *C. muscae* NRRL 1355/1363/2403).
   Added them to FastTree trees only (`second_copies.py`, `ft/`).
5. Removed the group rule: restored `mat_gene_gate.py`, `pipeline.py`,
   `report.py`, `order.schema.yaml` and `Mucoromycota/order.yml` to a8863f3;
   deleted `tests/detect/test_mat_gene_gate_groups.py`.
6. Reran the regression panel and LCG on v2. Diffed v2 against the baseline
   and against v1. Rescored labels (`label_score_v2.tsv`).

## Results
**Records** (reviewed_by: null):
- `64656_rsa-1403_MAT_Minus`: *Zychaea mexicana* RSA 1403,
  GCF_025766255.1. sexM XP_052979473.1 (191 aa, HMG box aa 46-113) and rnhA
  XP_052979475.1.
- `101103_nrrl1351_MAT_Plus`: *Circinella umbellata* NRRL 1351,
  GCA_025093555.1. sexP KAI7847721.1 (310 aa, HMG box aa 117-177) and rnhA
  KAI7847723.1.
- All proteins are 100% identical to their exon translations (table 1).
  VMA1, pntA and gpmI are in `extended_flank` only (not searched).

**Classifier.** LOO 88/88 → 90/90. Worst correct margin 13.6 → 12.6 bits.
Gate 98.4 → 99.7 bits. Rebuild is byte-identical. Leave-one-out scores:
*C. umbellata* sexP 115.8 vs sexM 32.0 (margin 83.8), P1 0.0; *Zychaea*
sexM 132.9 vs sexP 70.7 (margin 62.2), P1 41.2.

**Placement, 107-tip trees.** "In clade" = inside the MRCA of all
known-type tips of that type, with 0 tips of the other type inside.

| Tree | Plus record in sexP clade (support) | Minus record in sexM clade (support) |
|---|---|---|
| IQ-TREE, HMG box (Q.INSECT+I+G4) | yes (98.2/87) | yes (89.4/80) |
| IQ-TREE, full length (Q.MAMMAL+F+I+R4) | yes (98.5/99) | yes (95.2/70) |
| RAxML-NG, HMG box | yes (42) | yes (19) |
| RAxML-NG, full length | yes (83) | yes (35) |

- Neither record falls in the other type's clade in any tree.
- With the tier-2 genome-derived references (13706, 41833, 44442) and the
  labelled strains removed from the clade definition, the Plus record is
  still inside the sexP clade in both full-length trees (IQ 98.5/99, RAxML
  83). It is outside in both HMG-box trees, where the established-reference
  sexP MRCA has only 37/65 (IQ) and 9 (RAxML).
- The nearest known-type tips of the Plus record are the tier-2
  *S. racemosum* NRRL 2496 Plus record and the labelled *S. racemosum*
  CBS 440.59 (IQ full length 87.9/93; IQ HMG box 0/39).
- The Minus record's nearest known-type tip is the *Phycomyces* NRRL 1555
  sexM reference (IQ HMG box 88.5/69; IQ full length 94.7/94; RAxML 25 and
  47).
- The earlier 105-tip trees agree (IQ HMG box: Plus 98.5/86, Minus
  87.3/83; RAxML HMG box: 42, 24).

**Group rule removed (94d9d37).** Tests: 962 passed (972 minus the 10
group-rule tests).

**Regression panel, v2.**
- Against v1: 0 changed loci in every group.
- Against the baseline: Ascomycota, Basidiomycota and record_sources 0.
  Mucoromycota 28 call-touching loci, including 3 gains (*C. umbellata*
  Cirumb1, *C. minor* GCA_016758965.1, *Phascolomyces*). Zygo23: 24
  call-touching loci (classifier shifts and model changes), no label
  changed. These are the same changes as v1.

**LCG, v2.**
- Against v1: 3 calls lost, all of them calls the group rule kept:
  *F. verticellatus* RSA 2442 and RSA 2446 (undetermined) and
  *T. nigricans* RSA 1405 (Plus). Nothing else changed.
- Against the baseline: 415 call-touching loci, 0 calls lost, 16 gains,
  3 undetermined → Minus (*M. hiemalis* ×2, *Syncephalastrum* sp.
  NRRL 1507).
- The 16 gains:
  - 13 Circinella loci, *T. lucknowense* RSA 1093- and *R. microsporus*
    NRRL A-17693. All are sexP + rnhA, all reference
    `101103_nrrl1351_MAT_Plus`, and all pass the gate on score
    (163-176 bits; *T. lucknowense* 99.7, exactly at the gate).
  - *Absidia* sp. NRRL 3163 (see below).
- The *R. microsporus* NRRL A-17693 locus is nucleotide-identical
  (7,705 bp) to *C. umbellata* NRRL 2417 and *C. minor* NRRL 1365. The
  label tree already lists this genome as = *C. minor* by phylogeny.

**Label score, v2** (8 non-provisional strains): 6 correct, 2 uncalled
(RSA 1415, RSA 2442). This is the same as the baseline. RSA 2442 was
undetermined in v1.

**Circinella-group gains vs tree.** 15 loci: the 14 counted in v1 (13
Circinella + *T. lucknowense*) plus *R. microsporus* A-17693 (= *C. minor*).
- 11 gained loci have their own genome's sexP as a tree tip. All 11 tips
  are in the sexP half (IQ HMG box 0/19 for their clade, IQ full length
  83/93 in the label tree).
- The 4 gained copies without a tip were placed by FastTree. All 4 are
  inside the sexP MRCA (HMG box 0.99, full length 0.897).
- No gained call disagrees with the tree.

**Out-of-group gains.**

| | *Phascolomyces articulosus* RSA 2281 (GCA_025677815.1) | *Absidia* sp. NRRL 3163 (LCG) |
|---|---|---|
| Locus | JAIXMP010000004.1:1,101,267-1,109,036 | scaffold_198:27,667-35,354 |
| Genes | sexP(+), rnhA(+) 2.0 kb downstream; tptA, algA, glrA absent | sexP(−), rnhA(−) 2.2 kb away; tptA, algA, glrA absent |
| Scored sexP model | miniprot, 341 aa, single exon, Met start | miniprot, 267 aa, Met start |
| Polished model | exonerate, 70 aa, 62.9% to the Plus record | exonerate, 98 aa, 60.4% |
| HMG box (PF00505) | E 1.4e-13, HMM 2-54 (C-terminus missing) | E 6e-14, HMM 2-56 |
| rnhA | 85.2% to the Plus record | 82.2% |
| Baseline call | withheld by gate: 48.5 bits (margin 22.3), rnhA only flank | withheld by gate: 41.7 bits (margin 18.9) |
| Same protein, baseline classifier | sexP 87.6, sexM 25.0 (below 98.4) | sexP 92.8, sexM 24.0 |
| v2 call | Plus, medium, partial_locus; sexP 116.9, sexM 24.2, margin 92.7 | Plus, medium; 116.4 vs 23.4, margin 93.0 |
| Gate route | score ≥ 99.7 (not flank) | score ≥ 99.7 |
| P1 | 9.0 (locus); 0.0 (model) | 4.2; 0.0 |
| Other calls in genome | none. Withheld: P1-class sexM-like locus (contig 14, P1 108.1) and HSP fragments | none. Withheld: P1-class locus (scaffold_430, P1 93.6) and fragments |
| Tree: in sexP MRCA | yes in all 4 ML trees (98.2/87, 98.5/99, 42, 83) | yes in all 4 |
| Tree: clade with Circinella sexP | 15 tips: IQ HMG 81.8/62, IQ full 97.1/100, RAxML full 91; RAxML HMG 73 (20 tips) | same clade |

- The two gain proteins are 81.4% identical to each other (263 aa) and
  50-53% identical to the Circinella sexP. Their rnhA proteins are 86.9%
  identical.
- Two things drive both gains:
  1. The new Plus record gives a better sexP model. The baseline model
     scored 41.7-48.5 bits.
  2. The rebuilt classifier raises the same protein from 88-93 to
     116-117 bits.

**Verdict on the out-of-group gains: real Plus loci, not paralogs.**
- Each is a full-length sexP-type gene 2 kb from rnhA, the same layout as
  the Circinella group.
- Each sits inside the sexP clade in every tree, with the Circinella sexP
  (full-length support 91-100).
- Each has a P1 score below 10 and is the only call in its genome.
- *Absidia* sp. NRRL 3163 may be misidentified. Its sexP and rnhA are
  closest to *Phascolomyces*. This was not tested further.

## What changed in detection (v2 vs baseline)
- Two new tier-2 references.
- The rebuilt classifier (gate 99.7 bits).
- No gate code change. The group rule is gone.

## Limits
- The Plus clade assignment of the Circinella sexP-type rests on full-length
  trees and on tier-2 / labelled *Syncephalastrum* neighbours. HMG-box-only
  support for the immediate clade is weak (IQ 0/39; RAxML 8).
- RAxML-NG bootstrap support for the sexM and sexP MRCAs is low on the
  HMG box (19-42), as in the label tree.
- *T. lucknowense* RSA 1093- passes the gate at exactly 99.7 bits.
- The v2 panel reused the v1 baseline outputs (same commit, a8863f3) rather
  than rerunning them.

## Curator decisions
1. Sign off `101103_nrrl1351_MAT_Plus` and `64656_rsa-1403_MAT_Minus`.
2. Accept the *Phascolomyces* and *Absidia* NRRL 3163 Plus calls.
3. *Absidia* sp. NRRL 3163 identity (possible *Phascolomyces* relative).
4. *R. microsporus* NRRL A-17693 = *C. minor*: already noted in the label
   tree; decide whether to log it in ANNOTATION_ERRORS_FIXED_REPORT.md.

## Files
`verify_sources.{py,tsv}`, `make_records.py`, `rescore.{py,tsv}`,
`manifest_before.yaml`, `build.log`, `hmm_md5_build1.txt`,
`score_labels.py`, `label_score.tsv` (v1), `label_score_v2.tsv`,
`placement/` (`gains_extract.py`, `gains.{tsv,faa}`, `second_copies.{py,faa}`,
`build_aln.sh`, `tree_set_plus.faa`, `aln_hmm.afa`, `aln_mafft_g50.afa`,
`place.py`, `mrca.py`, `group.py`, `ml/`, `ft/`), `lcg/diff_v2_vs_base/`,
`lcg/diff_v2_vs_prev/`, `regression_v2/` (panel diff summaries).
