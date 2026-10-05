# 025. Mycotypha africana: sexM 153 kb from sexP, encoded by an unusually long 531-aa ORF

- **Category:** curation-first; biology; homothallism
- **Status:** candidate (record awaiting curator sign-off; gene model not checked with RNA)
- **Lineage:** Mucorales, Mycotyphaceae

## Summary
In the homothallic Mycotypha africana NRRL 2978 the sexM gene lies 153 kb from
the sexP locus on the same scaffold, as Schulz et al. 2016 reported (~150 kb).
MATPredict did not call it: it sits alone, with no MAT flank gene within 60 kb,
and it had been placed in the classifier's non-MAT negative set. Its open reading
frame is 531 aa, more than twice the length of other curated sexM proteins
(150-249 aa in the 11 other records): an N-terminal HMG box followed by a long C-terminal region
with two weak HMG-like repeats. The sexM-like part (aa 26-218) aligns along the
full length of the other sexM proteins, including the C-terminal conserved
motifs.

## Evidence
- Record `64632_nrrl-2978_MAT_combined` (db/Mucoromycota/Mucorales/), GCF_025528875.1
  (JGI Mycafr1), NW_026515730.1: sexP XP_052935775.1 at 1,357,729-1,358,472 (+)
  with rnhA XP_052935776.1 and glrA XP_052935777.1; sexM 1,203,454-1,205,049 (-),
  single exon, Met to stop, 531 aa.
- PF00505 HMG box aa 33-98 (i-E 6.8e-16); HMG-like repeats aa 221-249 and
  367-395 (i-E 5.4e-4 each).
- Classifier (2026-10-01 build) on the 531-aa ORF: sexM 171.6, sexP 70.4, P1 34.6
  bits; leave-one-genus-out in the rebuild: sexM margin +50.6 (aa 26-218).
- Gene tree: the protein groups with sexM in the 2026-10-03 full-length tree
  (analysis/2026-10-04_unclassified-hmg-in-gene-trees.md).
- Annotations disagree: RefSeq has only an 83-aa partial HMG model
  (XP_052935737.1, 1,204,699-1,204,947); the BFD re-annotation
  (F5387F0C_005515) spans the same ORF as this record.
- Alignment: results/2026-10-04_mycotypha_alignment/ (sexM_training.png,
  sexM_fulllength.png).
- Literature: Schulz E, Wetzel J, Burmester A, Ellenberger S, Siegmund L,
  Wöstemeyer J 2016. Endocytobiosis Cell Res 27(4):39-57, p. 46 (no DOI/PMID).

## Method that found it
Placement of non-MAT HMG genes in the 2026-10-03 sexP/sexM gene trees;
tblastn/exonerate of curated sexM against the scaffold; ORF extension.

## Limits
- The 531-aa length rests on one ORF in one assembly; no transcript evidence.
  The C-terminal repeats could be real or an assembly or annotation artefact.
- The classifier trains on aa 26-218 only, so the long tail does not shape the
  sexM HMM (curator ruling 2026-10-04).
