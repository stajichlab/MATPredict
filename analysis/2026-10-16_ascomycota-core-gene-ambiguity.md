# Ascomycota MAT1-1-1 and MAT1-2-1 ambiguity (2026-10-16)

Data: `results/2026-10-03_ascomycota_v060` (v0.6.0), 12,628 `Ascomycota:MAT` loci. Per-locus gene evidence
from the 19,380 reports in `reports_all.tar.zst`. Code and lists: `results/2026-10-16_ascomycota_core_ambiguity/`
(`core_diag.py`, `core_tsv.py`, `ambiguous_loci.tsv`, `label_disagrees.tsv`). Descriptive; no ground truth.

Question: where the locus carries both core genes, is the idiomorph label well supported, and would a profile
HMM for the alpha box (MAT1-1-1) and the HMG box (MAT1-2-1) help?

## Findings
- 4,883 loci (39%) have both MAT1-1-1 and MAT1-2-1 in the gene evidence; only 201 are flagged
  `homothallic_candidate`.
- Most have a clear winner: alpha modelled only 2,241 (median gap 180 bits; the label follows the higher
  bitscore in 2,224); HMG modelled only 458 (453); both modelled 2,028 (1,992, of which 175 have a gap under
  20 bits); neither modelled 156 (the label is the higher-bitscore gene in only 10).
- 204 loci have a label that is not the higher-bitscore gene: 191 MAT1-1 where the HMG gene scores higher, and 13
  undetermined. None are MAT1-2 against the alpha box. Sordariomycetes 104, Dothideomycetes 99.
- Three patterns in the 204:
  - 87 Pleosporales (alpha hit 32.2% identity, 31 bits; HMG hit 40%, 55 bits; genes found include MAT1-1-3 and
    MAT1-1-4). This is the known *Parastagonospora* HMG-box-named-MAT1-1-3 defect (see
    `2026-10-05_dothideomycetes-full-run.md`: 178 MAT1-1 in v0.6.0, 91 and 91 with records). Not re-checked
    against the records run.
  - 58 exact bitscore ties (gap 0), 45 of them Sordariomycetes tblastn-only hits (31.6% and 44.1% identity,
    55.8 bits each): the label looks like a tie-break to MAT1-1.
  - 59 with gaps of 6 to 40 bits, mostly Sordariomycetes.
- `ambiguous_loci.tsv` (341): both modelled with gap under 20 bits 170; neither modelled with label
  disagreeing 141; both weak (identity under the 5th percentile of single-gene loci, alpha 31.6, HMG 35.3) 20;
  ten loci fall in two groups.

## Caveats
- A small gap or a disagreement means uncertain, not wrong.
- This run counts every gene-evidence entry for the two genes (4,883 loci); an earlier count from the
  `genes_found` column gave 3,333.
- The identity threshold for "weak" is lenient (5th percentile of single-gene loci).

## Open
- Whether to build alpha-box and HMG-box profiles (about 117 loci would benefit: the 58 ties and 59 small
  gaps), or fix the tie-break first (a tie is not evidence for MAT1-1).
- Whether the 87 *Parastagonospora* loci are already corrected in the records run.
