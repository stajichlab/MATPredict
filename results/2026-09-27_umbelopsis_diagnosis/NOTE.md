# Umbelopsidaceae calls lost, and two unmodelled mislabels (2026-09-27)

Read-only diagnosis. Scans compared: results/2026-09-26_early_diverging/Mucoromycota
(code 634dda4) vs results/2026-09-27_btbA_homothallic/Mucoromycota_4457c3a.

## Lost calls: 5 genomes, one cause

Lost: U. nana GCA_054906655.1 (was Plus), U. sp. WA50703 GCA_964291815.1 and
GCA_053572175.1, U. vinacea GCA_977110975.1 and GCA_016758895.1 (all were Minus).
U. sp. AD052 GCA_025677805.1 was never called (0 modelled genes in both).

Every lost report says: "best cluster was flank-carried -- no core gene could be
modelled -- and its strongest core hit is weaker than E=1e-05 or lies outside the
flank span (+-20000 bp)". This is the flank-carried rule (3684c60, revised
b0898d9). In 634dda4 these were medium calls whose only modelled genes were
tptA/algA(/glrA); sexM/sexP were unpolished at 28-36%.

Direct tblastn of the 15 curated sexM/sexP proteins against a 100 kb region
(core_hits.txt): the core hit is an HMG-box fragment (59-87 aa, 28-36%, 33-44
bits) about 1.5 kb from tptA, with algA adjacent and glrA 10-15 kb away -- the
Mucorales arrangement, well inside the 20 kb window. Region E-values are
2.2e-4 (U. nana) and 6e-8 to 2e-6 (the others). Detect computes E against the
whole genome (~250x larger than the region), so ~40-bit hits fall above 1e-5.
The floor, not the window, withholds them.

No full model forms: miniprot built none of the six core genes. The genes are
too divergent from the curated set (no Umbelopsidales record exists).

## Labels the classifier gives (HSP translation vs sexM.hmm/sexP.hmm)

| locus | sexM | sexP | margin P-M |
|---|---|---|---|
| U. nana (was Plus) | 25.4 | 19.1 | -6.3 (below the 25-bit floor) |
| U. sp. WA50703 x2 | 32.4 | 12.8 | -19.6 |
| U. vinacea gz | 30.7 | 20.7 | -10.0 |
| U. vinacea WA | 38.2 | 22.9 | -15.3 |
| Umbelopsis sp. M5902 (labelled Minus) | 21.7 | 97.0 | +75.3 |
| Mucor hiemalis gzMucHiem1 (labelled Plus) | 97.2 | 50.9 | -46.3 |

M5902: best tblastn hit is sexP (58.2 bits, E 6.8e-12 over 131 aa) -- Plus.
M. hiemalis: best hits are sexM (53-55 bits over 158-194 aa) -- Minus. Both were
labelled by the fallback rule (short-hit identity) because no core protein was
modelled, so the classifier never ran.

## Recommended fixes (not implemented)

1. src/MATPredict/detect/flank_carried.py: the E<=1e-5 floor depends on genome
   size. Use a bitscore floor (or E scaled to a fixed search space), or accept a
   core hit within a few kb of a conserved flank (tptA) for Mucoromycota.
2. src/MATPredict/detect/classifier.py / pipeline idiomorph decision: when no
   core protein is modelled, score the tblastn HSP translation with the
   classifier instead of the identity fallback; keep min_margin.
3. Curate a tier-2 Umbelopsis record so its sexM/sexP genes can be modelled.
