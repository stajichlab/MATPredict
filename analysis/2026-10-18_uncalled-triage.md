# Triage of uncalled Ascomycota genomes in ten lineage-routed orders (2026-10-18)

Code: `scripts/triage_uncalled.py`. Rows and summary: `results/2026-10-18_uncalled_triage/`. Question: for genomes the
detector did not call, is there a MAT-like gene, and is it next to the flank genes (SLA2, APN2, COX13)? The search is
independent of the pipeline: all six frames of each contig are searched with the alpha-box (MAT1-1-1) and HMG-box
(MAT1-2-1) profiles (db/Ascomycota/classifiers/MAT_core), and miniprot places the flank genes from the curated proteins.
Classes: C3 a MAT-like domain (own-profile score of at least 40 bits) within 30 kb of a flank gene on the same contig;
C2 MAT-like domain present but no flank gene within 30 kb; C1 no MAT-like domain, flanks found; C0 no flank gene found.
"Strong" means a score of at least 100 bits. 754 uncalled genomes (v0.6.0 status `uncalled`) and 120 called controls
(12 random per order, seed 1). Thresholds are provisional. Descriptive; the profiles come from 22 and 26 curated sequences.

## Controls (called genomes, should be C3)
Botryosphaeriales 12 of 12, Hypocreales 12, Magnaporthales 12, Chaetothyriales 11, Microascales 11, Mycosphaerellales 11,
Pleosporales 11, Amphisphaeriales 8: the triage reproduces the pipeline's calls where the profiles apply. It does **not** in
two orders: Pezizales 0 of 12 (called by the MATtub family, whose genes the Pezizomycotina MAT profiles do not fit, so
the triage is uninformative there) and Xylariales 1 of 12 (see below).

## Uncalled genomes
| order | uncalled | C3 linked (strong) | C2 unlinked | C1 none | C0 no flank |
|---|---|---|---|---|---|
| Xylariales | 184 | 0 (0) | 182 | 0 | 2 |
| Pleosporales | 150 | 11 (9) | 137 | 0 | 2 |
| Hypocreales | 79 | 5 (4) | 39 | 3 | 32 |
| Chaetothyriales | 73 | 2 (2) | 64 | 1 | 6 |
| Mycosphaerellales | 62 | 16 (1) | 42 | 0 | 4 |
| Botryosphaeriales | 49 | 0 (0) | 47 | 1 | 1 |
| Microascales | 47 | 19 (17) | 27 | 1 | 0 |
| Amphisphaeriales | 45 | 6 (0) | 39 | 0 | 0 |
| Pezizales | 35 | 1 (1) | 28 | 4 | 2 |
| Magnaporthales | 30 | 3 (3) | 26 | 0 | 1 |
| all | 754 | 63 (37) | | | |

## Reading
- **Most uncalled genomes have no MAT-like gene next to the flanks.** 63 of 754 (8%) do, and 37 of those are strong
  (at least 100 bits). Those 37 are the candidate pipeline misses: Microascales 17, Pleosporales 9, Hypocreales 4,
  Magnaporthales 3, Chaetothyriales 2, Pezizales 1, Mycosphaerellales 1. The remaining 691 either lack a linked MAT-like
  gene or have flanks the search could not place.
- **Xylariales supports the known absence.** 0 of 184 uncalled Xylariales genomes carry a MAT-like domain next to the
  flanks (182 have an HMG-like domain elsewhere, as expected from non-MAT HMG proteins). This agrees with Robinson and
  Natvig 2019 and with the NC1011 interval work. The controls add a caution: 11 of 12 called Xylariales genomes also show
  no linked MAT-like gene, so some existing Xylariales calls may rest on a non-MAT HMG protein. 12 genomes is a small
  sample and the calls were not inspected.
- **Microascales is the clearest recall target** (19 linked, 17 strong), consistent with the Ceratocystis finding that the
  references for this order are thin.
- **Where flanks are not found** (Hypocreales 32, Chaetothyriales 6) the class is C0: the flank search, not the
  MAT search, is the limit.

## Caveats
- A 40-bit floor is permissive, so C2 includes ordinary HMG-box proteins; C3 requires linkage, which is the stronger
  signal. A divergent MAT gene below 40 bits would be missed.
- Linkage depends on miniprot placing a flank at 30% identity and 40% coverage; a missed flank turns C3 into C2.
- The profiles are Pezizomycotina gene-identity profiles and are not informative for MATtub (Pezizales).
- The 37 strong candidates have not been inspected individually.
