# Genera that need a curator ruling on MATsc against MTL (draft, 2026-10-29)

From the 354-genome replay with the genome-level rule (`genus_table.tsv`, one row per genus: genomes, species, loci, family kept, mean best core-gene identity per family, gap, loci with MTL flank support). 156 `family_ambiguous` loci, 140 genomes, 27 genera. Identity gap = MATsc minus MTL, median over the genus's loci, in percentage points; the scores are best-hit identities of one gene (see README: about 6.5 points of noise per locus), so a gap under about 5 points is weak.

Nothing here is a ruling. The choice per genus is one of: **leave** (the genome-level identity rule decides, with `family_ambiguous`), **scope** (the genus is added to the `taxonomic_scope` of a family, as for MATkp, MATog and MATpk; needs at least one reference record from the genus or evidence that the references transfer), or **hold** (leave unresolved and report both).

## A. Consistent in the data, large gap (MATsc side): leave

| Genus | genomes (species) | kept | gap | note |
|---|---|---|---|---|
| Williopsis | 9 (6) | MATsc 9 | +9.7 | |
| Taiozyma | 3 (3) | MATsc 3 | +12.8 | |
| Millerago | 2 (2) | MATsc 2 | +15.7 | |
| Sporopachydermia | 3 (2) | MATsc 3 | +8.3 | |
| Kregervanrija | 4 (3) | MATsc 4 | +7.3 | one near tie called different idiomorphs under the old rule |
| Hansenula | 6 (6) | MATsc 6 | +3.4 | gap small |
| Candida (sensu lato, 3 species) | 3 (3) | MATsc 3 | +5.5 | |
| *Pichia* other than *P. kudriavzevii* (11 species) | 27 (11) | MATsc 28 | +4.3 | not routed to MATpk (species scope); MATpk call at 31-60% identity in the genus-scope trial, MATsc `medium` before |

Proposed: **leave**. A ruling would need a reference record from the genus.

## B. Lean MTL, no flank support: curator to decide

| Genus | genomes (species) | kept | gap | loci with PAP1/OBP1/PIK1 flank |
|---|---|---|---|---|
| Saccharomycopsis | 20 (6) | MTL 22, MATsc 5 | -5.0 | 0 |
| Ambrosiozyma | 13 (11) | MTL 15, MATsc 4 | -3.3 | 0 |
| Citeromyces | 3 (3) | MTL 3 | -3.1 | 0 |
| Saturnispora, Wenyingozyma | 1 each | MTL | -5.3, -11.2 | 0 |

The family is split inside the genus for Saccharomycopsis and Ambrosiozyma (a few loci keep MATsc), which the per-locus rule produced and the genome-level rule only partly removes (a genome with a single cross-called locus is decided by it). No genome has two MTL flank genes, so MTL here is a choice by identity alone, the same weakness as the MTL-in-Serinales finding of #188. Proposed: **hold** and check synteny (SLA2/DIC1 against PAP1/OBP1/PIK1 neighbourhood) for one genome per genus before any scope.

## C. Split or tied genera: needs the evidence

| Genus | genomes (species) | kept | gap | flank loci |
|---|---|---|---|---|
| Petasospora | 11 (7) | MATsc 6, MTL 5 | 0.0 | 2 |
| Wickerhamomyces | 9 (9) | MATsc 8, MTL 1 | +4.3 | 0 |
| Barnettozyma | 8 (5) | MTL 6, MATsc 2 | -6.2 | 6 |
| Gotozyma, Komagataea, Ruyongia, Xingzhongia | 2 each | split | -0.2 to +6.4 | 0-1 |
| Saccharomycodes | 3 (2) | MTL 3 | -2.3 | 3 (one gene each) |

Barnettozyma is the one genus where MTL flank genes show up in most loci (6 of 8, never two): worth a synteny look first. Saccharomycodes is close to Saccharomycetaceae by phylogeny, so MATsc would be the expected family although the identity leans MTL. Proposed: **hold**.

## What would settle each

1. Synteny of one genome per genus B and C (a figure with the flanks of both families).
2. A profile score of each locus's genes against both reference sets (alpha-box and HMG-box profiles), instead of best-hit identity.
3. A curated reference record from a genus (as for MATpk) before it is scoped.

## Questions for the curator

1. Accept **leave** for the genera of A?
2. For B and C: **hold** now, with a synteny check as the next step, or scope any genus now?
3. Is the order of work: Barnettozyma, Saccharomycopsis, Ambrosiozyma, Petasospora?
