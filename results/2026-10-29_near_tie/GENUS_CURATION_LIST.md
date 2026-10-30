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

## Synteny check (curator's order: Barnettozyma, Saccharomycopsis, Ambrosiozyma, Petasospora; all 156 loci were run)

`flank_synteny.py` cuts each ambiguous locus +/- 40 kb and searches 35 flank proteins from the curated records with miniprot (SLA2, CHA1, APN2, DIC1 from the MATsc, MATyl and MATog records; PAP1, OBP1, PIK1 from the MTL records). A hit counts at 40% identity or more; `flank_synteny.tsv` has the nearest hit of each gene (identity@distance to the locus edge, 0 = inside the locus span). Cluster job 31430330.

| Genus | loci | kept MATsc / MTL | SLA2 or CHA1 within 10 kb | DIC1 within 10 kb | PIK1, OBP1 or PAP1 within 20 kb |
|---|---|---|---|---|---|
| Barnettozyma | 8 | 2 / 6 | 7 | 7 | **0** |
| Saccharomycopsis | 27 | 5 / 22 | 21 | 0 | **0** |
| Ambrosiozyma | 19 | 4 / 15 | 13 | 13 | **0** |
| Petasospora | 11 | 6 / 5 | 11 | 11 | **0** |
| Wickerhamomyces | 9 | 8 / 1 | 9 | 9 | **0** |
| all 27 genera | 156 | 91 / 65 | 129 | 71 | **0** |

Findings:

1. **No MTL flank gene lies near any of the 156 loci, including the 65 where MTL was kept.** The MTL architecture has PAP1, OBP1 and PIK1 beside the locus; none is found within 20 kb at 40% identity or more. (The detector's own flank rule found one flank gene beside 38 MTL calls, never two; a hit below 40% identity is not counted here, and a divergent PIK1 in a distant genus could be missed. This is the weak point of the check.)
2. **SLA2 or CHA1 sits beside the locus in 129 of 156 loci**, and **DIC1 beside 71**, the same side-by-side as the *Komagataella* (MATkp) and *Ogataea* (MATog) cassettes: SLA2 on one side, DIC1 on the other. DIC1 is present in the Pichiales and Phaffomycetales genera (Barnettozyma 7 of 8, Ambrosiozyma 13 of 19, Petasospora 11 of 11, Wickerhamomyces 9 of 9) and absent in Saccharomycopsis and Saccharomycodes (Ascoideales, Saccharomycodales), which have SLA2 only.
3. So the identity preference for MTL in Saccharomycopsis, Ambrosiozyma, Barnettozyma and Citeromyces is not backed by the neighbourhood: these loci sit in an SLA2 (+ DIC1) neighbourhood, a MATsc, MATog or MATkp-like arrangement. The MTL calls are probably the closest reference by sequence only.

### Proposed (for the curator; not implemented)

- **Genera with SLA2 + DIC1 beside the locus** (Barnettozyma, Ambrosiozyma, Petasospora, Wickerhamomyces, Hansenula, Williopsis, Ruyongia, ...): a MATog/MATkp-type architecture. Candidate action: a reference record from one genus (a *Barnettozyma* or *Petasospora* species with a manual locus check) before scoping; until then the MATsc/MATog calls stand.
- **Rule change to consider:** an MTL call outside Serinales with no flank support should not win a cross-call against a family whose own flank genes are present (MATsc: SLA2). Applied to this replay it would keep MATsc in all 65 loci where MTL won. It makes the flank rule decide the family, the asymmetric evidence that the evaluation warned about, but here the asymmetry is complete (0 of 156). This needs the curator's decision.
