# #189: the 80 genomes with a secondary-system MAT locus (2026-11-04)

`decisions.tsv`: one **proposed** decision per genome (all `status=proposed`; nothing is final until the curator signs off). The 80 genomes are those flagged `curate` or `curate+reroute` in `../2026-10-28_idiomorph_crosscall_flags/crosscall_genomes.tsv` (`crosscall_v070.tsv`, a copy of the v0.7.0 table).

## Method
Re-ran the 80 genomes with the current detector (`run_80.sh`: v0.8.1 preset, the v0.8.0 database; the merged `resolution_flags` code) and extracted the called loci with their per-gene identity, contig-edge distance and flags (`evidence_report.py` -> `evidence.tsv`, called loci only). Two tiers (curator, 2026-11-04):
- **Tier 1, policy (70 genomes).** v0.8.x reports one family. The decision records the kept family (`single_family_MATsc` 50, `MATbz` 9, `MATog` 6, `MTL` 4, `MATss` 1) and what v0.7.0 had. The routing changes (#216, #220, the MATog scope) and the cross-family policy (#227) settled these; review a sample.
- **Tier 2, by hand (10 genomes).** v0.8.x still reports two families at two positions, always complementary idiomorphs except one. Flank genes (SLA2, DIC1, SUI1, PIK1, OBP1, PAP1, TGL1) were searched within 50 kb of each locus (`flank10.py` -> `flank10.tsv`, miniprot against `flank_refs_all.faa`).

## Tier 2 result
| decision | genomes | reading |
|---|---|---|
| `second_cassette_supported` | *Wickerhamomyces edaphicus* | MATsc alpha locus (SLA2, SUI1, DIC1) and an MTL A locus with SLA2 584 bp upstream: a SLA2-MATa cassette end, the same arrangement as the *Ogataea* two-cassette genomes |
| `second_cassette_candidate_fragment` | *Peterozyma xylosa* and *toletana*, *Williopsis suaveolens*, *Saccharomycopsis* sp. KB-2023a, *Kuraishia capsulata*, *Ambrosiozyma monospora* | MATsc locus with SLA2 plus a short (1.5-2 kb) MTL locus of the complementary idiomorph, mostly within 1.1 kb of a contig end and with no flank in range: probably a truncated second cassette, unconfirmed |
| `second_locus_unresolved` | *Ambrosiozyma philentoma*, *A. platypodis* | MATsc MATa locus with SLA2, SUI1, DIC1; an MTL alpha locus (genes at 41-49% identity) with no SLA2 or DIC1 flank |
| `same_idiomorph_two_families` | *Saccharomycopsis fibuligera* | both loci are alpha: not complementary, so a second copy or a cross-call, not a cassette pair |

## Acceptance check
- No genome has two loci of different families overlapping or at the same start (checked on `evidence.tsv`): 0 overlaps, 0 high-confidence loci with an identical start.
- The *Pichia* reference part of #189 is covered by `MATpk` (#216). The genera *Ambrosiozyma*, *Peterozyma*, *Williopsis*, *Kuraishia*, *Wickerhamomyces* have no family or record of their own.
- All MAT genes in tier 2 are at 22-55% identity (partial_core flag); none of this changes a call.

## Open points for the curator
- Whether the tier 2 `second_cassette_*` decisions should become a label in the output (a second-cassette note) or stay as curation notes.
- Whether to build candidate records for these genera, as for *Ogataea* (#256); *W. edaphicus* is the best candidate.
- Review a sample of tier 1.
