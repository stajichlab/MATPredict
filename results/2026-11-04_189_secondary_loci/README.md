# #189: the 80 genomes with a secondary-system MAT locus (2026-11-04)

`decisions.tsv`: one **proposed** decision per genome (all `status=proposed`; nothing is final until the curator signs off). The 80 genomes are those flagged `curate` or `curate+reroute` in `../2026-10-28_idiomorph_crosscall_flags/crosscall_genomes.tsv` (`crosscall_v070.tsv`, a copy of the v0.7.0 table).

## Method
Re-ran the 80 genomes with the current detector (`run_80.sh`: v0.8.1 preset, the v0.8.0 database; the merged `resolution_flags` code) and extracted the called loci with their per-gene identity, contig-edge distance and flags (`evidence_report.py` -> `evidence.tsv`, called loci only). Two tiers (curator, 2026-11-04):
- **Tier 1, suppressed by policy (70 genomes).** v0.8.x reports one family, but this is the cross-family policy and routing suppressing the other family, **not an independent resolution** (the decision value is `kept_family_<F>_by_policy`; corrected after the Fable review). The decision records the kept family (`single_family_MATsc` 50, `MATbz` 9, `MATog` 6, `MTL` 4, `MATss` 1) and what v0.7.0 had. The routing changes (#216, #220, the MATog scope) and the cross-family policy (#227) settled these; review a sample.
- **Tier 2, by hand (10 genomes).** v0.8.x still reports two families at two positions, always complementary idiomorphs except one. Flank genes (SLA2, DIC1, SUI1, PIK1, OBP1, PAP1, TGL1) were searched within 50 kb of each locus (`flank10.py` -> `flank10.tsv`, miniprot against `flank_refs_all.faa`).

## Tier 2 result (corrected after an independent Opus review of four exemplars, 2026-11-04)
| decision | genomes | reading |
|---|---|---|
| `second_cassette_supported` | *Ambrosiozyma philentoma* | MATsc MATa locus (DIC1, a2, a1, SLA2, SUI1) and an MTL alpha locus (alpha2, alpha1, SLA2 fragment) on contigs that share no sequence: a complementary pair at two distant positions |
| `second_cassette_candidate_fragment` | *Peterozyma xylosa* and *toletana*, *Williopsis suaveolens*, *Saccharomycopsis* sp. KB-2023a, *Kuraishia capsulata*, *Ambrosiozyma monospora* | MATsc locus with SLA2 plus a short MTL locus of the complementary idiomorph near a contig end with no flank in range; not an allele in *P. xylosa* (no shared sequence); unconfirmed |
| `second_locus_unresolved` | *Ambrosiozyma platypodis* | MATsc MATa locus with flanks plus an MTL alpha locus 4.6 kb from a contig start, no flank |
| `allelic_pair` | *Wickerhamomyces edaphicus* | one heterozygous MAT locus: both ends of the 3.5 kb MTL A contig (3,507 bp; the A locus itself is 2.4 kb) match the MATsc alpha contig at 97.2% and 97.3% on either side of the idiomorph (SLA2 and DIC1); **corrects the first draft, which called it a second cassette** |
| `duplicate_haplotype` | *Saccharomycopsis fibuligera* | one locus on two collinear haplotype contigs (92-98%) of a 29.4 Mb partly separated diploid or hybrid assembly; **corrects `same_idiomorph_two_families`** |

Lessons recorded here: (1) MTLA1 is MATa1 and MTLA2 the HMG MATa2, so the two families find the same genes under different names, and the MATsc `MATA2` label can land on the alpha2 homeodomain; (2) complementary idiomorphs from two families do not by themselves mean two cassettes: compare the loci at the nucleotide level (shared flanks at 95% or more = allelic pair; collinear along the whole contig = duplicate haplotype; no shared sequence = candidate second cassette). Only the exemplars *W. edaphicus*, *P. xylosa*, *A. philentoma* and *S. fibuligera* were reviewed independently; the other six tier 2 rows have not had the nucleotide comparison and may include further allelic pairs.

## Acceptance check
- No genome has two loci of different families overlapping or at the same start (checked on `evidence.tsv`): 0 overlaps, 0 high-confidence loci with an identical start.
- The *Pichia* reference part of #189 is covered by `MATpk` (#216). The genera *Ambrosiozyma*, *Peterozyma*, *Williopsis*, *Kuraishia*, *Wickerhamomyces* have no family or record of their own.
- All MAT genes in tier 2 are at 22-55% identity (partial_core flag); none of this changes a call.

## Open points for the curator
- The cross-family `cassette` label (branch `cassette-cross-family`) labels *W. edaphicus* too; a shared-flank check (allelic pair) is proposed before it is merged.
- Whether to build candidate records for these genera, as for *Ogataea* (#256); *A. philentoma* is the best candidate (a clean two-position pair); *W. edaphicus* would need both idiomorph alleles as one locus record.
- Review a sample of tier 1.

## Nucleotide comparison of every complementary pair (added after the Fable review; `flank_pairs.py`, `flank_pairs.tsv`, column `flank_check` of `decisions.tsv`)
`shared_flank_check` (PR #264) run on all 22 complementary-idiomorph pairs of loci in the 80 genomes, same-family pairs included:
- **Real cassette systems share sequence, so `shared_flank` does not mean allele.** *Saccharomycopsis schoenii* (MATss, homothallic, mirrored SLA2-MAT cassettes) shares 9.7-31 kb at 95.8-100% between its cassettes in six pairs. The check as built in #264 would call these `allelic_pair`.
- **Tier 1 genomes that carry a complementary pair and were not examined by hand before:** *Ambrosiozyma oregonensis* (MTL, 362 bp at 95.0%, a borderline hit on 2.5 kb contigs), *Saccharomycopsis javanensis* (two pairs, 279 and 478 bp), *Xingzhongia strasburgensis* (1,846 bp at 99.7%), *Ruyongia chambardii* (collinear, 6,131 bp at 97.0%): possible allele pairs or duplicate haplotypes; *Barnettozyma vustinii* and *Ogataea pini* share nothing.
- **Borderline:** *Williopsis suaveolens* (a single 125 bp hit at 99.2%) is `shared_flank` only by the 100 bp threshold; the tier 2 decision stays `second_cassette_candidate_fragment` with the check recorded.
Nothing here is confirmed by read depth; ploidy of the genomes is not known.
