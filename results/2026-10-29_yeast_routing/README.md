# Yeast routing experiment (#188 / #189)

Two candidate routing changes, run on the 354 yeast genomes that the v0.7.0 run called
(181 cross-called, 22 MTL-only, 40 MATsc-only control, 25 none control, 99 Ogataea; the groups
overlap by 13 genomes, so the group rows sum to 367). Same code as v0.7.0 (frozen worktree
`run-b317ece`), only `db/Ascomycota/order.yml` and the family registry differ.

- **v1**: MTL gets `fallback_searchable: false`, so a genome outside Serinales is not searched with MTL.
- **v2**: v1, plus Ogataea (taxid 461281) added to the MATog taxonomic scope.

`per_genome_calls.tsv` has the family calls per genome for v0.7.0, v1 and v2. `cmp_routing.py`
makes it from the three result folders on HPCC. Both jobs finished (354 of 354 genomes each).

## Result

| group | n | v1 same as v0.7.0 | v2 same as v0.7.0 | two or more loci: v0.7.0 / v1 / v2 | no call: v0.7.0 / v1 / v2 |
|---|---|---|---|---|---|
| cross-called | 181 | 0 | 0 | 181 / 35 / 42 | 0 / 1 / 1 |
| MTL only | 22 | 0 | 0 | 0 / 0 / 0 | 0 / **22** / 20 |
| MATsc only (control) | 40 | 40 | 31 | 3 / 3 / 6 | 0 / 0 / 0 |
| no call (control) | 25 | 25 | 25 | 0 / 0 / 0 | 25 / 25 / 25 |
| Ogataea | 99 | 96 | 43 | 31 / 31 / 54 | 6 / 9 / 7 |

1. **v1 removes the cross-call duplicates** (181 genomes with a MATsc and an MTL locus at the same
   place become 146 with one MATsc locus). The 35 that still have two loci have two MATsc loci
   (Ogataea 12, Pichia 8, Ambrosiozyma 6, and others); those are the genomes for the curation list in #189.
2. **v1 loses every MTL-only call** (22 of 22: Saccharomycopsis 6, Barnettozyma 2, Ogataea 3,
   and others). Those genomes had no MATsc locus, so MTL was the only family that found anything.
   Whether those are real mating-type loci or MTL matching a distant homolog is not known; the
   reference records of MTL were not checked here. As written, v1 trades a duplicate for a missing call.
3. **v2 moves Ogataea to MATog**: all 35 Ogataea genomes with a cross call and all 9 Ogataea in the
   MATsc-only control change from MATsc to MATog. The "MATsc-only control" is therefore not a pure
   control for Ogataea. Two of the MTL-only genomes gain a MATog call (Ogataea).
   The count of genomes with two loci goes up in the Ogataea group (31 to 54); two MATog loci in one Ogataea genome
   may be real (not checked here), so they need a reference check before they are called duplicates.
4. The no-call control stays empty in both versions, so routing did not introduce calls where
   there were none.

## The 22 MTL-only calls: flanking-gene check

`flank_search.py` takes 60 kb either side of each v0.7.0 MTL locus and searches it with miniprot for the
Candida albicans MTL neighbours PAP1, OBP1 and PIK1 (`flank_refs.faa`, from the database records
5476_sc5314_MTL_A, 5476_sc5314_MTL_alpha and 36911_cl143_MTL_A). Results are in `flank_search.tsv`
(contig length and distance to each contig end are included). The search reproduces every flank gene the
detector had already found and adds none.

| flank genes found | genomes | what it means |
|---|---|---|
| 3 of 3 (PAP1, OBP1, PIK1) | 5 | the C. albicans arrangement (PAP1, OBP1, PIK1, then the MTL genes, within about 8 kb) in five independent genomes: GCA_030705185.1 and GCA_030571215.1 ([Candida]), GCA_029290875.1, GCA_029300705.1 and GCA_029300795.1 (Saccharomycetales). Genuine MTL loci. |
| 1 of 3 (PIK1) | 2 | GCA_029296565.1 and GCA_029293515.1 (same contig length and coordinates, so near-identical assemblies). PIK1 then MTLA2, MTLA1, 232 bp from the contig end, so PAP1 and OBP1 would be off the contig. Consistent with a real locus; cannot be confirmed. |
| 0 | 15 | no PAP1, OBP1 or PIK1. In 7 of these the locus is within about 600 bp of a contig end or the contig is under 7 kb, so the check says nothing. In 8 there is room (12 to 40 kb on at least one side) and no flank was found: GCA_047716205.1, GCA_031125705.1 (Saccharomycopsis), GCA_004195225.1 (Petasospora), GCA_003706265.2 (Nakazawaea), GCA_030574195.1 (Groenewaldozyma), GCA_030583325.1 (Saccharomycopsis), GCA_003706705.3 (Barnettozyma), GCA_030564985.1 (Saccharomycopsis). |

The 8 with room and no flank may carry a different layout (these genera are outside Serinales, where the
neighbours need not be PAP1, OBP1 and PIK1), or may be false matches from MTLA1 and MTLA2 homology; the
data here does not separate the two. They are the ones to check against a reference genome of the same genus.

**Decision (curator, in this conversation): all 22 stay as calls but are marked provisional**, with the five
three-flank loci eligible for promotion by the flank rule in `RULE_DRAFT.md`.

## What this does not decide

- Which of the 15 MTL-only calls without flank genes are real (see above). The routing rule in
  `RULE_DRAFT.md` keeps MTL searchable and marks them provisional.
- Whether two MATog loci in Ogataea are expected. Needs the Ogataea polymorpha reference layout.
- The 35 genomes still showing two MATsc loci under v1.

The routing change itself is a database change and needs a curator ruling; nothing here edits `db/` or `src/`.
