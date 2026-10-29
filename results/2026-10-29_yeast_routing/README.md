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
   26 or so MTL reference records are all from Serinales. As written, v1 trades a duplicate for a missing call.
3. **v2 moves Ogataea to MATog**: all 35 Ogataea genomes with a cross call and all 9 Ogataea in the
   MATsc-only control change from MATsc to MATog. The "MATsc-only control" is therefore not a pure
   control for Ogataea. Two of the MTL-only genomes gain a MATog call (Ogataea).
   The count of genomes with two loci goes up in the Ogataea group (31 to 54); Ogataea carries the
   MATa and MATalpha cassettes at separate positions, so two MATog loci may be correct biology
   and need a reference check before they are called duplicates.
4. The no-call control stays empty in both versions, so routing did not introduce calls where
   there were none.

## What this does not decide

- Whether the 22 MTL-only calls are genuine. If they are, MTL has to stay searchable for them and
  the cross-call needs a different rule (prefer the family whose taxonomic scope contains the genome,
  fall back to MTL only when the scoped family has no call).
- Whether two MATog loci in Ogataea are expected. Needs the Ogataea polymorpha reference layout.
- The 35 genomes still showing two MATsc loci under v1.

The routing change itself is a database change and needs a curator ruling; nothing here edits `db/` or `src/`.
