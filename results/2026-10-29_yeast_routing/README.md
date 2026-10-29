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
| 0 | 15 | no PAP1, OBP1 or PIK1. In 8 the locus is within about 1.3 kb of a contig end or the contig is under 7 kb, so the check says nothing (GCA_030558145.1, GCA_047716225.1, GCA_003706705.3, GCA_030267685.1, GCA_964187715.1, GCA_030564985.1, GCA_031125185.1, GCA_019232045.1). One more (GCA_030267585.1, Ogataea) has 1.3 kb on one side and 5.8 kb on the other, which is limited room. In 6 there is room (12 to 40 kb or more on at least one side) and no flank was found: GCA_047716205.1, GCA_031125705.1 and GCA_030583325.1 (Saccharomycopsis), GCA_004195225.1 (Petasospora), GCA_003706265.2 (Nakazawaea), GCA_030574195.1 (Groenewaldozyma). |

The 6 with room and no flank may carry a different layout (these genera are outside Serinales, where the
neighbours need not be PAP1, OBP1 and PIK1), or may be false matches from MTLA1 and MTLA2 homology; the
data here does not separate the two. They are the ones to check against a reference genome of the same genus.

**Decision (curator): all 22 stay as calls, flagged `family_unconfirmed`** (a flag separate from confidence;
"provisional" is already a reads-side word). The flag becomes `family_confirmed` when two of PAP1, OBP1 and PIK1 lie
within 15 kb of the MTL genes, which promotes the five three-flank loci (`RULE_DRAFT.md`).

## What this does not decide

- Which of the 15 MTL-only calls without flank genes are real (see above). The routing rule in
  `RULE_DRAFT.md` keeps MTL searchable and flags them `family_unconfirmed`.
- Whether two MATog loci in Ogataea are expected. Needs the Ogataea polymorpha reference layout.
- The 35 genomes still showing two MATsc loci under v1.

The routing change itself is a database change and needs a curator ruling; nothing here edits `db/` or `src/`.

## Replay of the implementation (branch `family-unconfirmed`, `--preset v0.8`)

The same 354 genomes with the rule implemented in the detector (`v3_replay_compare.py`, output `v3_replay_summary.txt`):

| group | n | two or more loci: v0.7.0 / v0.8 | no call: v0.7.0 / v0.8 | confidence changed |
|---|---|---|---|---|
| cross-called | 181 | 181 / 45 | 0 / 0 | 0 |
| MTL only | 22 | 0 / 0 | 0 / 0 | 0 |
| MATsc only (control) | 40 | 3 / 3 | 0 / 0 | 0 |
| no call (control) | 25 | 0 / 0 | 25 / 25 | 0 |
| Ogataea (unchanged until the reference check) | 99 | 31 / 31 | 6 / 6 | 0 |

- No MTL-only call is lost; 5 are `family_confirmed` and 17 `family_unconfirmed`, as predicted from the flank search.
- The no-call control stays empty and no confidence value changes.
- Of the 181 cross-called genomes: 114 end with one MATsc locus, 22 with one MTL locus, 23 keep a MATsc and an MTL locus on
  different positions (none overlap; they are the secondary-locus cases for #189), and 22 have only MATsc loci or only MTL loci in two copies.
- **For the curator**: in the 22 genomes where MTL is the only family left, MTL won the cross-call on score because neither family is
  scoped to the genome. v1 (MATsc always kept) would have kept MATsc there. They carry `family_unconfirmed` and the `cross_call` note names MATsc.

## Cassette indication (curator's decision, 2026-10-29)

Two non-overlapping MATog calls of different idiomorphs in one genome are each indicated as a cassette (`MATa_cassette`,
`MATalpha_cassette`, role `cassette_undetermined`; which one is expressed is not determined). Replay on the 99 Ogataea genomes with
Ogataea added to the MATog scope in a throwaway copy of the database (`cassette_replay.py`, `cassette_replay_summary.txt`):
53 genomes carry both cassettes, 38 one locus (not marked), 7 no call, and the hybrid assembly with four MATalpha loci is not marked.
With Ogataea routed to MATog no MATsc call remains in those genomes. The scope change itself is not in PR #214; it is a routing decision.

## MTL against MATsc where MTL wins (25 genomes)

`mtl_vs_matsc_side_by_side.txt` lists both calls of each genome (same coordinates in every case). MTL wins on the tier
(MTL `high`, MATsc `medium`); neither family has a profile classifier (`idiomorph_classifier` is null for both), so no classifier scores exist to compare.
By the best core-gene identity of each family's call: MATsc is higher in 13 genomes (Ogataea 5, Williopsis 2, Taiozyma 2, Wickerhamomyces, Barnettozyma, Kregervanrija, [Candida]), MTL in 10 (Saccharomycopsis 5, Barnettozyma 2, Petasospora, Gotozyma, Ambrosiozyma), 2 are equal. Identities are measured against different reference sets, so they are not strictly comparable.

## Replay with the decided rules (branch `family-unconfirmed`, `--preset v0.8`; `v5_replay_summary.txt`)

The curator's decisions of 2026-10-29 are in: cross-call ties go to the higher best core identity (within 3 points the family without a
confirmation rule, MATsc, wins and the call is flagged `family_unconfirmed`), and Ogataea (461281) is in the MATog scope.

| group | n | two or more loci: v0.7.0 / now | no call: v0.7.0 / now | confidence changed |
|---|---|---|---|---|
| cross-called | 181 | 181 / 51 | 0 / 0 | 7 (Ogataea, now MATog) |
| MTL only | 22 | 0 / 0 | 0 / **1** | 0 |
| MATsc only (control) | 40 | 3 / 6 | 0 / 0 | 2 (Ogataea) |
| no call (control) | 25 | 0 / 0 | 25 / 25 | 0 |
| Ogataea | 99 | 31 / 54 | 6 / 7 | 13 |

- All 22 genomes whose confidence changed are Ogataea genomes that MATog now calls instead of MATsc; no other genome changes confidence.
- The 25 genomes where MTL previously won on tier: 8 now keep MATsc, 3 keep MATsc flagged as a near tie, 8 keep MTL (flagged), 5 are Ogataea and
  now MATog, and 1 keeps both families as separate loci.
- **One call is lost**: GCA_019232045.1 (Ogataea, Wolfe_Odeg_UCD465). It had an MTL-only alpha locus 73 bp from a contig end; with Ogataea routed to MATog,
  MTL is no longer searched and MATog finds nothing. One genome of 99, but it is a real loss of a call that the Ogataea scope change causes.
- **The flag is common**: in the 181 cross-called genomes 103 of 239 loci are `family_unconfirmed` (MTL 61, MATsc near-tie 42). The near-tie flag
  marks more calls than the MTL-outside-scope flag did (the earlier replay: 51). The curator may want a stricter near-tie margin than 3 points.
- Ogataea with two or more loci goes from 31 to 54: these are the MATa and MATalpha cassette pairs.
