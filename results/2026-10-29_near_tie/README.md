# The near-tie margin: replaced by a genome-level decision (curator decision of 2026-10-29)

`--cross-family-policy on` (PR #214) kept, for two families calling one locus where neither is scoped to the genome, the call with the higher best core-gene identity, except that within `NEAR_TIE_POINTS` = 3 points the family without a `confirm_outside_scope` rule (MATsc) was kept and flagged `family_unconfirmed`. The 3 was a draft value, never validated. The curator asked for it to be evaluated.

## Evidence

`cross_call_pairs.tsv`: the 156 cross-called locus pairs of the 354-genome replay (v5), 140 genomes, all MATsc against MTL; identity gap = kept minus other best core-gene identity. An independent evaluation (a second model, read-only, 23 tool calls, from this table and the full reports) found:

- **No natural break at 3 points.** Per species the gaps run continuously from -3.3 to +1.7. The cutoff splits congeners: *Ambrosiozyma oregonensis* (-2.6) is kept as MATsc, *A. angophorae*, *A. neoplatypodis* and *A. ambrosiae* (-3.2 to -3.3) as MTL.
- **The score is too coarse for a 3-point rule.** The kept call's best alignment has a median 59 aa at 38% coverage, so one residue is about 1.7 points and the standard error of one identity is about 6.5 points. The two identities come from different reference proteins, often different genes (41 of 105 kept MATsc calls take their best identity from MATA2 or MATALPHA2, which the registry marks `idiomorph_informative: false`).
- **The MATsc tie-break is biased.** `sla2` is searched only by MATsc and PAP1/OBP1/PIK1 only by MTL, so gene counts and flank evidence are not comparable; no MTL call reached two flank genes. In 4 near ties the two families call different idiomorphs, so the tie-break silently sets MATa against alpha.
- **Per-locus decisions give genomes with two families** (4 genomes in v5), which cannot be right for one genome.
- **The flag measured "MTL was involved", not uncertainty**: every MTL call outside Serinales is flagged by the flank rule anyway, and every margin option flags 65-108 calls.

The evaluation's counts for 156 pairs: margin 0: 5 near ties; 1: 19; 2: 34; 3: 42; 5: 72; genome-level mean with no margin: 0 mixed genomes and 1 lower-identity call kept, against 4 and 12 at 3 points.

## Decision (option A)

1. No margin and no preference for either family.
2. For every group of overlapping cross-family calls with no family scoped to the genome, the family is decided once for the genome: the higher mean, over those groups, of each family's best core-gene identity (equal means: the smaller family name).
3. Every such kept call carries `family_ambiguous` (identities of both families, their gap, the genome means and locus counts): GFF3 attribute `family_ambiguous=<family:identity|...>` and a YAML block. `family_unconfirmed` is now only the PAP1/OBP1/PIK1 flank rule for MTL outside Serinales.
4. The set-aside call's gene evidence (gene, role, identity, coverage, status) is kept in its `cross_call` note, so the two families can be compared gene for gene after a run.

## Replay with the new rule (job 31429403; `compare_v5_v8.py`, `compare_v5_v8.txt`)

354 genomes, `--preset v0.8`, the database of `main` (no MATpk):

| | v5 (3-point margin) | v8 (genome-level) |
|---|---|---|
| genomes whose cross-called loci keep two families | 4 | 0 |
| `cross_call_near_tie` flags | 42 | 0 |
| `family_ambiguous` loci | - | 156 in 140 genomes (MATsc kept 91, MTL 65) |
| kept call has the lower identity | 12 | 1 (gap -1.5) |
| `mtl_outside_scope` unconfirmed / confirmed | 75 / 5 | 89 / 5 |
| genomes whose calls change | | 14 |

The 14 genomes change family (and so idiomorph names, MATa against A and so on); for most the per-locus rule had kept MATsc where the genus leans to MTL. The extra 14 `mtl_outside_scope` flags are those MTL calls, which have no PAP1/OBP1/PIK1 support.

## Open

- There is still no ground truth: the family is a choice between two reference systems by best-hit identity. The evaluation suggests per-genus curator rulings (Ambrosiozyma, Saccharomycopsis, Pichia, the Wickerhamomycetaceae, as for MATkp, MATog and MATpk), HMG/alpha-box profile scoring against both reference sets, and synteny of a sample per genus.
- With the genome-level rule a genome with only one cross-called locus is decided by that locus alone (as before).
