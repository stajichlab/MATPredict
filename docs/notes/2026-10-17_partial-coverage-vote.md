# Partial-coverage idiomorph vote (2026-10-17)

Branch `partial-coverage-vote` = `origin/ascomycota-core-hmms` (main + PR #86's alpha-box/HMG-box profiles) + `origin/botrytis-curation`
(record `40559_sas405_MAT_MAT1-2`, merged at the curator's request) + this change. Flag `--partial-core-policy {off,coverage_profile}`,
default `off`, recorded in `run.parameters.partial_core_policy`. Code: `src/MATPredict/detect/partial_core.py` (pure rule) and the glue in
`pipeline.py`. Data and scripts: `results/2026-10-17_partial_coverage_vote/`.

## What the rule does

1. **Reference coverage for polished genes.** `GeneEvidence.coverage` (percent, 0-100, as for raw hits) was null for exonerate/miniprot
   models. It is now `min(100, 100 * model_aa / reference_aa)`, where `model_aa` = sum of the canonical model's exon lengths / 3 and
   `reference_aa` is the length of the reference record the model was built from (`record_protein_lengths`, the table `full_length_models` already
   uses). It is a LENGTH ratio, not the aligned query span (the model does not carry the span), so an internal insertion would read slightly high.
   Written to the report and GFF3 attribute `coverage=` and shown in the HTML as before (they already handled a non-null value). Nothing exonerate
   emits changed. Raw hits keep their own search coverage. Only that field is new when the flag is off.
2. **Trigger (coverage rule).** The vote's winner W has a best non-optional core gene (modelled, not `lost_to`) with coverage < PARTIAL (80), W has
   no core gene >= FULL (90), and the other idiomorph O has a modelled non-optional core gene >= FULL.
3. **Edge guard.** If W's partial gene lies within 500 bp of either end of its contig, the vote is kept and `idiomorph_vote_note` says why
   (`kept_edge_artefact`). Either end is tested, since the model does not say which end is missing; this errs toward keeping the old vote.
4. **Profile confirmation.** The alpha-box (MAT1-1-1) and HMG-box (MAT1-2-1) profiles of `db/Ascomycota/classifiers/MAT_core` score each side's own
   protein (`classifier.own_profile_scores`, best over both tools' models). If O's gene beats W's gene by >= 25 bits the call is O (`flipped`).
   If the profiles favour W, are within 25 bits, or cannot be run (no profile files, no pyhmmer, no protein) the call is `undetermined`, with both
   reasons in `idiomorph_vote_note`. The 25-bit margin is PROVISIONAL, borrowed from the Mucoromycota classifier's `min_margin`; the Ascomycota
   profiles are gene-identity profiles (alpha box vs HMG box), not allele profiles, so no margin has been calibrated for them. The profiles are loaded
   here, not through `order.yml` (the manifest says they are deliberately not referenced there), so the `idiomorph_classifier` path is untouched.
5. **Side effects of a change.** When the rule flips or sets `undetermined`: confidence is capped at `medium` (provisional, see questions),
   `idiomorph_margin` is null (the vote's margin describes the set-aside idiomorph), `idiomorph_candidates` keeps the original vote. The rule
   applies in the strict and in the relaxed pass. It does nothing where a Mucoromycota-style HMM classifier already decides (`verdict`), nor to
   optional genes (MAT1-1-2, MAT1-1-3, ...), which stay `--idiomorph-policy`'s business.

`idiomorph_vote_note` (written only when the rule triggered or the guard fired): `decision`, `reason`, winner/other gene with coverage and bitscore,
thresholds, `profile` {scores, margin_bits, call}, `winner_edge_distance`.

## Choice of PARTIAL / FULL (data)

Estimated coverage (exon aa / reference aa, 19,400 v0.6.0 Ascomycota reports, `scan/`) of all modelled non-optional core genes, 10-point bins from 0:

| gene | 0 | 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100+ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| MAT1-1-1 | 8 | 251 | 353 | 1294 | 1829 | 149 | 279 | 150 | 163 | 806 | 1322 |
| MAT1-2-1 | 22 | 493 | 1610 | 1503 | 689 | 149 | 525 | 118 | 69 | 685 | 1852 |

Most genes cover < 50% of their reference (a domain-only model against a long reference from another clade), so a low coverage alone is
ordinary and carries no information; the rule can only fire when the OTHER idiomorph has a full-length gene in the same locus. Full-length genes sit
at >= 90 (the 80-90 bin is the trough: 163 and 69). Contested loci (both idiomorphs modelled, a vote with two scores): 2,015. Trigger counts
(loci = genomes) for the grid of thresholds, v0.6.0 reports, before the edge guard and profiles:

| PARTIAL below | FULL >= 0.8 | FULL >= 0.9 | FULL >= 0.95 |
|---|---|---|---|
| 0.6 | 15 | 12 | 11 |
| 0.7 | 56 | 53 | 17 |
| 0.8 | 57 | **54** | 18 |
| 0.9 | 57 | 62 | 25 |

PARTIAL 0.8 / FULL 0.9 sits on the plateau (0.7-0.8 changes 1 locus; 0.6 loses the Botrytis remnant at 62%). FULL 0.95 would lose most
true cases because reference lengths differ by a few residues. Not tuned on the truth set (the choice was made from this scan before any run).

## Evaluation (same code, `off` vs `on`)

Run: this branch, current `db/` (includes the SAS405 record), `detect --taxid <assembly-summary taxid>`, MATPREDICT_HTML=0, no held-out records
(none of the genomes is a curated record). SLURM `pcv-*` arrays on `preempt`, 2 CPUs, 16 GB, idempotent. A genome's call is its best
`Ascomycota:MAT` locus (high > medium > low, then modelled genes); `undetermined` is its own outcome; `none` = no locus.
Per-genome tables: `truth_off_on/per_genome.tsv`, `campaign_off_on/per_genome.tsv`, all notes: `vote_notes_on.tsv`.

### (a) Botrytis, 34 genomes (17 stated MAT1-2 + T4 expected MAT1-2; 15 + B05.10 expected MAT1-1)

| | correct | wrong | undetermined | no locus |
|---|---|---|---|---|
| off | 16 | 18 | 0 | 0 |
| on | 33 | 1 | 0 | 0 |

17 genomes change MAT1-1 (high) -> MAT1-2 (medium): 16 of the 17 BioProject genomes plus T4. Profile margins 95-98 bits (HMG 220-223 vs alpha
118-127). The 16 stated-MAT1-1 genomes and B05.10 are unchanged. The remaining wrong one is Rf2 (GCA_015148645.1): its MAT1-2-1 is on another
contig than the 222-bp MAT1-1-1 remnant, so no full-length gene shares the locus and the rule cannot see it (expected, noted in the curation note).

### (b) the 203 BioSample truth genomes (+ T4 and B05.10 = 205 runs)

| | correct | wrong | undetermined | no locus |
|---|---|---|---|---|
| off, all 205 | 179 | 22 | 0 | 4 |
| on, all 205 | 196 | 5 | 0 | 4 |
| off = on, 171 non-Botrytis | 163 | 4 | 0 | 4 |

Gains 17 (all Botrytis), losses 0, no non-Botrytis genome changed label or confidence: the rule did not trigger on any of the 171 (the 4 wrong are
Penicillium roqueforti x2, Podospora tetraspora, a Magnaporthales sp.; the 4 no-locus are Cryphonectria x2, P. roqueforti, Morchella conica; all
unchanged). So it was not exercised on Xanthoria/Penicillium/Magnaporthales truth genomes; this is absence of triggers, not evidence of safety there.

### (c) campaign-wide (v0.6.0 reports, 19,400 Ascomycota genomes)

Scan (`scan/`): 2,015 contested loci; **54** meet the trigger at 80/90 (53 labelled MAT1-1 by v0.6.0, 1 MAT1-2), all with one locus per genome.
Re-run flag off vs on on those 54 genomes with the current code and db (the current db differs from v0.6.0's, so not every scan trigger recurs):

| outcome with `on` | genomes |
|---|---|
| MAT1-1 (high) -> MAT1-2 (medium), flipped | 44 |
| MAT1-2 (high) -> MAT1-1 (medium), flipped | 1 (Monilinia vaccinii-corymbosi) |
| -> undetermined (medium) | 2 (Polytolypa hystricis, Calogaya cf. arnoldii) |
| vote kept, contig edge (30 bp) | 1 (Aspergillus fumigatus GCA_020498755.1, a truncated true MAT1-1-1) |
| no note (trigger no longer present with current code/db) | 6 |

By genus (49 notes over the `on` runs of all 243 genomes): Botrytis 31 flipped to MAT1-2 (the 17 above plus 14 other Botrytis assemblies),
Ophiostoma 5, Sclerotinia 4, and one each of Botryotinia, Leptographium, Sporothrix, Intubia, Stromatinia -> MAT1-2; Monilinia -> MAT1-1; Polytolypa and
Calogaya `undetermined`; Aspergillus fumigatus edge-guarded. Profile margin of the 46 flips: 25-50 bits 1 (Monilinia, 40.5), 50-90 bits 5,
>= 90 bits 40; both `undetermined` have |margin| < 15. A 50-bit margin would turn Monilinia into `undetermined`; 70 bits would also do so for
Stromatinia (66).

**These campaign flips have no truth.** They are not in the BioSample truth set, so 44 of the 45 flips are unverified as correct or wrong; the
profile agrees with the rule in each (that is how it flips), so they are not independent evidence. Cases to look at first: Sclerotinia/Stromatinia/
Monilinia (Sclerotiniaceae contain homothallic fused-locus species where both genes are real; the 188/228-bit Monilinia pair and the 159/242-bit
Sclerotinia pair show a strong alpha and a strong HMG profile score in one locus) and Calogaya/Polytolypa (both profiles strong).

### Combinations with `--idiomorph-policy` (optional-gene policy, branch origin/optional-gene-policy, not edited)

The two flags are independent and compose. Unit test on a scratch merge (`tests/detect/test_partial_core_with_idiomorph_policy.py`, skipped on this
branch until the other lands; both branches change the same call sites, the merge conflicts were all keep-both):

| `--idiomorph-policy` | `--partial-core-policy` | call (partial MAT1-1-1 + full MAT1-2-1 + overlapping MAT1-1-3 cross-hit) |
|---|---|---|
| legacy | off | MAT1-1 (MAT1-1-3 supersedes MAT1-2-1) |
| legacy | coverage_profile | MAT1-1 (MAT1-2-1 is superseded, so no full gene remains for this rule) |
| core_first | off | MAT1-1 |
| core_first | coverage_profile | MAT1-2 |

Real data, 34 Botrytis genomes (same scratch merge): legacy/off 16 correct, legacy/on 33, core_first/off 15, core_first/on 32. `core_first` alone
does not fix the 17 (they stay MAT1-1 but drop high -> medium) and **makes Sl8 (GCA_015148665.1, stated MAT1-1) wrong**: its MAT1-2-1 relic (16% of
reference, 87 bits) is a "core" vote once the optional MAT1-1-5 (389 bits) is ignored, so the locus is called MAT1-2. That is the other branch's
behaviour, found here, not caused by this change; this rule did not fire because that locus holds no full-length MAT1-1-1.

## Tests

`tests/detect/test_partial_core_vote.py` (22): the rule (flip when the profiles agree; `undetermined` when they disagree, are inconclusive or
unavailable; a truncated true MAT1-1-1 at a contig edge does not flip; two full-length genes, two partial genes, unknown coverage, ties unchanged;
default off unchanged; coverage values; report field; relaxed pass) through `run_pipeline` with stubbed polishing and profile scoring, plus the pure
function. `tests/detect` + `tests/db` on the HPCC `test` pixi env (pyhmmer available), excluding `test_draw.py`/`test_synteny.py` (no pygenomeviz):
1082 passed, 4 skipped (the four-combination test skips without `--idiomorph-policy`; it passes on the scratch merge, 947 tests/detect passed there).

## Assumptions and what was not verified

- The coverage is a length ratio against the record the model was built from. A reference that is itself a fragment, or a clade whose reference
  proteins are much longer than the real gene, lowers coverage; the "full" side then fails to trigger (safe direction), the "partial" side may fire
  wrongly only if the other idiomorph is also full length in the locus.
- PARTIAL 80 / FULL 90 / EDGE 500 bp / 25 bits are all provisional. Edge distance uses either contig end.
- The `medium` cap after a change is my choice, not a ruling.
- Truth is the submitters' statement; Rf2 and the other truth disagreements outside Botrytis were not examined. Nothing on the 45 campaign flips was
  checked against literature or reads. The Botrytis calls rest on the curation notes' biology (not re-verified here).
- v0.6.0 scan coverage is an estimate (exon length over the current db's reference length); the campaign re-run used the current db and code.
- Not run: leave-one-record-out held-out versions (the Botrytis SAS405 record is in the db for every run; without it the scores were 222 vs 416, the
  rule's trigger does not depend on it), distant-clade truth other than the 203, Xanthoria (Teloschistales) beyond one Calogaya locus.
- `off` and `on` for the final code differ only by the flag; the `on` arrays were re-run after the relaxed-pass change, the `off` arrays (and the
  `core_first` Botrytis runs) predate it and are unaffected by it (the relaxed-pass hook is inert when the flag is off).
