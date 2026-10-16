# Alpha-box and HMG-box profiles: what they show (2026-10-16)

Profiles: `db/Ascomycota/classifiers/MAT_core/` (MAT1-1-1 alpha box, 22 curated sequences, median pairwise identity
22%; MAT1-2-1 HMG box, 26 sequences, 26%). Built by `scripts/build_asco_core_hmms.py` from curated records only,
deterministic, leave-one-genus-out 47/48 (the one failure scored 0 on both). Not referenced by `order.yml`: no call
changes. Scoring: `scripts/score_core_profiles.py`; a locus is called by the profile when one profile beats the other
by at least 25 bits (provisional, the Mucoromycota `min_margin`). Descriptive; the profile is not ground truth.

## Controls (profile call against campaign label)
- 200 random alpha-only loci: 199 MAT1-1, 1 undetermined. 200 HMG-only: 198 MAT1-2, 2 undetermined.
- 200 clear both-gene loci (gap of 100 bits or more): the profile agrees with the label in 188.

## The 204 loci whose label is not the higher-bitscore gene (profile call)
| pattern | n | profile call |
|---|---|---|
| Pleosporales, alpha 31 bits and HMG 55 bits | 87 | MAT1-2 in all 87 |
| exact tie between MAT1-1-1 and MAT1-2-1 | 58 | undetermined 50, MAT1-1 8 |
| small gap (6 to 40 bits) | 59 | MAT1-2 in 52, undetermined 7 |

The label is MAT1-1 in 191 of the 204 and 13 undetermined; the profile agrees with the label in 5.

## Parastagonospora (186 genomes, 182 with a call in each run)
| run | label | n | profile MAT1-1 / MAT1-2 / undetermined |
|---|---|---|---|
| v0.6.0 | MAT1-1 | 178 | 91 / 87 / 0 |
| v0.6.0 | MAT1-2 | 4 | 0 / 4 / 0 |
| with records | MAT1-1 | 91 | 91 / 0 / 0 |
| with records | MAT1-2 | 91 | 0 / 91 / 0 |

The profiles confirm the records run: every MAT1-1 and every MAT1-2 call agrees with the profile, giving 91 and 91.
The 87 v0.6.0 MAT1-1 calls that the profile calls MAT1-2 are the 87 changed by the records run.

## Where the wrong labels come from (not a tie-break)
`assign_idiomorph` already returns `undetermined` on an exact tie. The vote takes the best bitscore over every
informative gene of each idiomorph, not just the core genes. In the 204 loci the MAT1-1 side wins on an accessory
gene: MAT1-1-3 in the 87 Pleosporales (the HMG gene named MAT1-1-3), and MAT1-1-2 in 99 of the other 104 (for
example 142 bits at 37.5% identity against a real MAT1-2-1 at 55.8 bits and 44%). MAT1-1-2 carries an HMG-like
domain, so an HMG-box protein hits it.

`replay_vote.py` replays the vote from the reports (it reproduces the campaign label in 12,534 of 12,628 loci) with
accessory genes removed from the vote:
| vote over | loci whose label changes | profile-scored loci agreeing with the profile |
|---|---|---|
| all genes (campaign) | – | 685 of 974 |
| without MAT1-1-2 | 135 (128 to MAT1-2) | 793 (108 gained, 0 lost) |
| without MAT1-1-2 and MAT1-1-3 | 288 (242 to MAT1-2, 45 to undetermined) | 871 (186 gained, 0 lost) |
| core genes only | 305 | 872 (187 gained, 0 lost) |

## Caveats
- The scored set is the ambiguous list plus controls (974 loci), not the whole campaign, so the 0-lost figure is
  for these loci. The 12,628 labels outside it were not scored.
- Dropping the accessory genes could cost a locus whose only evidence is MAT1-1-2 or MAT1-1-3: the replay leaves
  those unchanged (no core gene to vote on), but the effect on such loci is not measured.
- 45 loci become undetermined (the tie group): both core hits are weak and untranslated, so the honest call is
  that there is no evidence.

## Open
- Make MAT1-1-2 and MAT1-1-3 (and the other accessory genes) non-voting for the idiomorph call, or vote only when no
  core gene is found. Requires a replay on the whole campaign and a curator decision.
- Whether to activate the profiles as a classifier for Ascomycota (`idiomorph_classifier` in `order.yml`).
