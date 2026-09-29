# sexM-like paralogs, strain-name labels and Absidia blakesleeana
Status: decided (R4 built; disputed labels set; Absidia known gap); some items open

## Question
1. Can strain-name mating-type labels ("+", "-", "Plus", "Minus") serve as
   known answers, and where do calls disagree?
2. Are the disagreeing Minus calls HMG paralogs, and can a rule remove them?
3. Why is Absidia blakesleeana never called?

## Data and code version
- Reports: LCG held-out (`results/2026-09-28_lcg_holdout/runs2`, code 076afe4),
  Jena, the Mucoromycota scan and Zygo 23 (code a9fb69c / 9a458de).
- Label convention: Plus/Minus (curator, 2026-09-29).

## Method
1. Parse labels from strain names; exclude leaked strains; score calls.
2. Replay stricter flank variants of the MAT-gene gate (V1-V3, V2b).
3. Extract the called "sexM" proteins, compare across opposite-labelled
   strains and genera, check synteny, and replay candidate rules.
4. Trace the A. blakesleeana loci with direct tblastn of the flank genes.

## Results

### Strain-name labels
Source: `results/2026-09-29_strain_labels_and_absidia/NOTE.md`.
- 25 labelled genomes, all in LCG (the one BFD label is a curated-record strain;
  Jena IDs carry none). Clean set n = 21: agree 10, disagree 7, uncalled 4.
- By confidence: high 6 agree / 2 disagree; medium 4 agree / 5 disagree.

### Stricter flank rules do not help
Source: `results/2026-09-29_paralog_partial_check/NOTE.md`.
- V1/V3 (>= 2 flanks, not rnhA alone): 0 of 5 targets withheld, 0 collateral.
- V2b (any partial call below margin 75/100): 3 targets withheld, but 46-56 LCG
  calls withheld and clean labels fall to 8 agree / 2 disagree / 9 uncalled.
- The partial targets pass on the absolute-score route (109-134 bits).

### Two different causes
Source: `results/2026-09-29_sexM_like_paralog/NOTE.md`.
- **Partial sexM+rnhA calls (Circinella, Helicostylum, Thamnostylum) are real
  sexM.** The same protein is 100% identical between a Plus-labelled and a
  Minus-labelled strain, and the "Plus" strains carry no sexP. This points to
  label or strain problems, not detection. (This supersedes the earlier
  "leans detection" verdict in the strain-label note.)
- **Weak second Minus calls are at least three non-MAT HMG genes.** P1 (73 aa)
  is present in both mating types (85% full-length copy in the Plus-only
  M. aligarensis NRRL 3099), sits with tptA/glrA/algA, scores 77-84 bits as
  sexM, and is the only call in 2 LCG and 3 Jena genomes.
- Rule replay: R4 (a P1 class) withholds exactly the P1 calls (17 LCG, all 6
  Jena) with Zygo 23/23; R1 (flank route must include rnhA) removes every
  Umbelopsis call; R2 (core >= 100 aa) loses Zygo loci. R4 was built (see
  classifier-builds report).
- Z. heterogamus NRRL 1489: sexP is present (scaffold_2088, 4.8 kb; annotated
  protein scores sexP 296.9) but its cluster is core-only and never admitted;
  sexM sits on another scaffold. The assembly splits the locus.
- S. megalocarpus SC16: the Minus locus (sexM 100%, score 228.3) is withheld at
  the 0.5 fraction floor. A strong-core floor rescue would add 28 LCG loci (21
  as a second idiomorph): on hold.
- Strain identity: 10 nominally different Mucor species with identical profiles
  are one M. indicus lineage; "R. arrhizus" NRRL 1470 is C. umbellata;
  "T. repens" NRRL 6240 is Circinomucor circinelloides.

### Absidia blakesleeana
Source: `results/2026-09-29_strain_labels_and_absidia/NOTE.md`.
- Correction: 2 of 3 genomes have a strong annotated sexP (NRRL 1300 margin
  107.8; NRRL 1303 margin 146.0); NRRL 1301 does not.
- Every flank gene is present genome-wide (tptA 74.0%, glrA 71.9%, rnhA ~47%,
  algA 42.0%, btbA ~34%), but none lies within 100 kb of sexP. The sexP cluster
  holds only core genes and is never admitted.
- This matches Schulz et al. 2016 (rnhA on another scaffold in heterothallic
  Absidia). Cause: divergent gene order (biology).

### Disputed and misidentified strains (curator ruling 2026-09-29)
`results/2026-09-29_strain_labels_and_absidia/disputed_labels.tsv` (origin of
labels unknown): Ellisomyces RSA_581-, Gilbertella CBS_442.64-, Pirella
RSA_622-, Circinella angarensis RSA_198_Plus, C. umbellata RSA_505_Plus,
Thamnostylum repens RSA_459_Plus, Backusella lamprospora NRRL_6044_Plus.
`misidentified_strains.tsv`: B. ctenidia NRRL 6239, R. arrhizus NRRL 1470,
T. repens NRRL 6240 (plus the M. indicus lineage).

## What changed in detection
41bd471 (R4 P1 paralog class). No flank-rule variant adopted.

## Limits
- Only 25 labelled genomes; many labels are now disputed.
- P1 trained on one sequence; P1b and the Z. exponens gene are not covered.
- Absidia conclusion rests on 2 genomes with a strong sexP.
- The LCG `genus_taxonomy.tsv` mapped "Absidia" to a beetle taxid (255787); fixed
  to 4828 on 2026-09-29 (scoring was unaffected).

## Curator decisions
- Made: Plus/Minus convention; disputed labels excluded; misidentified strains
  excluded from species-level scoring; build R4; hold the floor rescue;
  A. blakesleeana is a known gap (option c); core-only admission (option b) and
  the rnhA-flank hypothesis go to the Lichtheimiaceae exploration.
- Open: how to handle the misidentified strains in BFD metadata.

## Files
`results/2026-09-29_strain_labels_and_absidia/` (lcg_label_scores.tsv,
score_output.txt, disputed_labels.tsv, misidentified_strains.tsv,
absidia/flank_positions.txt); `results/2026-09-29_paralog_partial_check/`
(replay_output.txt); `results/2026-09-29_sexM_like_paralog/` (models.tsv,
cross_type.tsv, replay_rules.tsv, p1_replay.tsv, floor_rescue.tsv).
