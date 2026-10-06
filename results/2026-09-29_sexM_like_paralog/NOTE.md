# sexM-like calls: true sexM, paralog, or strain problem?
Status: open (rule design awaits curator decision)

## Question
1. Are the partial sexM+rnhA calls with identical scores across opposite-labelled
   strains, and the weak (~30-bit margin) second Minus calls in the two-idiomorph
   genomes, one conserved non-MAT HMG gene misread as sexM?
2. If so, which rule would reject it, and what would the rule cost?
3. Why were the Zygorhynchus heterogamus sexP and the Syzygites megalocarpus Minus
   locus missed?
4. Are the 10 nominally different Mucor species with identical classifier profiles
   one species?

## Data and code version
- Reports: LCG results/2026-09-28_lcg_holdout/runs2 (code 076afe4); Jena
  results/2026-09-28_mucor_jena_holdout/runs_scaffolds (076afe4); Mucoromycota scan
  results/2026-09-28_next_fixes/Mucoromycota_a9fb69c; Zygo 23
  results/2026-09-28_next_fixes/zygo23_a9fb69c/{scaffold,contig}.
- Genomes: ZyGoLife LCG `genomes/<Org>.sorted.fasta`; Jena `annotate_results/*.scaffolds.fa`;
  BFD `input_clean_genomes/<asm>.fa.gz`.
- Classifier HMMs: db/Mucoromycota/classifiers/MAT/{sexM,sexP}.hmm at the remote
  polish-scope-cuts head (copied as clf_*.hmm).
- Read-only. No code, database or worktree was changed.

## Method
1. extract_models.py: translate every polished gene model in every LCG call from the
   genome and the report's exon coordinates (2,665 models, 536 genomes; models.tsv).
2. cross_type.py: tblastn each called core protein (sexM for Minus, sexP for Plus)
   against up to 4 confident opposite-type genomes of the same genus (only that
   idiomorph called, high, mat_locus). A true idiomorph gene has no near-identical
   copy there (1,860 comparisons, 92 genomes; cross_type.tsv).
3. Pairwise identity (local, BLOSUM62) among target proteins; tblastn of LCG sexP
   proteins against the target genomes.
4. replay_rules.py and p1_replay.py: replay candidate rules on existing reports.
5. Markers (sexP, rnhA, glrA, tptA protein identity) for strain identity.
Tree placement in the final ML tree was NOT done; the cross-type test answers the
"present in both mating types" question directly.

## Results

### 1a. Partial sexM+rnhA calls (Circinella, Helicostylum, Thamnostylum): true sexM
- The called "sexM" models are 120 aa (one exon, HMG box), 48-57% to curated sexM.
  They form one clade: C. umbellata = C. angarensis at 90%, Helicostylum 63-66%,
  Thamnostylum 57-61%; only 38-45% to the full-locus sexM of Helicostylum pulchrum
  and Pilaira.
- They are Minus-specific: 1 of 34 partial Minus calls has a >=95% copy in a
  confident Plus genome of its genus (cross_type.tsv). In Circinella the confident
  Plus genomes (C. rigida, C. simplex) carry no copy above 40%.
- The Plus-labelled carriers carry no sexP: best tblastn hits of LCG sexP proteins
  are 34-40% HMG background in C. umbellata RSA_505_Plus, C. angarensis
  RSA_198_Plus, T. repens RSA_459_Plus, Backusella NRRL_6044_Plus and Pilaira
  RSA_1997_Plus.
- The same protein is 100% identical in C. umbellata NRRL 1366 and RSA_505_Plus,
  and 99.2% in C. angarensis RSA_618- and RSA_198_Plus.
- Verdict: not a both-type paralog. These are true sexM alleles. The disagreement
  with the "Plus" labels points to label or strain problems (see 4), not detection.

### 1b. Weak second Minus calls: at least three different non-MAT HMG genes
- They are not one family. P1 (73 aa; Mucor hiemalis/indicus group), P1b (76-98 aa;
  M. globosus / M. racemosus A-19185, 93-95% to each other) and a Zygorhynchus
  exponens gene (119 aa) are 32-43% to each other and 42-46% to a real M. indicus
  sexM (paralog_candidates.faa).
- P1 is present in both mating types: an 85% full-length copy in the Plus-only
  M. aligarensis NRRL 3099 (cross_type.tsv). It always sits with tptA/glrA/algA and
  no rnhA. Its absolute score is 77-84 bits, below the gate threshold; it passes on
  flank support.
- P1 is also the ONLY call in 2 LCG genomes (M. indicus NRRL 13468, 13081) and in
  3 Jena genomes (CBS221_71, CBS223_63, CBS763_74; 92-97% to P1; CBS763_74 also has
  a 219-307-bit sexP-like annotated protein elsewhere). Those Minus labels come from
  a paralog.
- P1b and the Z. exponens gene: not shown in opposite-type genomes (best 34-40%).
  Unresolved.

### 2. Rule replay (replay_rules.tsv, p1_replay.tsv)
Flank route = classifier score below 99.9 bits, or fragment-typed.

| Rule | LCG drop (genomes losing all calls) | Jena | Mucoro293 | Zygo 23 scaffold / contig |
|---|---|---|---|---|
| R1: flank route must include rnhA | 23 (5) | 6 (3) | 14 (12) | 23/23 / 23/23 kept |
| R2a: core model >= 100 aa | 67 (38) | 10 (7) | 22 (18) | loses 7 / 2 |
| R3: R1 or R2a | 20 (4) | 6 (3) | 14 (12) | kept / kept |
| R4: third class "P1 paralog" HMM | 17 (2) | 6 (3) | n/a (training genome) | kept / kept |

- R4's HMM was built from ONE non-held-out BFD P1 copy (GCA_000697295.1, Mucor
  indicus B7402). It flags exactly the P1 calls: 17 of 570 LCG core proteins and all
  6 Jena weak calls, with P1 scores 145-159 against sexM 77-84. No other call is
  flagged; no Zygo, labelled strain or literature positive is affected.
- R1 also drops Absidia inflata (rnhA sits on another scaffold in Absidia, Schulz
  2016), Mucor (= Umbelopsis) ramannianus and every Umbelopsis call in the
  Mucoromycota scan (rnhA is absent from the Umbelopsis locus). R1 needs lineage
  exemptions; R4 does not.
- R2 is unsafe (loses Zygo loci).

### 3. The two misses
- Z. heterogamus NRRL 1489: sexP IS present, on scaffold_2088:2246-3163 (4,830 bp),
  61% over 307 aa to Z. moelleri sexP; the annotated FBX90_004320-T1 scores sexP
  296.9 vs sexM 72.5. Its cluster was core-only (sexP/sexM hits of one gene) and was
  not admitted. The sexM call sits on a different 5,955-bp scaffold_1293. Cause:
  the assembly splits a locus Schulz 2016 describes as one (5.3 kb), and core-only
  clusters are never admitted. Same failure class as Absidia blakesleeana.
- S. megalocarpus SC16: the Minus locus is on scaffold_38:1-9448 (sexM 100%, score
  228.3, margin 160.4), sexM+rnhA+btbA, withheld at the 0.5 fraction floor (tptA,
  algA, glrA not at the locus; Syzygites uses glrA in place of tptA, Idnurm 2011).
  The split-locus rule does not run because the family already has a call (Plus).
  Its sexM is 100% identical to the Syzygites sp. MES 3091 Minus call, so both
  Syzygites genomes carry both idiomorphs.
- A rescue of below-floor loci with a model-typed core >= 99.9 bits and margin >= 25
  would recover S. megalocarpus and Zygorhynchus sp. NRRL 3102, but adds 28 LCG loci
  (21 as a new idiomorph in an already-called genome), 9 Jena, 9 Mucoro293, 0 Zygo
  (floor_rescue.tsv). That would roughly double the two-idiomorph set; not
  recommended without the two_idiomorphs checks.

### 4. Strain identity
- The 10 identical-profile genomes (M. hiemalis NRRL 3140 and A-26125, M. indicus
  NRRL 13132 and 555, M. racemosus NRRL 1427, M. rouxianus NRRL 1430, M. rouxii
  NRRL 1894, Mucor sp. NRRL A-25783 and A-25793, Backusella ctenidia NRRL 6239) are
  100% identical at sexP, rnhA, glrA and tptA. tptA is 100% to M. indicus
  NRRL 13468 and 64-72% to reference M. hiemalis NRRL 1419, M. racemosus NRRL 1504
  and M. circinelloides. They look like one M. indicus lineage (possibly one clone),
  so the other names are likely misidentifications or synonyms. B. ctenidia
  NRRL 6239 is the clearest case.
- Rhizopus arrhizus NRRL 1470 carries C. umbellata's sexM (100%) and rnhA (100%),
  and is 57-58% to other R. arrhizus rnhA: the genome is Circinella umbellata, not
  R. arrhizus.
- Thamnostylum repens NRRL 6240 is 100% identical to Circinomucor circinelloides
  NRRL 22899 at sexP, rnhA, tptA and glrA: the genome is C. circinelloides.

## What changed in detection
Nothing. Candidate rules only.

## Limits
- R4 was trained on one sequence; P1b and the Z. exponens gene are not covered.
- Cross-type test depends on confident opposite-type genomes existing in the genus;
  8 of 42 partial Minus calls had none.
- Models of 73-120 aa are HMG-box fragments; protein identity on them is noisy.
- No tree placement was done. No read or ITS data were used for strain identity.

## Curator decisions
Open: adopt R4 (a curated "sexM-like paralog" class) and whether to extend it to
P1b; how to treat the misidentified genomes (R. arrhizus NRRL 1470, T. repens
NRRL 6240, the M. indicus lineage) in scoring and BFD metadata; whether to pursue
core-only admission (Z. heterogamus, Absidia) and a floor rescue (Syzygites) with
the two_idiomorphs checks.

## Files
models.tsv, models.faa, cross_type.tsv, target_sexM.faa, paralog_candidates.faa,
P1_train_bfd.faa, jena_weak.faa, replay_rules.tsv, p1_replay.tsv,
floor_rescue.tsv; scripts extract_models.py, cross_type.py, replay_rules.py,
p1_replay.py.
