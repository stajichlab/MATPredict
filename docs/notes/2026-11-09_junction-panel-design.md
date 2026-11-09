# Junction panels and "cannot be typed" species: design for #195 and #196 (2026-11-09)

Status: **proposal for the curator to rule on.** No analysis was run for this note. The numbers are from the committed wave tables on main
(`datasets/reads/sra_survey/pilot/wave*/gates_final.tsv`, cap 700, as of 2026-11-08).

## Why presence typing fails here

`reads-type` asks whether a run's reads contain the k-mers unique to one idiomorph. That gives the mating type only if each nucleus carries one idiomorph sequence. Four kinds of genome break that assumption:

- **Silent cassettes.** Every haploid has the sequence of both idiomorphs: the expressed locus plus silent copies.
  - Saccharomycetaceae: MAT plus HML and HMR (*Saccharomyces*, *Kluyveromyces*, *Lachancea*, *Zygosaccharomyces*).
  - *Komagataella* and *Ogataea*: two loci with switching.
  - *S. pombe*: mat1 plus mat2-P and mat3-M.
- **Homothallic, one locus.** Both MAT genes sit in one locus (*Sclerotinia sclerotiorum*).
- **Heterozygous diploids.** MTL a and alpha sit on two homologues (*Candida albicans*, *Pichia kudriavzevii*).
- **Switching cultures** (HO+ yeasts, *S. pombe* h90). A culture is a mix of cells of both types, so "the" active MAT of an isolate is not defined.

## (1) What a junction panel types, and from which reads

The idiomorph-specific region (Ya/Yalpha in *Saccharomyces*, P/M in *S. pombe*) is shared between the expressed locus and the silent copies. So are its immediate borders: X and Z1 in *Saccharomyces*, and the H1/H2 boxes in *S. pombe*. Even a k-mer spanning the Y|Z1 border also exists at a cassette.

Only sequence outside the homology boxes is unique to the expressed locus. In *S. cerevisiae*, the Z side of MAT reaches unique chromosome-III sequence about 330 bp past Y (Z1 plus Z2). Two designs follow, cheapest first.

**A. Copy-number ratio (no new read structure; recommended first).** Count depth on Ya-unique and Yalpha-unique k-mers, which the panels already report as `depth_1` and `depth_2`. With the standard cassettes (HMLalpha, HMRa), each state has its own ratio:

| State | Ya copies | Yalpha copies | Ya:Yalpha |
|---|---|---|---|
| haploid MATa | 2 | 1 | 2:1 |
| haploid MATalpha | 1 | 2 | 1:2 |
| diploid a/alpha | 3 | 3 | 1:1 |

*S. pombe* is the same with P and M. An h−S strain, which lacks mat2/mat3, shows M only.

*C. albicans* has no cassettes, so the ratio is a straight zygosity test: 1:1 is a/alpha, and one side near zero is a homozygous strain after loss of heterozygosity (LOH).

The ratio test uses single reads and works on single-end runs. Where a panel already types, it needs no new download: `depth_1` and `depth_2` are already in the results tables. It breaks on non-standard cassettes (such as HMLa), aneuploidy, mixed libraries, and switching cultures. A switching culture reads as a mixture, so it lands near 1:1 or between the states.

**B. Linked-pair junction panel (only if A fails on the controls).** Count read pairs with one mate on Ya- or Yalpha-unique k-mers and the other mate on MAT-locus-unique outer flank, the Z side beyond Z2. That types the expressed locus directly, whatever the cassettes are. The costs:

- It needs paired runs with inserts of about 400 bp or more.
- It needs new pair-aware counting in `reads-type`.
- The reads are not kept after typing, so every run has to be downloaded again.

Either design goes behind an opt-in policy flag (default `off`), per the detector rules in `CLAUDE.md`.

## (2) What counts as "cannot be typed" (#195), and how it is reported

The gate today gives PROVISIONAL_HOMOTHALLIC to any species with at least 50% `both` calls. That label wraps four different causes. The proposal is one gate reason per cause, each with its own page wording. The curator approves the text.

| Reason (gate word) | Species now | Site wording (draft) | Frequencies shown? |
|---|---|---|---|
| `homothallic_locus` | *S. sclerotiorum* (HOLD:panel_L1_error; `panel_homothallic.py` from the pilot is not wired in) | "Homothallic: both MAT genes in one locus in every genome." | No |
| `silent_cassettes` | released *S. paradoxus* (0 of 299 typed), *S. eubayanus* (2/301), *S. pastorianus* (11/91), *S. mikatae* (0/8), *O. haglerorum* (0/23); **also now labelled PROVISIONAL_HOMOTHALLIC:** *K. marxianus* (40 of 49 both), *K. phaffii* (50/50), *S. kudriavzevii* (20/22); *S. pombe* (HOLD: no unique k-mers) | "Silent mating-type cassettes: both idiomorphs are present in every cell, so reads cannot show the active mating type by presence." With design A or B in place: the active MAT or "switching/diploid" per run. | No (until A or B passes) |
| `zygosity` (exists) | *C. albicans* UNVALIDATED (609 of 695 typed counting a/alpha; 4 of 7 assembly comparisons agree), *P. kudriavzevii* PROVISIONAL (470/475, 3/3) | "(MTL zygosity)", as now; add the share of homozygous strains once the ratio test runs | Zygosity shares |
| `cause_unchecked` | *Ceratocystis fimbriata* (HOLD:panel_L1_error; unidirectional switching is a hypothesis) | "Withheld: cause not yet checked." | No |
| technical (as now) | the 16 `panel_build_failed` species of #196, `panel_L1_error`, low depth | Withheld, with the technical reason | No |

Three species need a ruling on wording alone, with no new analysis. *K. marxianus*, *K. phaffii* and *S. kudriavzevii* are cassette yeasts, but they are published today as "provisional (homothallic)". Their `both` calls come from the cassettes, which is a structural effect, not a biological finding.

**Puzzle to check (no cause known):** *S. uvarum* types 263 of 294 runs with single calls at flank depth 1.7. That is unexpected in a cassette species, and could mean the panel picks up cassette-divergent k-mers. It should stay UNVALIDATED until the cause is known.

## (3) The yeast panels (#196): species, truth sets, checks

**Species and the first step for each:**

| Group | Species | First step |
|---|---|---|
| *Saccharomyces sensu stricto* | *S. paradoxus*, *S. eubayanus*, *S. pastorianus*, *S. mikatae*, *S. uvarum*, *S. kudriavzevii*, *S. arboricola* | Panels give `none` at flank depth 1.1–1.9, not `both`, so the panel itself is suspect (few loci, divergent flanks). Rebuild Ya/Yalpha panels from curated tier A/A1 DB records (#196 point 2), then apply ratio test A. |
| *Zygosaccharomyces* (8) | *bailii*, *favi*, *mellis*, *parabailii*, *pseudobailii*, *rouxii*, *sapae*, *siamensis* | `panel_build_failed`. Several are allodiploid hybrids with two MAT sets. A cause TSV first (#196 point 1). |
| Other failed builds | *Lachancea fermentati*, *L. thermotolerans*, *L. kluyveri*, *K. lactis*, *D. hansenii*, *C. dubliniensis*, *C. viswanathii*, *L. orthopsilosis*, *S. japonicus* | Cause TSV first |
| Diploids | *C. albicans*, *P. kudriavzevii* | Ratio test A on the existing results (`depth_1`, `depth_2`) |
| *S. pombe* | | Panel from the P and M cassettes themselves, not the 36 kb three-locus representative. Then ratio test A. |

*S. cerevisiae* stays skipped on the site by design. It is still the best positive control.

**Truth sets** (the curator picks; see the open questions):

- **Simulated mixes** (the acceptance already written in #195/#196). Reads simulated from references of known state are mixed at 0, 10, 25, 50, 75, 90 and 100%, in 3 replicates:
  - pass if the minor idiomorph is detected at a share of 10% or more in 3 of 3 replicates;
  - pass if the ratio test classifies haploid a, haploid alpha and diploid a/alpha correctly at 10x or more.
- **Public lab strains of known state:**
  - *S. cerevisiae*: BY4741 (MATa), BY4742 (MATalpha), BY4743 (a/alpha), S288C (MATalpha).
  - *S. pombe*: 972 h−, 975 h+, 968 h90. The h90 strain is the expected "switching" result.
- ***C. albicans*:** published MTL zygosity calls for sequenced isolates, in place of assemblies. Assemblies collapse the two haplotypes, which is why only 4 of 7 agree.
- **Positive control for `homothallic_locus`:** *S. sclerotiorum*, with at least 95% of runs `both` and 0 single calls.

**How results are checked:** first the controls above, then agreement with any assembly or literature call, reported with disagreements listed by run. This is the same rule as for the other reads species.

## (4) What is done, what is not, and what each step costs

**Done** (PRs #282, #283, #298, #310):

- **Zygosity gate:** `wave_driver.ZYGOSITY` for *C. albicans* and *P. kudriavzevii*, with words PASS_ZYGOSITY and PROVISIONAL_ZYGOSITY.
- **Cassette species released and measured:** `wave_driver.RELEASED` (the 5 species in the table above) as `UNVALIDATED:released_cassettes`, run to up to 300 runs each. The result is essentially no calls, which confirms that presence typing does not work for this group.
- **Yeasts with a/alpha panels from waves 3–5:**
  - *Brettanomyces bruxellensis* and *Diutina catenulata*: PROVISIONAL.
  - *Hansenula anomala* and *Metschnikowia pulcherrima*: PROVISIONAL_HOMOTHALLIC (to re-check under the cassette rule).
  - *Williopsis jadinii*: UNVALIDATED.

**Not done**, with estimated costs. Compute is on the HPCC; "dev" means code and tests.

| Step | Issue | Dev | Compute | Needs |
|---|---|---|---|---|
| 1. Gate reasons `homothallic_locus` / `silent_cassettes` / `cause_unchecked`, page wording, relabel the 3 species | #195 points 1–2 | ½ day | none (regate from tables) | wording ruling |
| 2. Cause TSV for the 16 failed builds plus the *Saccharomyces* holds | #196 point 1 | ½ day | minutes (read build logs, BLAST a few loci) | — |
| 3. Ratio test A on *C. albicans* and *P. kudriavzevii* | #195 point 3, #196 point 3 | 1 day (opt-in policy, tests) | none: uses `depth_1`/`depth_2` already in the tables | MTL truth set |
| 4. Rebuild *Saccharomyces* and *S. pombe* panels from curated records | #196 point 2 | 1 day | minutes per panel | — |
| 5. Simulated-mix acceptance per rebuilt panel | #195/#196 acceptance | ½ day | about 60 short jobs per panel, about 1 h on `short` | step 4 |
| 6. Re-type the cassette species with ratio test A | | — | download-bound: about 1,000 *Saccharomyces* runs, about a day on the cluster | steps 4–5 |
| 7. Linked-pair junction panel B, only if A fails the controls | #195 point 3 | 3–5 days (pair-aware counting) | re-download every paired run of the group | steps 4–6 |
| 8. *C. fimbriata* cause on the assemblies | #195 point 4 | ½ day | minutes | — |

## (5) Open questions for the curator

1. **Is an active mating type a wanted output for switching (HO+) species?** Wild *Saccharomyces* are mostly homothallic and often diploid, so the honest answer for most isolates may be "a/alpha or switching" rather than a frequency. If "silent cassettes: not typable" is the final answer for this group, steps 4–7 are not needed for v0.9.0.
2. **Approve or edit the site wording** for each reason in the table in (2).
3. **Relabel *K. marxianus*, *K. phaffii* and *S. kudriavzevii*** from "provisional (homothallic)" to `silent_cassettes` now? Also re-check *H. anomala* and *M. pulcherrima*?
4. **Ratio test A first, with linked-pair panel B only if A fails the controls?** This is recommended, because A needs no new reads where a panel already types.
5. **Truth sets:**
   - May *S. cerevisiae* lab strains (BY4741/2/3) serve as controls, although *S. cerevisiae* stays off the site?
   - The *S. pombe* strains 972, 975 and 968?
   - Which published *C. albicans* MTL zygosity table should be used?
6. **What closes v0.9.0?** Suggested: steps 1, 2 and 8, plus this ruling. Steps 3–7 then move to v0.10 as one "active-MAT typing" milestone.
