# Regression gate judgement: clean run `cleanrun-2026-10` (run-8b1a80d) vs v0.6.0 (run-7c7ed99)

Date: 2026-10-22. Inputs: `results/2026-10_cleanrun_regression/{base,cand,candflags,diff,diff_flags}`.
Diagnostics (targeted re-runs, scripts, logs): `/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10_cleanrun_regression_diag/`
(`runs/<sha>/<genome>/`, `zygo/<sha>/scaffold/`, `bisect_summary.py`, `tasks.tsv`). Bisect worktrees: `.claude/worktrees/reg-bisect-<sha>` (24, detached, unedited).
No repo was committed to; the frozen worktrees run-7c7ed99 and run-8b1a80d were not touched.

## Verdict in one paragraph

**The gate does not pass yet. One code change is responsible for all the Mucoromycota and Zygo call regressions, and it was never approved for detection:**
`0d216c4` "detect: the locus verdict scores the winners of resolved pairs only". It came in inside PR #64 (`wasm-spike`, merged 2026-10-14). Its own message says "the broader effect is to be measured on the comparison panel", and this panel is that measurement.
It costs 2 real Plus calls (R. pusillus JY54, Radiomyces spectabilis), turns 1 Minus/high call into undetermined (Mucor irregularis B50), and lowers 5 Cunninghamella Minus calls from high to medium. Bisection pins each of these to exactly `e26d0e4` (its parent) -> `0d216c4`.
Everything else is explained by approved records, record corrections or curated rules (expected or improvement). A few minor items need a curator look but do not block the run.
One factual correction to the brief: the clean-run flags do change something on this panel (3 loci; see below).

## Preliminary checks

- **Baseline identity.** `git diff 7c7ed99 247c0f5 -- src db` is empty, so the first-parent main history from 247c0f5 to 8b1a80d is a valid bisection axis. There are 56 first-parent merges touching src/db.
- **"Flags change nothing" is not quite true.** A sorted diff of `diff/*/regression_diff.tsv` against `diff_flags/*/regression_diff.tsv` has the same rows and counts. But 3 rows differ in content, all from `--homothallic-policy both_genes` (expected category 5 of the clean-run plan):
  - `GCA_009982855.1_Lizem1` (ascomycota): MAT1-1 -> **undetermined**, `homothallic_candidate`.
  - `GCA_023646395.1_KEW_Xsp1` and `GCA_023647215.1_KEW_Xsp2` (record_sources): MAT1-1 -> **undetermined**, `homothallic_candidate`. This is the right representation for these homothallic Xanthoria loci.
  - Basidiomycota and Mucoromycota diffs are identical between the two runs.
- **Zygo23 truth** (results/2026-09-23_zygo_bar/zygo_truth.tsv): base and candidate both score 23/23 locus and 23/23 idiomorph correct on scaffold and contig inputs. No Zygo idiomorph is wrong; the Zygo changes are confidence or classifier-provenance only.
- **The brief's premise is partly wrong.** It said no Mucoromycota-specific change was approved. Two Mucoromycota changes are in range:
  - **#25** (`b14b421`/`356b9f3`): the Mycotypha combined record, a classifier rebuild and a cleaned negative set, with the gate moving from 99.7 to 97.8 bits. Its BFD effect was measured on 2026-10-04 (`results/2026-10-04_mucoro_bfd_rebuild`: 286 same, 1 changed, 2 gained).
  - **f870eb1** (#64): the sexM/sexP pair rule on profile bitscore. Approved 2026-10-16 (memory: pair-rule-profile-bitscore).
  - **0d216c4** is the third Mucoromycota-relevant change, and the one that was not approved.
  - The reference FASTA, the 5 kb intron cap (#92) and the polish changes are **not** the cause: every one of the Mucoromycota call changes is identical before and after #92.
- **PR B (#133 curation-corrections) is in the candidate.** It was merged at 5d1ef89, although the clean-run plan says PR B "rides a later release". It is the cause of the Xanthoria and many liq-reference model changes (below). Re-word the plan, or accept it explicitly.

## Root cause of the Mucoromycota/Zygo regressions (0d216c4)

0d216c4 removes "set-aside" genes from the profile-HMM classifier input. The set-aside set is built **by gene name** from *every* pair event in the cluster, `{e.loser for e in idiomorph_events_by_cluster[cluster]}` (pipeline.py ~l.2854); superseded tblastn HSPs are dropped from the HSP fallback in the same way. This goes wrong in three ways:

1. **Both genes set aside (JY54, Radspe1).** First-pass identity events on overlapping HSPs go both ways at different positions (JY54 has 6 events: sexM>sexP and sexP>sexM). So both `sexM` and `sexP` are set aside, and their good polished models are thrown away.
   The classifier then falls back to HSP fragments. The MAT-gene gate accepts only model-based scores (or 2 flanks at >= 40%), so the locus is withheld (`mat_gene_gate`), even though the HSP score is 231/320 bits Plus with a 97/177-bit margin.
2. **Wrong winner protected (B50).** The first-pass "winner" is a 48-aa sexP tblastn fragment (37.5% identity, 15% coverage) inside the sexM gene. The polished 2-exon sexM model (77.5% coverage), which the classifier used to score Minus 92 vs Plus 49, is set aside as the "loser".
   The HSP-only verdict then has a margin of 23.4, below 25, so the call becomes undetermined/medium, and a spurious `sexP` joins the gene set. This undoes the 2026-09-26/27 ruling that an identity comparison of short hits must not decide.
3. **Proteome-mode loci lose the classifier entirely (Zygo scaffold inputs).** The winner is a `diamond_proteome` hit with no polished model, and `diamond_proteome` proteins are never classifier input. The only scored protein was the loser's model or HSP, which lies inside the winner's gene.
   Now `idiomorph_classifier: null`. Cunninghamella Minus calls fall back to an identity margin of -3.1 / -7.7, so they become medium and the idiomorph rests on coverage alone. Mucor/Absidia Plus calls keep high on identity margins of 33 / 10, but lose classifier provenance.

Bisection evidence, from `runs/<sha>/` and `zygo/<sha>/` (identical at every earlier point back to 3fe45c9, apart from #25 margin shifts):

| genome | e26d0e4 (= 0d216c4^) | 0d216c4 |
|---|---|---|
| B50 KK076501.1 | Minus/high, model, margin 39.6 | undetermined/medium, hsp_fragment, 23.4 |
| JY54 JAKNBM010000127.1 | Plus/medium, model, 87.3 | withheld mat_gene_gate (hsp 97.1) |
| Radspe1 NW_026251940.1 | Plus/medium, model, 173.9 | withheld mat_gene_gate (hsp 177.0) |
| C. elegans 1390 scaffold_20 | Minus/high, hsp 49.6 | Minus/medium, no classifier, margin -3.1 |
| C. bertholletiae 1376 scaffold_248 | Minus/high, hsp 45.5 | Minus/medium, no classifier, -7.7 |
| Mucor sp. A-21230 scaffold_522 | Plus/high, model 88.2 | Plus/high, no classifier, 33.2 |

Truth: there is no curated record for M. irregularis, R. pusillus or Radiomyces. v0.6.0 called all three (Minus/high, Plus/medium, Plus/medium) on full-length models with HMM margins of 43 to 167 bits. Zygo23 truth confirms the Cunninghamella idiomorphs.
The Mucoromycotina pilot cannot show this. Its 25 genomes contain none of these 3, and they run genome-only with no proteome (the LCG/Jena inputs have no proteome). So **24/25 identical is consistent with the panel**: the defect hits about 3 of 80 genomes with a taxid (about 4%), and proteome inputs only in the Zygo scaffold mode.

## Per-panel tables

Verdicts: **E** = expected (an approved change behaving as designed), **I** = improvement, **A** = needs attention.

### Ascomycota (22 call-touching loci, 30 genomes)

Genomes listed under a commit-range bracket were bisected; the others were attributed from the reference record now used in the cand report (a record new since base, or changed in range). Cause commits are main merges.

| locus | change | cause | verdict |
|---|---|---|---|
| WaNRRLY366 (Hansenula) MATsc JH119216.1, JH119233.1 | span -18/-30 bp, call same | flank/model trim, polish code (not bisected) | E |
| GCA_000498155.1 Shiraia AXZN01000055.1 | MAT1-2-1 279 bp/35% -> 564 bp/42%, 3 exons; margin 14.6 -> 74 | new ref 5022_unknown-2 (Parastagonospora, #36, present before 4ebd2f0) | I |
| GCA_000613145.2 YALI0 MATyl | gain B/high, matb1/matb2 99.7-100% | new MATyl B record 4952_clib122 + scope (#106; bisected 60d7651 -> ff55ca2) | I |
| GCA_000966595.1 Z. brevis LAFY01004296.1 | MAT1-1-1 252 bp/46% -> 888 bp/95%, 3 exons | new Z. tritici ref 1047171_ipo323 (#36) | I |
| GCA_001510565.1 Knoxdaviesia | MAT1-2-1 291 -> 294 bp, 47% -> 52% (ref 61372_suf959, #29); MAT1-2-7 HSP added | #29 F. oxysporum refs (before 4ebd2f0); MAT1-2-7 from Ceratocystis record #98 (seen at c9c0ec6) | E |
| GCA_001557505.1 Reticulascus | **high -> medium**; margin 5.2 -> 2.8 | opposing MAT1-2-1 HSP improves with ref 5022_unknown-2 (**#36**, bisected 3fe45c9 -> 84361bc); MAT1-1-1 ref -> 27337_nbrc6119 (#62/#63) | E (conservative; real locus, MAT1-1 kept) |
| GCA_002924885.1 P. oryzae | span 48-65 kb -> 48-86 kb, call same | new MAT1-2-1 HSP from 5158 Ceratocystis record (#98) | E (span only) |
| GCA_003705895.2 Lipomyces japonicus MAT + MATtub | **call lost** (MAT1-2 high, MAT1-2 medium) | Saccharomycotina subphylum guard, 2fbacb1 (**#106**; bisected 60d7651 -> ff55ca2) | E: pre-registered loss (plan section 6.4); real SLA2/APN2-anchored locus, needs a Lipomycetales family later |
| GCA_004918105.1 A. pullulans | **high -> medium**; margin 26.9 -> 0.4 | MAT1-1-1 now modelled 2 exons/375 bp with ref 1873960_unknown-1 (Pseudocercospora, **#36**, bisected 3fe45c9 -> 84361bc); MAT1-2-1 1 -> 2 exons from #133 liq146 fix | E. Locus carries both MAT1-1-1 and MAT1-2-1, so a near-tie is honest; Aureobasidium is a known problem clade |
| GCA_007556565.1 Elsinoe | MAT1-2-1 1 -> 2 exons, 43% -> 47% (ref -> 115814) | #133 liq80 MAT1-2-1 correction (bisected 0b1f66f -> 5d1ef89) | I |
| GCA_009982855.1 Lizonia (Lizem1) | MAT1-1-1 357 -> 696 bp (62%), MAT1-2-1 313 -> 648 bp (54%); margin 44 -> 91; with flags: homothallic_candidate/undetermined | new Cochliobolus C4/C5 refs (#36; already at 5f01190, so not #132); flag = homothallic policy | I (both complete genes; homothallic class is right) |
| GCA_010015785.1 Sacpr1 | MAT1-2-1 402 -> 366 bp, 34% -> 39%, 2 exons (ref 1873960_unknown-2) | #36 refs | I |
| GCA_022813925.1 Varicellaria | span 6-13 kb -> 6-66 kb (MAT1-1-1 24% HSP at 31 kb, MAT1-1-3 HSP at 66 kb); MAT1-2-1 1 -> 2 exons | span: ipo323/mo44 HSPs (#36, before 0b1f66f); model: #133 | E; span over-extension to 66 kb on weak HSPs is cosmetic (A, low) |
| GCA_022814105.1 Ptychographa | MAT1-2-1 1 -> 2 exons, 45% -> 50% (ref -> 27339) | #133 (bisected) | I |
| GCA_025802115.2 Bertia JAOYNC020017679.1 | **gain** MAT1-2/medium, partial_locus, margin 1.7 (APN2 + MAT1-2-1 2-exon 36% model; MAT1-1-3 HSP lost to it) | first appears between 162d324 and 4ebd2f0 = **#62/#63** V. dahliae records (bisected) | **A (low)**: same genome already has a MAT1-1 call (JAOYNC020110414.1), so this gives an opposite-idiomorph second call on an HMG region with a 1.7-bit margin. Curator: real second idiomorph or HMG paralog? |
| GCA_025802115.2 JAOYNC020110414.1 | span start -114 bp | flank model | E |
| GCA_026023515.1 Claussenomyces | MAT1-2-1 156 -> 168 bp (ref ipo94269) | #36 | E |
| GCA_937854325.1 Physciaceae MAG | **high -> medium**; margin 0.4 | new opposing MAT1-2-1 HSP (26% identity, 18% coverage) from ref 5016_c4 (Cochliobolus, #36; present at 5f01190) | E/A (low): conservative, call kept; the margin rule lets a 26%-identity fragment veto high |
| GCA_964254465.1 Pectenia | MAT1-2-1 653 -> 453 bp, 1 -> 2 exons, 44% -> 48% (ref -> silveira) | #133 liq80 correction (bisected) | I |
| GCA_964255695.1 Leprocaulales | MAT1-2-1 917 -> 615 bp, 2 -> 3 exons, 39% -> 46% | #133 (bisected) | I |
| GCA_022814105.1 / GCA_007556565.1 / GCA_022813925.1 / GCA_004918105.1 (models) | intronic aa removed from the reference, so exonerate now places the intron | #133 (c553380, 1de3509, 77f9f9d): liq80xsp/liq146xsp MAT1-2-1 keeps the published 350-aa start | I |

Core-model changes (14): all are shorter or longer models that now agree with the curated reference structure. Identity rises in 12 of 14. Two cases go down: Reticulascus MAT1-1-1 (39.6% -> 35.6%, a different best reference) and A. pullulans MAT1-1-1 (51.7% over 178 bp -> 39.5% over 375 bp, i.e. twice as long).
Of the 14: 5 come from #133 (the liq protein fix), the rest from new references (#29, #36, #62/63). None comes from #92, #119 or #126. The 5 kb cap hit none of the Ascomycota changes; every Ascomycota bisect pair 221aaf1/295df47 is identical.

### Basidiomycota (5 loci, call never changes)

| locus | change | cause | verdict |
|---|---|---|---|
| S. johnsonii redHD BAABVI010000025.1 | HD1 12.6 kb/2 exons -> 4.6 kb/4 exons, 42.9% -> 43.7%; span -1.3 kb | 5 kb intron cap **#92** (bisected 221aaf1 -> 295df47) | I (chimera removed) |
| Gelatoporia HD CAUPSP010000001.1 | HD2 7.6 kb/3 exons -> 222 bp/2 exons (174 bp, 10% coverage) | #92 (bisected) | E: the old model spanned HD1; the new one is a fragment. Call unchanged |
| C. cinerea CC3 Aalpha NW_003307544.1 | HD2 2961 bp/100%/8 exons over 16 kb (12.7 kb intron) -> 2310 bp/88.2%/6 exons, last exon wrong | #92 (bisected) | **A (low)**: on the reference strain's own assembly the uncapped model reproduced the 987-aa record protein at 100%. CC3 seems to carry a 12.7 kb insertion in HD2 (the record's gene is 3.5 kb in JAAGWA). The cap is right in general; this one model got worse. Call unchanged |
| C. cinerea CC3 PR NW_003307535.1 (+ record_sources GCA_016772295.1) | pheromone 369 bp/2 exons -> 216 bp/1 exon, 100%; span +1.2 kb | B43 record curation 6983433 (**#56**: phb3.1 fix, rcb names, 4 pheromones) | I |
| Serpula PR NW_006763292.1 | receptor 1158 -> 939 bp, 60% -> 65% | #92 (bisected) | I |

### Mucoromycota (80 genomes; calls changed in 3; 11 classifier shifts)

| locus | change | cause | verdict |
|---|---|---|---|
| **GCA_000587855.1 M. irregularis B50** KK076501.1 | Minus/high -> **undetermined/medium**; sexP fragment added to genes | **0d216c4** (#64) | **A: regression** |
| **GCA_023512895.1 R. pusillus JY54** | Plus/medium -> **withheld mat_gene_gate** | **0d216c4** (both genes set aside) | **A: regression (call lost)** |
| **GCF_025331425.1 Radiomyces spectabilis** | Plus/medium -> **withheld mat_gene_gate** | **0d216c4** | **A: regression (call lost)** |
| Mycafr1 NW_026515730.1 | sexP 525 bp/46% -> 741 bp/100%; margin 227 -> 176 | Mycotypha combined record + rebuild (**#25**, bisected 3fe45c9 -> d869561) | I (own record) |
| RhiStoB9770, CunBer175, CunBerB7461, GCA_001683725.1, GCA_024139405.1, Chocucu1, Zycmex1 | classifier margin shifts (-48 to +8 bits), all still decisive (>= 128 bits Plus, >= 88 Minus); CunBer spans +15 bp | classifier rebuild **#25** (bisected) | E |
| GCA_000697015.1 CunEleB9769 | span -15 bp, margin 67.5 -> 63.7 | #25 | E |
| GCA_024139395.1 | loser sexP model 255 -> 246 bp; Minus call same | #59 core-span (162d324) | E |

### record_sources (6 call-touching loci)

| locus | change | cause | verdict |
|---|---|---|---|
| GCA_000150145.1 A. fumigatus A1163 | `homothallic_candidate` -> `mat_locus`; MAT1-2-1 750 bp/100% -> 363 bp/93% fragment | remnant MAT1-2-1 dropped from search refs (83c3939, **#103** remnant status) | I (A1163 is heterothallic MAT1-1) |
| GCA_016772295.1 C. cinerea | pheromone model 369 -> 216 bp | #56 B43 record | I |
| Xanthoria Xsp1 / Xsp2 | MAT1-2 -> MAT1-1 (margin 82 -> 10); MAT1-2-1 1 -> 3 exons; with flags: homothallic_candidate/undetermined | **#133** liq80xsp/liq146xsp MAT1-2-1 correction (bisected 0b1f66f -> 5d1ef89) | E: accepted by curator; homothallic, the record expects both |
| GCA_964255735.1 Xanthocarpia | MAT1-2-1 1 -> 3 exons, 79% -> 82% | #133 | I |
| GCF_000011425.1 A. nidulans | gain MAT1-1/medium (was withheld `modelled_gene_bar+polish_capped`) | polish-cap core-pair fix **#104** | I (record self-call now succeeds) |

### Zygo23 (scaffold+proteome and contig inputs; truth 23/23 both sides)

| locus | change | cause | verdict |
|---|---|---|---|
| scaffold: C. bertholletiae 1376, 1380; C. elegans 1390, 1393, 2310 | Minus/high -> **Minus/medium**; classifier hsp_fragment -> none; identity margin negative | **0d216c4** (bisected e26d0e4 -> 0d216c4) | **A: regression** (confidence; idiomorph still correct) |
| scaffold: Absidia RSA 623, Actinomucor, Circinomucor, M. alternans, Mucor sp. A-14906/21230/21232/21236/25970 (9) | Plus/high unchanged, but classifier input model -> none (decided by identity margin 10-33) | **0d216c4** | A (same defect; no visible call change) |
| contig: Absidia, C. echinulata x2, C. japonica | margin shifts -6 to -15 bits, still >= 176 | #25 rebuild | E |
| contig: C. bertholletiae 1376, 1380 | span +9/+15 bp, margin -4 | #25 | E |

## Needs attention

1. **0d216c4 must be fixed or reverted before the run (blocking).** Evidence: the bisect table above. Each of the 3 Mucoromycota call losses or changes and the 5 Cunninghamella confidence drops starts exactly at 0d216c4.
   Fix options:
   - (a) Revert 0d216c4. This restores v0.6.0 classifier behaviour; f870eb1 stays.
   - (b) Make set-aside **positional** (exclude a loser only where it overlaps its winner, and only when the winner itself has a scored protein).
   - (c) Score `diamond_proteome` winner proteins in proteome mode.
   - Never let a tblastn-fragment "winner" from first-pass identity protect itself from a polished model.
   - Re-run the Mucoromycota panel, Zygo23 and the 288 BFD Mucoromycotina after the fix.
   - **Who decides:** it is a Mucoromycota rule change, so it needs the curator's ruling.
   - **Production impact:** in the full Mucoromycotina run (genome-only, `--phylum`) mechanisms 1 and 2 still apply. The 25-genome pilot just did not sample them.
2. **GCA_025802115.2 Bertia (low).** A new second, opposite-idiomorph MAT1-2 call with a 1.7-bit margin appeared after #62/63. Check the region (MAT1-1-3 vs MAT1-2-1 HMG).
3. **CC3 HD2 under the 5 kb cap (low, model-only).** The reference strain's HD2 is truncated with a wrong last exon. Consider a reference-guided cap exception or a note.
4. **Margin rule vs weak opposing fragments (low).** GCA_937854325.1 dropped to medium because of a 26%-identity, 18%-coverage MAT1-2-1 HSP; GCA_001557505.1 is similar. They are conservative rather than wrong. Track this on the pilot.
5. **Plan bookkeeping.** PR B (#133) is in the tag, contrary to plan section 2. The "flags change nothing" statement should read "3 loci become homothallic_candidate (expected category 5)".

## Recommendation

**Gate: fails on one item. Resolve 0d216c4 first, then the gate passes.**
- Every other call-touching change maps to an approved record, record correction or rule, and the Lipomyces loss is the pre-registered guard cost.
- After the fix, re-run only the Mucoromycota group and Zygo23 of the panel (about 1 h on highclock) and confirm:
  - B50 Minus/high;
  - JY54 and Radspe1 Plus called;
  - Cunninghamella scaffold calls back to high with a classifier verdict;
  - the other Mucoromycota rows unchanged against this cand.
- Ascomycota and Basidiomycota can proceed as they are: 0d216c4 had no effect on any Ascomycota or Basidiomycota bisect pair. Those families have no classifier, so the affected block is never reached.
- Items 2-4 can go on the curator's list without holding the run.

## What I could not verify

- **Whether 0d216c4 was ever signed off.** It is on main via the wasm-spike PR, and I found no curator ruling. The pair-rule memory covers f870eb1 only.
- **Biological truth for B50, JY54 and Radspe1.** There are no curated records. I rely on v0.6.0 model-based HMM verdicts with large margins. The LCG/Jena held-out truth covers neither those genomes nor the panel's BFD rows; only Zygo23 truth was usable.
- **Exact introducing commit for the Ascomycota changes before 4ebd2f0.** GCA_025802115.2 bisected to the #62/#63 window. A. pullulans and Reticulascus bisected to #36. The other new-reference attributions come from the reference record named in the report, not from a re-run.
- **The two small Hansenula span trims and the P. oryzae span were not re-run.**
- **Which pilot genome is the 1/25 that changed.** I did not check, and I did not diff the pilot against v0.6.0 myself.
- **Determinism.** Re-runs of base-equivalent code (3fe45c9) reproduced the baseline calls for the genomes tested, but I did not run a determinism check.
- **Record self-call misses.** The candidate self-call table (`diff/record_selfcall_candidate.tsv`) lists 5 misses: b11221 alpha, 5325 PR_B1, 418459 rustHD, 984962 PR_B1, 1301609 combined. I did not compare it with a baseline self-call table.
