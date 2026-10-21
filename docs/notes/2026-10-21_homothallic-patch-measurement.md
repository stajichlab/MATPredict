# `--homothallic-policy both_genes`: measurement (2026-10-21)

Patched run: 3,181 genomes (sets in `results/2026-10-18_homothallic_patch/sets.tsv`), `--partial-core-policy coverage_profile`, otherwise defaults.
Re-run locally from the committed reports with `compare_patch.py` and `analyze_patch.py` (outputs in `tables/`). `why_not.py` had already been run on the HPCC (6,792 core models, `why_not.tsv`); section 'Why the 298 are not recovered' uses it.

## Result

- **No side effects on legacy output.** 2,595 Dothideomycetes loci appear in both this run and the 10-18 audit run: 0 differ in idiomorph, confidence,
  legacy class or genes. The patch changes only the locus class (`mat_locus` 63 or `idiomorph_gene_only` 5 -> `homothallic_candidate`; 68 loci, 68 genomes).
- **Controls untouched:** 0 of 1,430 baseline heterothallic Dothideomycetes (*Alternaria*, *Parastagonospora*, *Zymoseptoria*, *P. teres*, ...),
  0 of 40 *A. fumigatus*, 0 of 15 *Neurospora*/*Sordaria*, 0 of 47 two-locus genomes, 0 of *Macrophomina* (466) and *Botryosphaeria dothidea* (141).
- **Recall on the 298 conflicting genomes: 29 (10%).** *Didymella* 14, *Stemphylium* 5, other 10. **Not recovered: *Pyrenophora* 0 of 62 and
  *Aureobasidium* 0 of 119**, the two groups that dominate the conflict set. Why they fail is not measured here.

## Problem: false positives on heterothallic *Botrytis*

16 genomes whose BioSample states MAT1-2 (15 *B. cinerea*, 1 *B. fabae*) are re-classed `homothallic_candidate` (idiomorph MAT1-2, medium):
a `MAT1-1-1` gene scores 116-124 bits on the alpha profile and `MAT1-2-1` 221-223 bits on the HMG profile, 1,579 bp apart. This is the failure tracker's F3:
the `MAT1-1-1` is a 61% remnant (see `docs/notes/2026-10-17_botrytis-mislabels.md`), and the alpha floor (min re-classed value 116 bits) does not separate a remnant from a gene.
Of the 17 truth-set genomes re-classed, 1 is correct (*Sclerotinia minor*, stated homothallic or diploid) and 16 are wrong. Adopting the policy as it stands would
call heterothallic *B. cinerea* homothallic.

Other re-classed sets: flip-review 35 of 48 loci (the fused-locus verdicts: Polytolypa, Monilinia, S. nivalis, S. sativa re-classed; Calogaya not);
holdout 2 of 74 (`KEW_Xsp1`, `KEW_Xsp2`, both genes strong at 360-393 bits, 3,354 bp apart; not judged here).

## Why the 298 are not recovered (`why_not.tsv`, first failing gate per modelled core gene)

Across the 298 conflicting genomes: `short` 548 models, `polish_status` 425 (optional genes, not modelled), `ok` 152, `profile_floor` 2. The gate that blocks them is the
**200 aa length floor** (`MIN_PROTEIN_AA`), not the profile gates.

- ***Pyrenophora* (62 genomes):** `MAT1-2-1` is ok (median 236 aa); `MAT1-1-1` is `short` at median **196 aa**, with 204 bits on the alpha profile. Four residues under the floor.
- ***Aureobasidium* (119):** both models are fragments: `MAT1-1-1` median 124 aa (158 bits alpha), `MAT1-2-1` median 103 aa (159 bits HMG). Length 100-125 aa; either a
  real remnant or a modelling failure, not distinguishable here.

Sensitivity (protein-level estimate from `why_not.tsv`, not a pipeline re-run: loci with an alpha model and an HMG model that pass the profile floors, margin 25 and the ORF fraction;
ignores separation and `lost_to`):

| length floor (aa) | loci with both genes | of the 298 conflicting | *Botrytis* MAT1-2 false positives |
|---|---|---|---|
| 200 (current) | 69 | 30 | 16 |
| 190 | 132 | 93 | 16 |
| 150 | 134 | 94 | 16 |
| 120 | 168 | 107 | 16 |
| 100 | 245 | 168 | 16 |

A floor of 190 recovers the *Pyrenophora* set (+63) with no new *Botrytis* false positives, but the *Botrytis* false positives are the same 16 at every floor, so the floor
does not fix them. Below 150 the gain is *Aureobasidium* fragments, which a length-based gate cannot tell from remnants. Whether a 196 aa *P. tritici-repentis* `MAT1-1-1` is
complete is not established here (the 23-of-26 and 24-of-29 curated-length basis for the 200 aa floor is in `homothallic.py`).

## Decision for the curator

Default stays off, and it should not join the clean run as is. Options: (a) require the alpha gene to be full-length against the roster (remnant status, `remnant-status-step2`
removes remnant proteins from references) and re-measure; (b) raise the alpha floor above ~125 bits, which also drops the other 28 true candidates' margin
(min 116, median 189) so needs the distribution by set first; (c) leave the policy off and keep `undetermined`.

## Length floor changed to 190 aa (2026-10-21) and the remnant check (design, not yet implemented)

`MIN_PROTEIN_AA` is now 190 (commit d785fa5); the 3,181-genome re-run is in `runs_out/runs190` on the HPCC (array 30213034), results below.

### What separates the *Botrytis* false positives (measured on the 68 re-classed loci of the 200 aa run, `tables/alpha_hmg_coverage.tsv`)

The detection report already stores, per core gene, the identity and the **coverage of its best reference gene** (`gene_evidence[].identity`, `.coverage`).
The *Botrytis* `MAT1-1-1` model is 209 aa, **92% identical to a curated *B. cinerea* `MAT1-1-1` but covers only 62% of it**: a truncated copy of a known gene, which is what F3 describes.

- **Coverage alone cannot be the gate.** 60 of the 68 loci cover under 90% and 55 under 70%. The divergent groups (*Didymella*, *Phaeosphaeria*, *Lizonia*, *Boeremia*) cover 56-62% at 50-62% identity because
  the references are from other genera, so low coverage there says nothing about completeness.
- **Identity plus coverage does separate them.** Rule: alpha-gene identity >= 85% and coverage < 80% -> not complete, no candidate. It flags 32 of the 68 loci:
  **29 *Botrytis*** (the 16 truth-set MAT1-2 genomes and 13 more with the identical 209 aa model, not in the truth set) and **3 *Sclerotinia* (*S. nivalis* x2, *S. sativa*; 219 aa, 94.9% identity, 65% coverage)**,
  the fused loci the flip review says stay `undetermined`. Thresholds 85/80 versus 90/90 give 32 versus 31 flagged; the result is not sensitive in that range.
- **Kept (36):** *S. minor* (330 aa, 98.5% coverage, truth: homothallic or diploid), `KEW_Xsp1`/`KEW_Xsp2` (352 aa, 100%), and the divergent genera with no close reference.

Limits: the threshold is fitted on this set and the false positives are one lineage; a divergent remnant with no close reference would not be caught (it is indistinguishable from a divergent gene by this evidence).
The rule also relies on the reference roster containing a full-length same-clade gene, which it does for *Botrytis* because of the new MAT1-1 record.

### Implementation plan
Pass `identity` and `coverage` from the gene evidence into `CoreModel` (optional, default none = unchanged behaviour for existing tests), add `REMNANT_MIN_IDENTITY = 85`, `REMNANT_MAX_COVERAGE = 80`
to `homothallic.py`, test on both sides, then re-measure on `runs190` (expect: *Botrytis* 29 -> 0, *Sclerotinia* nivalis/sativa 3 -> 0, *S. minor* and `KEW_Xsp` kept, *Pyrenophora* recovered).

### 190 aa run (3,149 of 3,181 genomes at the time of analysis; `tables190/`, `patched190_reports.tar.zst`)

| | 200 aa | 190 aa |
|---|---|---|
| re-classed loci | 68 | 131 |
| conflicting genomes recovered (of 298) | 29 | 92 |
| *Pyrenophora* | 0 of 62 | 62 of 62 |
| *Aureobasidium* | 0 of 119 | 0 of 119 (100-125 aa fragments) |
| controls re-classed (1,430 heterothallic Dothideomycetes, 40 *A. fumigatus*, 15 *Neurospora*/*Sordaria*, 47 two-locus) | 0 | 0 |
| truth-set false positives (*Botrytis* MAT1-2) | 16 | 16 |

Of 2,563 Dothideomycetes loci in both this run and the 10-18 audit run, none differ in idiomorph, confidence or genes (32 audit loci are missing because their genomes had not finished).
Open: 59 of the 62 *Pyrenophora* are loci where `MAT1-2-1` lost the idiomorph vote to `MAT1-1-3`; they are now labelled `homothallic_candidate` with idiomorph MAT1-1.

### Remnant check, final run (`tables190r/`, `patched190r_reports.tar.zst`; code 9deb233, 26b14e6; 20 homothallic tests pass on the HPCC)

Full re-run, 190 aa floor plus the remnant check: 3,181 of 3,183 genomes (array 30218629, all 128 tasks completed). Two genomes have no report: `GCA_986281625.1_DAH1041FM` (empty FASTA, cannot run) and
`GCA_016165855.1_IBUN_P_ulei_V1` (hit the 3,600 s per-genome timeout; also missing in the earlier runs).

| | 200 aa | 190 aa | 190 aa + remnant check |
|---|---|---|---|
| re-classed loci | 68 | 131 | **99** |
| conflicting genomes recovered (of 298) | 29 | 92 | **92** |
| *Pyrenophora* | 0 of 62 | 62 | 62 |
| truth-set false positives (*Botrytis* MAT1-2) | 16 | 16 | **0** |
| truth-set candidates | 17 | 17 | 1 (*S. minor*, stated homothallic or diploid: correct) |
| controls re-classed (1,430 heterothallic Dothideomycetes, 40 *A. fumigatus*, 15 *Neurospora*/*Sordaria*, 47 two-locus) | 0 | 0 | 0 |
| holdout | 2 of 74 | 2 | 2 (`KEW_Xsp1`, `KEW_Xsp2`) |

- The check removed exactly the predicted loci: 32 removed, 0 added, the other 99 identical in idiomorph, confidence and protein length. The 32 are 29 *Botrytis* (209 aa; 16 truth-set MAT1-2 plus 13 outside the truth set) and 3
  *Sclerotinia* (*S. nivalis* x2, *S. sativa*; 219 aa).
- All 2,595 Dothideomycetes loci match the 10-18 audit run in idiomorph, confidence, legacy class and genes: the patch changes only the locus class.
- Flip-review loci re-classed: 3 (Polytolypa, *Monilinia vaccinii-corymbosi*, `GCA_002162485.1`); the *Sclerotinia* fused loci stay as they were.

Caveats. The 85/80 thresholds were chosen after seeing the *Botrytis* result, so the clean outcome is partly fitted: the held-out test is the next genus where a truncated copy of a known gene occurs, and nothing tests that
yet. *Aureobasidium* (0 of 119) is not recovered and cannot be by a length gate. Whether the 59 *Pyrenophora* loci in which `MAT1-2-1` lost to `MAT1-1-3` should be labelled `homothallic_candidate`, and with idiomorph MAT1-1,
is a curator call. The remnant check covers the alpha gene only; an HMG remnant was not measured.
