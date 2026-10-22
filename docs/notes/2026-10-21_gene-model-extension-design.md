# Gene-model 3' extension and model assessment: design v2, blind-spot audit and evaluation plan (2026-10-21)

Status: **proposal; revised after an independent adversarial review (go-with-changes), then tested on held-out data. Curator decision 2026-10-22: rescue only (section 14). Nothing is implemented.**
**Reading guide:** sections 1-8 are the design as revised after the review; sections 11-15 are measurements and decisions made after it, in order. **Where an earlier section contradicts section 14 or 15, sections 14-15 win:** the false-complete ceiling of sections 9, 12 and 12b was withdrawn after the held-out test (section 13b), and the roster-growth effect was validated on the reserved genes (section 15b).
Tags: **[M]** measured by the author, **[V]** verified by the independent reviewer (how stated in `results/2026-10-21_gene_model_extension/design_review_2.md`),
**[R]** read in code, **[I]** inference or recollection not tested, **[X]** not checked. The reviewer's numbers are reported as the reviewer's; the author has not re-run them.

Base of the code cited: `origin/main` 2d14712, which includes #117 (`homothallic.py`, merged), #119 (roster exon structures, merged at 9ae22d2) and #120. The roster-record `record_version` bumps missed by the #119 merge are in #123, and the homothallic idiomorph follow-up is #124.

Origin: PR #117 found that the 196 aa *Pyrenophora tritici-repentis* `MAT1-1-1` model ends at the edge of the exonerate alignment, not at the gene's end.
A survey of all 3,181 patched-run reports (`results/2026-10-21_gene_model_extension/survey.tsv`) and a prototype chooser led to design v1; the review
(`design_review_2.md`) then reproduced the observation and changed the design as listed in section 10.

---

## 1. The problem

For each gene, the pipeline aligns **all** curated references of that gene with `exonerate --refine region` and with `miniprot`, and keeps the highest-raw-scoring alignment per tool
(`search.py:663-713, 763-771`, **[V]**). `polish.classify()` makes **exonerate canonical whenever it produced a model** (`polish.py:117-118`), whatever its length, stops or agreement with miniprot.
Exonerate (`protein2genome`, local) ends where its alignment ends: the model is the alignment, with no extension to a stop or start (`search.py:908-935`; v1 cited `855-870`, which is the command setup).

Survey, 3,181 reports, in-frame continuation after the canonical model's last exon **[M; reproduced exactly by an independent implementation on all 4,270 rows, [V]]**:

| 3' state of the exonerate model | MAT1-1-1 (n=1,721) | MAT1-2-1 (n=1,781) |
|---|---|---|
| internal stop in the model (as exons are concatenated; see below) | 142 | 350 |
| next codon is a stop (complete) | 182 | 444 |
| stop within 1-29 aa (edge slop) | 281 | 82 |
| open 30-99 aa | 184 | 47 |
| open >= 100 aa | 924 | 855 |

64% / 51% have >= 30 aa open **[V]**. Read these with three corrections:
- the denominator includes internal-stop rows; among clean models it is 70% / 63%;
- the rows are dominated by clonal lineages (Macrophomina 490 rows, Aureobasidium 357, Botryosphaeria 280, Alternaria 267; 223 genera): report per lineage, not per genome;
- **the "internal stop" counts are overstated.** The report does not serialise exonerate `cds_blocks`; 198 of 3,502 canonical models (5.7%) carry frameshifts and **188 of the 492 "interrupted" models (38%) are frameshift models** whose stops come from reading through the frameshift.
  Rebuilding the blocks by re-running exonerate on those 198 (same options, 170 reproduce the reported exons exactly), 87 of 190 become clean (46%) **[V]**. Every "interrupted" statistic in v1 must be redone with `cds_blocks`.

Other corrected facts: the *Botrytis* remnant's in-frame stop is at residue 10, the model is 5'-truncated and 3'-complete (209 aa after the stop, no M, upstream stop 18 aa before the model; 29 of 30 *Botrytis* genomes with a stop) **[V]**;
exonerate models do not end in a stop in 0 of 3,296 non-frameshift core rows (but 2 MATA1 yeast models do); miniprot alternates end in a stop in 27% (v1 said 19%) **[V]**;
the extended length is plausible: `MAT1-1-1` clean model + extension has median 356 aa (IQR 344-378) against curated references of median 359 (10-90%: 234-482). This tight reference-like distribution is the strongest single piece of evidence for the extensions **[V]**;
`_translate_model` (`pipeline.py:851-880`) returns the longest stop-free stretch over **three frames**, so the existing ORF-completeness test already passes 37 stop-interrupted models via a wrong frame **[V]**.

---

## 2. Definitions

| term | definition |
|---|---|
| `clean` | no internal in-frame stop in the **correct reading frame** (below), using the `cds_blocks` translation for frameshift models |
| `3'-closed` | the codon after the CDS is a stop |
| `3'-open` | not 3'-closed, and the frame continues >= `EXT_MIN_AA` before the next stop |
| `starts M` | first codon is a start |
| `n_term_truncated` | no M and an in-frame stop within `NTERM_STOP_WINDOW` aa upstream with no M between. **Window-dependent:** it flags 4.9 / 6.1 / 6.3 / 15.7 / 17.1% of models for windows of 10 / 20 / 30 / 60 / 100 aa **[V]**, so the window must be fixed on a truth set, not guessed |
| `complete` | clean, 3'-closed, starts M or flagged `n_term_unknown`, aligned coverage of the best reference >= 80% (profile envelope coverage once built) |
| `partial` | clean, aligned coverage 50-80% |
| `remnant` | clean, identity >= 85% to the best reference and aligned coverage < 80% (existing rule), or `n_term_truncated` with coverage < 80% |
| `interrupted` | internal stops remain after the `cds_blocks` correction, with the first-stop position and the number of reference columns that align after it |

Principle, refined: a clean, 3'-closed model is preferred to a shorter one only on **evidence of the right class** (section 3.4); **a model with internal stops is never promoted;** a clean model that is a truncated copy is flagged, not called complete.
"Longer without stops is better" is a ranking default and not a correctness test: a 419 aa miniprot model with 6 internal stops, a clean *Botrytis* remnant, and stop-interrupted *Aureobasidium* models whose reference aligns across the stops are all counterexamples.

---

## 3. The algorithm (deterministic code; no model in the loop at run time)

Ordinary Python in `detect`; inputs are the genome contig, the model's exons or `cds_blocks`, the genetic code and the curated proteins in the database. All steps are pure functions of those, hence reproducible and unit-testable.

### 3.1 Where it runs, and what Phase 1 does *not* do

New module `src/MATPredict/detect/model_selection.py`, called once per polished gene right after `polish.classify()` (`pipeline.py` ~2722-2726) and before `GeneEvidence` is built.

**Phase 1 does not swap tools and does not move the 5' end.** The canonical model stays exonerate's; the new step only *assesses* it and, when the evidence allows, proposes a 3' extension of that same model. Miniprot's model is reported with an `independent_support` flag.
(v1's chooser, as specified, made miniprot canonical in 56% of loci (`mini+ext` 1,157 of 3,147) and moved 5' ends in about 78% of those, contradicting its own statement that the 5' end is not moved. **[V]**)
A tool swap is a separate, named, replayed policy (section 9, stage 3b).

```
assess(outcome: PolishOutcome, contig_seq: str, genetic_code: int | None,
       references: Sequence[RefProtein], neighbours: Sequence[Span]) -> ModelAssessment
```

### 3.2 Frame

- **Exonerate, no frameshift:** frame 0 of the concatenated exons. The summed exon length is a multiple of 3 in 3,299 of 3,299 rows **[V]**; no frame search is needed or allowed.
- **Exonerate with frameshifts:** translate the **`cds_blocks`** (the aligned blocks, which leave frameshift bases out; `pipeline.py:837-850`). `PolishModel` already carries `cds_blocks` (`polish.py`), the report does not serialise them: **carry them into the shadow output.**
  Extension for these models starts from the **block CDS** end (the exon-concatenated CDS is 1-16 bp longer than the block CDS, so walking from it starts in the wrong frame) **[V]**.
- **Miniprot:** the GFF phase (`PolishModel` does not keep it: add `phase: int | None`); 6.4% of miniprot alternates have a CDS length not divisible by 3 (170 of 2,650) **[V]**.
- **No "frame with zero internal stops" search.** v1's fallback manufactures clean models from interrupted ones: of 492 exonerate canonicals with frame-0 stops, 16 have a stop-free *different* frame, in 3 of which the tool frame is certainly 0 (a spurious "rescue" that would extend 42-58 aa fragments by 7-107 aa), and 40 of 3,010 clean models are clean in two frames **[V]**.
- A model that is not clean in its tool frame is `interrupted` and is **not extended**.
- `genetic_code` null (6 reports): no extension, flag `genetic_code_unknown`. Only codes 1 and 12 occur in the reports (same stops); other stop sets need unit tests only **[X on real data]**.

### 3.3 3' extension (an unscored, flagged step)

```
cds   = block CDS (frameshift) or spliced exons, in the tool frame, 5'->3'
rest  = len(cds) % 3                      # leftover bases are the first bases of the next codon
walk downstream on the gene's strand from the CDS 3' end, codon by codon (leftover first), until
      a stop of `genetic_code`            -> closed, ext_aa = codons read
   or EXT_MAX_BP (3000) / contig end      -> unclosed ('contig_edge'; no claim; v1's "no stop in 3000 bp" rows were all contig edges [V])
   or the start of a neighbouring modelled CDS on the same strand (ext_overlaps_neighbour; stop there, no claim)
new last-exon end = old end + 3*ext_aa    (stop codon excluded, as exonerate does)
```
Strand, multi-exon order, partial codons and stop table were independently re-implemented and agree on 4,270 of 4,270 rows **[V]** (both implementations took the report's exons as the CDS, so frameshift handling was not covered; section 3.2).

**Extension is not scored by BLOSUM.** A local alignment score cannot decrease when residues are appended (`score(exo+ext) >= score(exo)` always); for tails >= 100 aa, 43% gain <= 10 raw units, median per-residue gain 0.1 **[V]**. So "ties within a tolerance go to the longer" reduces to "always extend".
v1's presentation of the extension as a scored choice is withdrawn; the decision is a rule on **evidence class** (3.4). (`TIE_TOL = 3` raw BLOSUM62 units is also effectively an exact tie, and the prototype's `round(score/3)` binning is not a tolerance; do not port either.)

### 3.4 Evidence tiers for an extension

The open frame is weaker evidence than v1 said. At the real genome GC, iid arithmetic gives P(open >= 100) = 0.4-1.4% and P(>= 200) <= 2e-4, but that is the **wrong null**: the open run in the other frames after a real stop is >= 100 codons in 1.4% (frame +0), 16% (+1), 5% (+2) of complete models, and in 78% (Alternaria frame +1), 19% (Pyrenophora +1); the +1/+2 shifted frames of the extension region are open >= 100 codons in 20% / 4% of all models and 22% / 32% in Aureobasidium **[V]**.

| ORF-only length | evidence rule |
|---|---|
| < 30 aa | not an extension (edge slop) |
| 30-199 aa | **requires independent support**: miniprot's own model ends in a stop at the same codon, or tail similarity beats the nulls, or profile envelope gain, or spliced RNA reads |
| >= 200 aa | acceptable on the open frame alone (<= 0.5% in all but one clade under the empirical null) |

Independent support available today: miniprot's **own** stop coincides with the exonerate extension's stop in 200 of 205 rows **[V]**, but only 10% of extensions have one (Ascochyta, Aspergillus, Didymella, Alternaria, Pyrenophora, Dothistroma), and **none in Macrophomina, Botryosphaeria, Aureobasidium, Hortaea, Venturia, Phyllosticta**, the clades that dominate the extension count.
The tail aligns better than 10 amino-acid and codon shuffles, +1/+2 frames and reversed tail in only 37% of extensions (Macrophomina 4%, Botryosphaeria 0%, Venturia 0%; Pyrenophora 94%, Didymella 100%, Alternaria 80%) (Smith-Waterman BLOSUM62, best of all curated references) **[V]**.
So for the largest clades the extension rests on ORF length and on the plausible-total-length argument (section 1) alone. The proposal must decide **in advance** what promotes those lineages (Macrophomina, 490 rows, all open) to `complete`; this note's answer: total length inside the curated range plus a >= 200 aa open frame, flagged `ext_evidence=orf_only`, and never reported as independently confirmed.

### 3.4b Rescue use: extension is for loci whose tail we cannot align (curator ruling 2026-10-21)

In production the extension has a second job beyond tidying models: **rescuing loci whose C-terminus evolves too fast to align to any curated reference**, where the homology alignment ends early and a full-length alignment is impossible (the reviewer found no alignable tail for *Pyrenophora*, *P. teres* or *Aureobasidium* by exonerate's own scoring, section 5). So:
- **Absence of tail homology is expected and is not evidence against an extension.** Tail similarity can only raise confidence; it is never required.
- An extension is always **computed and reported** with its `ext_evidence` class (`orf_only`, `orf+independent_stop`, `orf+tail_similarity`, `orf_weak`), so consumers can decide how much weight it carries.
- A gene counts as **`complete_orf_supported`** (usable for the completeness tests that decide whether a locus is called) when: the homology alignment covers the conserved domain (alpha box / HMG box present, profile bits above the existing floors), the frame stays open to a stop in the tool frame, the total length lies inside the per-order curated range (and under the 1.25 x guard), no neighbour overlap, and no `GT`-donor ambiguity at the model end. For >= 200 aa extensions that is sufficient; for 30-199 aa the gene is **`provisional`**: it is **not rejected** and it is flagged, so a consumer may accept it at lower confidence rather than discard the locus.
- What this protects against: calling a gene complete on a chance open frame (5-20% of unrelated frames are >= 100 codons open); it does not require homology, it requires the *combination* of domain + length range + stop + guards.
- **How rescue is measured** (stage 0b, section 9): take curated records whose deposited CDS is known, withhold the record, and give exonerate only the N-terminal part of a related reference (a simulated fast-evolving tail: the alignment ends where homology is made to end). Then check whether extension to the next in-frame stop recovers the deposited gene end, and how often it overshoots. This measures rescue accuracy against truth **without relying on tail homology**, at several truncation points (30%, 50%, 70% of the reference), per clade.

### 3.5 Guards (all new in v2)

1. **Neighbour guard:** stop at the next modelled CDS on the same strand; flag `ext_overlaps_neighbour`. 24 of 2,010 extensions (1.2%) overlap another reported gene (21 same-strand: *MAT1-1-1* into *MAT1-2-1* in *Didymella*, *Phoma*, *Elsinoe* ...); a lower bound, since non-MAT neighbours are not in the report **[V]**.
2. **Length guard:** total protein <= 1.25 x the longest curated record of that gene; 52 (5.8%) of extended `MAT1-2-1` exceed the curated 90th percentile, and the three homothallic loci below reach 391-454 aa **[V]**.
3. **Donor-site flag:** whether the two bases after the model end are `GT` (a 5' splice site), standard flag. Overall 4-7% (no enrichment, as expected for a model that runs out of similarity inside an exon), but **Aureobasidium 24.2% (n=194), Friedmanniomyces 39.6%, Phyllosticta 12.7%** versus Macrophomina 1.6%, Botryosphaeria 0%, Pyrenophora 0% **[V]**:
   in about a quarter of Aureobasidium models the exonerate alignment ends where an intron begins, and an extension would read intron as codons. **Do not extend a model that ends at a `GT` donor on ORF length alone**; require independent support. This matters for Aureobasidium, the second motivating clade and the one with the weakest ORF-only evidence.
4. **Contig-edge:** no claim when the walk reaches the contig end.

### 3.6 Scoring a model (assessment, not selection)

Computed for the canonical model and, for comparison, the miniprot model: `terminal_stop`, `internal_stops` (positions, from the correct frame), `starts_M`, `n_term_truncated`, `ext3_aa` and `ext_evidence`, `ext_overlaps_neighbour`, `donor_GT`, `aligned_coverage` (reference residues inside the local alignment span / reference length, **not** the existing length ratio), profile envelope coverage (pyhmmer against a per-gene profile built from the curated alignment: median 0.355 -> 0.603 on extension, gain >= 0.1 in 85% **[V]**), identity, and the class `complete | partial | remnant | interrupted`.
The raw BLOSUM62 local score (`helixer_rule.query_score`) is kept as a recorded number and used **only** for the exonerate-vs-miniprot comparison in 3.7. Re-scoring with a different reference changes that comparison in 3.2%, 3.1%, 2.2% of rows for the 2nd, 3rd, 4th best reference (all tool swaps, none extension decisions) **[V]**, so the "one reference per target" concern is minor.
Expected length range per order: 10th-90th percentile of the curated records of that gene in the order (global range for orders with < 3 records). A flat aa floor is wrong in both directions (curated *Pezizales* `MAT1-1-1` are 140-171 aa).

### 3.7 Exonerate versus miniprot: a named rule, not a silent swap

Phase 1 records, never swaps: `model_ambiguity` when the two structures differ by more than a tolerance at either end, with the end and size (miniprot is longer **and clean** in about 1,000 rows per gene, and 5'-longer by >= 30 bp in 60% of clean pairs **[V]**).
Known tool limits that a rule must respect: **miniprot places only GT-AG introns (99.8%)**; the *P. tritici-repentis* `MAT1-1-1` intron (53-56 bp, in all 76 *Pyrenophora* `MAT1-1-1` models) is **GC-AG**, so miniprot codes it as CDS, which is the "419 aa, 6 internal stops" model (v1 called this "miniprot skipping an intron"; the arithmetic checks: 1,259 bp = 588 + 53 + 615 + 3 **[V]**).
Exonerate's own introns are 93.4% GT-AG, 3.3% GC-AG (Pyrenophora 40.8%, Fulvia 23.5%) and 1.3% TA-AG (**Aureobasidium 16% of introns have a TA donor, which does not occur in nature**, so those introns are mis-placed) **[V]**. A splice-site check is therefore a concrete deterministic improvement (stage 3b), not a hypothetical.

### 3.8 5' end

Reported only (`starts_M`, `up_to_stop`, `up_first_M`, `n_term_truncated`); not moved in Phase 1. It is the **larger residual**: 79% of models lack M, and in 76% of those an in-frame M exists upstream before the stop; miniprot is longer at the 5' end in 60% of clean pairs; for `MAT1-2-1` the problem is mostly 5' (the *P. tritici-repentis* `MAT1-2-1` exonerate alignment starts at reference residue 109 of 343, is 3'-closed, and miniprot aligns from residue 9 to 335 aa clean) **[V]**.
"The problem is mostly 3'" holds for `MAT1-1-1` in Pleosporales and Botryosphaeriales; it does not hold in general. The `MAT1-2-1` median model is 121 aa (length ratio 31%), a worse problem than `MAT1-1-1`'s 189 aa (49%) **[V]**.

### 3.9 Staging and the policy flag

`--gene-model-policy {legacy, shadow, choose_extend}`; default `legacy`. `shadow` writes a `model_assessment` into `GeneEvidence` and the report and changes no call. Risks the review raised, with the mitigation adopted:
1. A new report field makes every replay diff light up (`detect/regression.py:134-150` compares model exons and length): give the regression harness a **whitelist mode** that ignores `model_assessment`; byte-identity of `legacy` guarantees only the old path.
2. Shadow must use `cds_blocks` or it silently measures the wrong thing for 188 models.
3. A three-way flag plus per-consumer flags multiplies states: cap the matrix; consumers read **one** field, `model_assessment.complete`, and are enabled together in a **named bundle**, each bundle replayed once.
4. Shadow output is the only new evidence and will be read as validation: pre-commit an **acceptance table** (section 9) before it is opened.

### 3.10 Code to be written (for review)
`model_selection.py` (`find_frame`, `extend_3p`, `donor_flag`, `neighbour_guard`, `assess`, constants in one place with their evidence); `phase` on `PolishModel`; `cds_blocks` serialised in shadow output; `GeneEvidence.model_assessment`; `tests/detect/test_model_selection.py`:
both strands, multi-exon, partial codon at a junction, contig edge, genetic code 12 and a synthetic non-standard stop set, a frameshift model, the *Botrytis* case (clean remnant stays `remnant`), the miniprot GC-AG case (stops rejected), complete genes not extended, a `GT`-donor model not extended on ORF length, the neighbour guard, determinism, and `legacy` reports unchanged.

---

## 4. Blind-spot audit

Tests that treat the canonical model as the gene. Strategy for all: **give each the right signal and switch it behind a policy flag after a replay; do not edit thresholds.** The reviewer verified the file:line citations except where marked.

| # | test | where | failure | strategy | replay that must pass |
|---|---|---|---|---|---|
| 1 | homothallic completeness | `homothallic.py:48-50,89` | passes vacuously for edge-ended models; the aa floor measures aligned extent; passes some stop-interrupted models through a wrong frame; 1,684 of 2,418 `short`/`orf_truncated` rows would reach >= 190 aa clean after extension (upper bound) | use `complete`; drop the aa floor; per-order range | 3,181 set: *Pyrenophora* and *Aureobasidium* recovered under the evidence tiers, *Botrytis* excluded, controls at 0; **the 3 loci (Macan1, ASM1692080v1, Decga1) whose HMG gene extends by 192-219 aa into the alpha gene on the same strand, taking separation to -317, -251, -224 bp, must be explained (the extension is itself suspect there)** |
| 2 | homothallic remnant check | `homothallic.py:81-84` | uses the length-ratio `coverage`, which extension raises | apply to the chosen model's `aligned_coverage`; add `n_term_truncated` | *Botrytis* (29) and *Sclerotinia* (3) still removed |
| 3 | `full_length_models`, `HOMOTHALLIC_MIN_MODEL_COVERAGE = 0.5` | `pipeline.py:539-571` | edge model reads 32-60%, straddling 0.5 | `aligned_coverage` | sample of 50 loci changing class |
| 4 | partial-core vote trigger | `partial_core.py:48-49,112-121` | edge model reads partial and can fire a flip | **measured effect is small, not the highest risk v1 assumed:** of the 47 loci with the coverage-rule note (41 flipped, 6 undetermined), 38 of the 41 flips have a winner the extension cannot change (30 *Botrytis*, 1 *Botryotinia*, 1 *Sclerotinia*, 1 *Stromatinia* have an internal stop and are not extended; 5 *Ophiostoma* extend by 1-4 aa); of the other 3, only *Sporothrix* (22 -> 87% by length ratio) crosses 80% and would lose its flip (**1 of 41**); three `undetermined` (*Calogaya*, *Leptosphaerulina*, *Polytolypa*) would cross 90% and revert to the bitscore vote. "41 of 46 right" is dominated by *Botrytis* pseudoreplication (31 of 41), so it does not test the rule's value **[V]** | replay exactly those 47 under (a) length-ratio on the extended model, (b) `aligned_coverage` on the extended model, (c) `aligned_coverage` on the unextended model; require the 38 unchanged flips and *Sporothrix* unchanged under (b); curator adjudicates the three `undetermined`; add a negative-control set | the 47-locus replay with those expectations |
| 5 | confidence cap | `confidence_cap.py:96` | partial winner capped | chosen-model flags (cap is off by default) | rerun the `confidence-cap` evaluation |
| 6 | strong single gene | `strong_single_gene.py:57,79` | a correct edge-truncated lone gene fails the 90% bar; the *Cochliobolus* case | `complete` | the 12-locus check (decided to follow the clean run) |
| 7 | ~~idiomorph resolution coverage~~ | `idiomorph.py:55` | **removed:** the recorded `winner/loser_coverage` are tblastn hit coverages decided **before polishing** (`basis: first_pass_identity`); a longer model cannot change them; `idiomorph.py:55` says coverage "is recorded but not yet used". The 59 *Pyrenophora* MAT1-1-3-over-MAT1-2-1 resolutions are counted correctly but are not fixed by this work | none | — |
| 8 | `_translate_model` | `pipeline.py:851-880` | wrong frame can win | strict translator for the new module only; classifier scoring keeps the old function until a replay shows no change. Note `pipeline.py:2850-2870, 3050, 3081` translate **both** models for the profile classifier, so profile scores are already a max over tools; extending the protein would change them, which is **not** "no change" | profile bits unchanged on the 3,181 set |
| 9 | helixer rule | `helixer_rule.py:2-25` | core genes keep the homology model | out of scope; held-out test pending | — |
| 10 | `mat_gene_gate` | `:68,87` | domain score | no change | bits unchanged |
| 11 | `two_idiomorphs`, `split_locus`, `_short_orf_genes` | `two_idiomorphs.py:58`, `split_locus.py:117`, `pipeline.py:650-700` | not length-blind | no change | — |
| **new** | **coordinate consumers** | `homothallic.py:107-108` (separation <= 0 discarded as "one region seen twice": 99/99 reproduced, **3 of 99 flip**); `core_span`/`supported_span`/`beyond_core_bp`; `contig_edge_distance`; the 500 bp `EDGE_BP` rule `partial_core.py:51`; `locus_merge.py`; `span_exceeds_plausible_bound`; `split_locus._edge_distance` | an extended model moves spans and edges | include in every replay; the neighbour guard (3.5) is the first mitigation | the three loci above |
| **new** | **model output and replay** | `report.py:213-214,228,446-453` (GFF `coverage=`, YAML exons), `report/svg.py:59-60,224-255`, `report/html.py:337-370`, `detect/regression.py:134-150` (exon and length diff), `track_report.py:44-68`, `tracks.py:66` | changed models show as "model changes" | whitelist mode in `regression.py` (3.9) | — |
| **new** | **benchmark** | `benchmark.py:379-440` re-translates reported exons with `codon_start=1`, `table=1`, `to_stop=True` and scores against curated proteins | hard-codes table 1; changes benchmark and holdout identity/coverage | use as the truth-set scorer (good), after removing the table-1 assumption | — |
| **new** | **polish window** | `pipeline.py:1731-1746` | an extension can end outside a clipped window | compute from the contig, not the window | — |
| 12 | `scoring.py`, `holdout.py` | — | no references to `.exons` or `coverage` (grep) | none | — |

Order of work: 1 and 2 (the motivating case, with the evidence tiers), then 3 and 8 (no vote change), then 5-6 (off by default), then 4 under its own replay.

---

## 5. Is a better HMM gene model needed?

Revised answer after the review: for the **3' end of a clean, non-frameshift, fully aligned model, extension is enough.** But the residual is **not small and not mostly 3'**: frameshift models (5.7% of models), the 5' end (79% without M), and extensions with no independent evidence in the dominant clades. A gene-finder-class candidate generator is the standard answer to the **5'** problem and to splice-site errors, so it is the likely next candidate generator rather than a late option; **decide after a 5' truth set (section 7), not after the 3' extension.** Do not build an in-house gene finder.

- **Profile HMMs already shipped** are domain models: they do not know where a protein ends. Their use here: **envelope coverage as a completeness signal** (cheap; does not discriminate real from spurious tails better than Smith-Waterman, since both are driven by the reference set).
- **AUGUSTUS with protein-profile extension (PPX)** [R]: the most credible "better gene model" for single loci (supplies start, stop and splice signals); needs trained species parameters (shipped models exist for a handful of model fungi [I; the installed AUGUSTUS species list was not checked, AUGUSTUS was not installed in the reviewer's environment]; none for Macrophomina, Aureobasidium, Pyrenophora) and a run-time dependency. A *candidate generator*, never canonical.
- **GeMoMa**: needs reference genomes and annotations per reference species; the database holds proteins and CDS coordinates, so a poor fit unless about 170 curated genomes are fetched. **Skip.**
- **GeneWise with a profile**: unmaintained, behaves like exonerate. **Skip.**
- **Helixer**: independent second opinion for 3' and 5' ends where available; not run here **[X]**; the core-gene use needs the pending held-out test.
- **exonerate `protein2genome:bestfit`:** exists in the installed exonerate 2.4.0, "with which the entire protein is included in the alignment", **only available with exhaustive alignment (`-E yes`)**; incompatible with `--refine region` and heuristic seeding, and slow **[V]**. **Result: it could not be run.** On exonerate 2.4.0 with `-E yes`, a 1.7 kb window on the Ptr134 `MAT1-1-1` locus hit a 40 min limit with no output; six loci on 1.6-1.8 kb windows produced nothing in 1.5 h; and a reduced probe (only the C-terminal half of the reference, 124-228 aa, against only the 650-900 bp extension window) timed out at 25 min on all four loci, including the smallest (124 aa on 650 bp). The production setting (`--refine region`) takes 0.1-0.4 s on the same inputs, so `bestfit` is at least 3-4 orders of magnitude slower and **not feasible per locus** **[V]**. By-product: local exonerate of the C-terminal reference half against the extension region gave **no alignment** for *Pyrenophora*, *P. teres* and *Aureobasidium* (the *Parastagonospora* control aligned at 100% over all 124 aa), so the tails of the lead clades are not alignable to the nearest reference C-terminus by exonerate's own scoring; a forced `bestfit` alignment there would exist by construction and would need a null. Not run: `bestfit` on any of the 21 loci, and any other global mode.
- **Several references** is already the status quo (all curated references of the gene are aligned): running exonerate against each reference separately closes the 3' edge in **0** of 5-19 references at each of 10 extension loci; extended totals are nearly reference-independent (Ptr 403-409 aa, Aureobasidium 360-386). So **the edge is about where homology ends, not which reference was used** **[V]**.
- **`--refine region` vs `--refine full`:** identical for *Pyrenophora*, *P. teres* and the controls; for 3 of 4 *Aureobasidium* `MAT1-1-1`, `--refine full` gives 192-224 aa instead of 124-149 (+68-75 aa at the 3' end), so `--refine region` is part of the edge problem there. `search.py` notes `--refine full` crashes deterministically on some windows; a crash fallback chain full -> region -> none is cheap and was not in v1 **[V]**.

**Splice-site scoring is a concrete, small, deterministic improvement** (section 3.7): the GC-AG and TA-donor findings are the main thing a trained gene finder would add here. Do it with a position-weight check on exonerate's and miniprot's introns after the extension step is shadowed.

---

### 5b. Can miniprot be made to place non-canonical splice sites? (tested 2026-10-21, miniprot 0.18-r281)

Yes, through four controls, **measured on the *P. tritici-repentis* Ptr134 `MAT1-1-1` locus** (a 3 kb window; query: the curated *Bipolaris* C5 `MAT1-1-1`; the intron is GC-AG, transcript-sense `GCGA...CCAG`, 53 bp) **[M, `results/2026-10-21_gene_model_extension/miniprot_splice_tests/`]**:

| setting | result |
|---|---|
| default (`-j 1`, general model) | 1 exon, 293 aa, **5 internal stops** (the intron is coded as CDS) |
| `-j 2` (vertebrate/insect) | same as default |
| `-j 0` (no splice model), or `-C 0` (ignore splice signals), or `-C 0.3` | 2 exons, 276 aa, 0 stops, but the intron boundary is **2 bp off** (`GA..GC`) |
| `-J 15` (intron-open penalty 29 -> 15) or `-F 50` (frameshift/stop penalty 23 -> 50) | 2 exons, 282 aa, 0 stops, a **different, canonical GT-AG** 35 bp intron |
| **`--spsc FILE`** scoring the junction (score 14) | **2 exons, 276 aa, 0 stops, intron exactly exonerate's: `GC..AG`, 53 bp, 1397174-1397226** |

- The manual says the general model "slightly prefers" `G|GC...YAG|` and `|AT...AC|`, so GC-AG and AT-AC are not forbidden; the default failure here is the surrounding-sequence penalty plus the stop-codon cost, not the dinucleotide alone.
- `--spsc` is a TAB file: `contig  offset  strand  D|A  score`. **Offset convention, found empirically (only this combination worked among 12 tried):** with the forward-strand intron at 1-based `[a, b]` in the target, the **acceptor offset is `a-1`** and the **donor offset is `b`** (0-based junction coordinates on the forward strand) for a minus-strand gene; **both** sites must be scored (donor-only and acceptor-only change nothing). The plus-strand convention was not tested and must be checked before use.
- `--spsc` also raises `-J` and `--J2` by 10 unless set, scores outside the file get `--spsc0` (-7), and scores are capped at `--spsc-max` (14). So **a splice-site table is a way to inject a clade-specific splice model into miniprot deterministically**: score every candidate donor (GT, GC, AT...) and acceptor (AG, AC) in the window from a position-weight matrix trained on that clade's introns.
- **Limits of this test:** one locus, one query, the junction was given by the answer (an oracle). It shows the mechanism works, not that a trained splice-score table would find the right junction without knowing it. The test of that is stage 3b: score all candidate sites with a clade PWM, run miniprot with `--spsc`, and compare with the deposited CDS of the curated records and with exonerate. Cost: one extra miniprot call per gene, the window is already cut.
- `-j 0`/`-C 0` alone are not a fix: they remove the penalty everywhere and shifted this intron by 2 bp, which is a frame-compatible but wrong boundary.

## 6. Where the problem comes from

Mostly **design** (exonerate is always canonical, then tested with a vacuous completeness definition), amplified by exonerate's local-alignment edge; the reference is a minor factor (above); miniprot helps (clean and longer in about 1,000 rows per gene) but is itself edge-limited (628 of 1,349 `MAT1-1-1` alternates end >= 30 aa before a stop), cannot place GC-AG introns, and is much more reference-sensitive than exonerate (single-reference miniprot models of one locus range 82-396 aa), so the pipeline's all-reference selection is doing real work.
Where miniprot reaches its own stop (all 4 *Pyrenophora* loci checked, 2 of 2 *P. teres*, 8 of 8 controls, **0 of 4 *Aureobasidium***), the motivating *Pyrenophora* clade does not need the extension to be right; *Aureobasidium* needs it most and has the least support.

---

## 7. Truth set and training

v1's primary truth (annotated proteins in BFD genomes) is **circular where it matters**: most annotations in this ecosystem are funannotate/Augustus/exonerate-evidence models that end at the same edges, so agreement is a coin toss and disagreement is ambiguous. Use, in order:
1. **The 171 curated records in `db/`** (26 `MAT1-1-1`, 29 `MAT1-2-1`, plus flanks): deposited NCBI CDS with exon structure, independent of the pipeline. The existing holdout machinery (`benchmark.py`, `docs/holdout-benchmark.md`, radii record/species/genus/family/order) already withholds a record and re-detects it; **add one metric: exact protein and exon structure, and 3' and 5' end offsets, of the detected model against the deposited CDS.** n is about 55 for the two main genes; the divergence gradient is built in.
2. RNA-seq-supported annotations only where they exist (the BFD tree has `do_annotation*` and PASA directories **[X: not opened]**); label the support (full-length cDNA/PASA vs ab initio).
3. UniProt/RefSeq proteins with manual or transcript-supported annotation (a few dozen across Pleosporales, Sordariomycetes, Eurotiales): length class and C-terminus checks, not exon truth.
4. Cross-species conservation of the extended C-terminus (e.g. `...NDAFRVGADEDATLPTFDGANNA` in *P. tritici-repentis* and *P. teres*): consistency, not truth.
5. Spliced RNA reads across a donor right after the model end: direct evidence of the intron failure mode.

Splits: leave-one-genus-out is the minimum, but the data are clonal (2,756 genomes, 223 genera): **weight per lineage cluster** (the BFD `do_ANI` machinery), report per clade, hold whole clades out for threshold selection (Macrophomina is a good one: largest, no independent corroboration), and fix flagging thresholds before opening the truth. Report 5' and 3' errors separately; primary harm metric: **wrong promotion rate** (a model extended or promoted that disagrees with the truth).

**A learned chooser is not warranted.** Four candidates, fewer than ten informative features, and tens of independent lineages: a logistic regression would overfit lineage identity (genus is collinear with GC, reference distance and extension class). Minimum viable evaluation: the hand rule, the curated-record holdout (n about 55) and about 100 lineage-balanced annotated or RNA-supported loci.

Cheaper and more valuable, ranked: (i) fix the `cds_blocks` handling and re-measure; (ii) add **verified C-terminal-complete proteins from well-supported genera** (*Pyrenophora*, *Alternaria*, *Ascochyta*/*Didymella*, where miniprot independently ends in a stop in 205 loci) to the roster as new clade references, reducing reference divergence at the source at no run-time cost (curator promotion only, never automatic, and not used to evaluate themselves); (iii) a per-order profile-HMM completeness instrument; (iv) splice-site PWMs after (i)-(iii).

### 7b. Planned: an intron and splice-site census of `MAT1-1-1` / `MAT1-2-1` (requested by the curator, not started)

Question: how many `MAT1-1-1` genes have introns, what is their splice-site contingent (GT-AG, GC-AG, AT-AC, other), and **does this explain where typical gene predictors fail on MAT genes**? Design:
1. **Intron counts and classes per gene and per clade** from (a) the deposited CDS of the 26 + 29 curated records (independent of the pipeline), (b) the exonerate canonical models in the 3,181 reports (already partly counted: 93.4% GT-AG, 3.3% GC-AG, 1.3% TA-AG over 2,785 models; *Pyrenophora* 40.8% GC-AG; *Aureobasidium* 16% TA donors, which do not occur in nature and so mark mis-placed introns, **[V]**), (c) later, RNA-supported annotations.
2. **A reference distribution:** the same classes over all introns of annotated genomes of the same clades, to test whether MAT genes are enriched in non-canonical or short introns (the `MAT1-1-1` intron of *P. tritici-repentis* is 53-56 bp and GC-AG in all 76 models).
3. **Predictor concordance by splice class:** for each intron class, the fraction recovered exactly by exonerate, by miniprot (default and `--spsc`), and by a typical ab initio/annotation predictor (AUGUSTUS and Helixer where run, and the genome's own annotation where it exists), measured against the curated CDS. The result answers whether non-canonical introns and short introns are *why* MAT genes are hard.
4. Report per lineage cluster (the data are clonal), with n per class; short introns (< 60 bp) and GC-AG counts are small, so state the counts.
This census also supplies the position-weight matrices for stage 3b.

---

## 8. Open questions, not verified

- `bestfit` was not runnable (section 5); no global-alignment mode was evaluated.
- The donor-site proxy (3.5) and the "intron after the alignment edge" failure mode are quantified only as a proxy; no RNA or curated-record check has been run.
- The **380 "clean miniprot alternates longer by a splice jump"** in v1 could not be reproduced by the reviewer (after extension 1,787 of 1,795 clean pairs end at the same stop); that claim is withdrawn pending evidence.
- The 5' residual has no survey beyond the figures in 3.8; `MAT1-2-1` has no chooser or extension run beyond the survey.
- Helixer, AUGUSTUS-PPX, GeMoMa and RNA-seq were not run; AUGUSTUS's installed species list was not checked.
- The reviewer's measurements are the reviewer's; the author has not re-run them. The reviewer ran as a model-agent; the failure message of its first session names `claude-sonnet-5-5`, so **which model produced the review cannot be confirmed from the harness output alone**, and independent re-checks of the headline numbers by a person are recommended before they are quoted.
- Run-to-run noise (one *Ogataea* coordinate, one intermittent test) could mask small effects in a replay.

## 9. Plan and acceptance criteria (pre-committed)

*Superseded in part by section 14: the "wrong-promotion rate 0" and false-complete criteria below were withdrawn after the held-out test; stage 0b/0c are done (sections 11, 13) and the shadow-stage acceptance now reads as in section 14 item 5.*

| stage | work | done when |
|---|---|---|
| 0b | **rescue test (3.4b):** withhold each curated record, align only a truncated reference (30/50/70% N-terminal), extend, compare with the deposited CDS end: rate of exact stop recovery, mean overshoot and undershoot (aa), false-extension rate on complete genes, per clade and per evidence class | the table exists before any threshold on the 30-199 aa tier is chosen; the provisional tier is justified or dropped by this table |
| 0 | holdout baseline on the curated records: exact-match and end-offset rates of today's canonical models, per clade, lineage-weighted | a table, split fixed in advance |
| 1 | `model_selection.py` (3' extension, guards, flags, `cds_blocks`, `phase`), shadow mode, `regression.py` whitelist, tests | `legacy` byte-identical; shadow extension lengths reproduce the survey (4,270 of 4,270) on non-frameshift models and the `cds_blocks` rebuild on frameshift models |
| 2 | shadow on the 3,181 set and the curated-record holdout | **pre-committed acceptance (approved by the curator 2026-10-21):** no promotion that disagrees with a curated record (wrong-promotion rate 0 on the holdout, reported with n); the 0b rescue table shows the provisional tier's overshoot rate; every extended model inside the 1.25 x curated length guard; no complete (3'-closed) model extended; per-lineage table; the three separation <= 0 loci explained |
| 3a | consumer switches (section 4 order), one named bundle, each replayed | each replay's stated criterion; no change to complete-gene calls |
| 3b | exonerate-vs-miniprot swap rule and splice-site check, as a **separate** policy | its own replay and curated-record result |
| 4 | residual analysis (5', interrupted); add a candidate generator only if the residual justifies it | a written decision |
| 5 | roster growth, per-order completeness profiles | held-out improvement |

Principles carried over: defaults stay unchanged until a replay is read; thresholds are provisional and recorded with the evidence they were fitted on; the curator merges PRs; nothing is called validated without a held-out test.

## 10. Changes from v1 (summary of the review)

Verdict go-with-changes. Changed: no tool swap and no 5' move in Phase 1; extension unscored with evidence tiers (ORF-only needs >= 200 aa; 30-199 aa need independent support); `cds_blocks` for frameshift models and no frame search; neighbour, length and donor-site guards; audit corrected (#7 removed, coordinate consumers, report/regression/benchmark added, base commit stated) and the partial-core fear shown to be smaller than stated; "interrupted" counts to be redone; truth set moved to the curated records first; learned chooser dropped; GC-AG and TA-donor splice findings added; the "one best reference", "380 splice jumps" and "miniprot skips an intron" statements corrected or withdrawn.
Retained: the diagnosis, the deterministic pure-function module, shadow mode with a byte-identical `legacy`, clean and 3'-closed as hard filters, the *Botrytis* remnant as 5'-truncated, per-consumer policy flags, no in-house gene finder.

---

## 11. Stage 0b results: rescue test on the curated records (2026-10-21)

`results/2026-10-21_gene_model_extension/stage0b_rescue/`. **Method.** For each of the 34 curated `MAT1-1-1`/`MAT1-2-1` genes with a deposited exon structure (16 + 18; 11 orders), the record (radius `record`) or its whole genus (radius `genus`) was withheld,
every remaining reference of that gene was **truncated to its first 100/70/50/30%** (a simulated fast-evolving tail: the alignment cannot reach beyond that), and the **real** exonerate (`protein2genome --refine region`, `--maxintron 5000`, all references in one query, highest-scoring model kept, as `search.py` does) was run on the record's own ±3 kb genomic window (fetched from NCBI).
The in-frame extension of section 3.3 was applied and its end compared with the deposited CDS end. 272 runs; a model was found in 242 (30 `none`, mostly at 30% truncation); no strand mismatches.
**191 rows had an alignment ending >= 30 aa before the true end** (from 32 distinct genes: rows repeat across truncations and radii, so per-gene numbers are given too); an extension reached a stop in 189 of them.

**Headline [M].**
- **No overshoot.** Over all truncated rows, 0 extensions ended more than 20 aa *beyond* the deposited end (and none beyond 5 aa in the accepted set below). The feared failure (an extension running on past the true stop) did not occur on these 34 genes: the **largest overshoot in any row is +2 aa** (12 rows from two *Teloschistales* genes, one of them probably a deposit-end artefact); the 7 accepted rows that are neither exact, within 5 aa nor > 20 aa short are one gene (*Sclerotinia* `5180_44ba1`) 13 aa short.
- **Exact end recovered in 120 of 191 (63%); within 5 aa in 134 (70%); within 20 aa in 141 (74%).** All 48 remaining failures are **undershoots** (extension stopped > 20 aa short of the true end).
- **By extension length:** >= 200 aa: 41 of 49 exact, 45 within 5 aa, 2 undershoots; 30-199 aa: 79 of 131 exact, 37 undershoots; **< 30 aa: 0 of 9 exact (all undershoot, median 84 aa short)**, consistent with not treating a < 30 aa continuation as an extension when the alignment ended far from the end.
- **The `GT`-donor flag is almost useless here** (4 of 187 truncated extensions end at a `GT`): the alignment edge is rarely at the donor itself, so the flag of section 3.5 does not detect an intron after the edge; do not rely on it for that.

**Why it fails: introns after the edge.** The extension cannot splice, so it stops at the first in-frame stop, which falls inside the next intron when one follows the edge.
| exons in the deposited gene | genes | rows | exact or within 5 aa | undershoot |
|---|---|---|---|---|
| 2 (one intron) | 18 | 108 | 98 (91%) | 10 |
| 3 | 12 | 71 | 36 (51%) | 28 |
| 4 | 2 | 10 | 0 | 10 |
Per gene: 2-exon genes are right in every truncation for 13 of 18 (mean per-gene fraction good 0.85); 3-exon 5 of 12 (0.56); 4-exon 0 of 2. By order, *Pezizales* 0 of 10 exact (both genes have 3-4 exons), *Helotiales* 10 of 33, *Mycosphaerellales* 14 of 24 versus *Pleosporales* 36 of 42 and *Hypocreales* 14 of 16 (small n; not a clean clade effect, confounded with exon count).

**A length-range screen helps but does not fix this.** Rule tested: accept when the ORF extension is >= 200 aa, **or** 30-199 aa **and** the total protein length lies within the curated 10th-90th percentile for that gene (`MAT1-1-1` 234-482 aa, `MAT1-2-1` 176-395 aa); everything else is `provisional`.
Accepted: 145 rows, **123 exact or within 5 aa (85%), 15 undershoots, 0 overshoots**; provisional/rejected: 44 rows, of which 33 are undershoots (75%) and 11 correct. So the screen catches 33 of 48 undershoots at the cost of 11 correct extensions held back as provisional.
**By exon count, the accepted set is 95 of 95 exact or within 5 aa for 2-exon genes** (0 undershoot, 0 overshoot); for 3-exon genes 28 of 48 exact with 13 undershoots and 7 "other"; for 4-exon 0 of 2. The 15 accepted failures are all 3-4-exon *Helotiales* or *Mycosphaerellales* genes whose short total length still falls inside the global range.

**Near-complete alignments (100% references, gap < 30 aa):** 37 rows; the extension is exact in 26 (12 of 16 genus radius, 14 of 21 record radius) and > 5 aa off in 6, so a short slop extension is right about 70% of the time.

**What this changes in the design.**
1. **The 3' extension is reliable for genes with at most one intron after the edge, and unreliable for genes with two or more.** That matches the original concern and gives it a number: for multi-intron genes a splice-aware method is needed, which is the case for a gene-finder-class candidate generator (section 5), now with a measured trigger (undershoot concentrated in >= 3-exon genes).
2. **The reference-structure guard was tested (section 11b): a partial fix.** Originally proposed as follows: use the **reference gene's exon structure** (the curated records carry exons): map the model edge to the reference residue and flag `intron_expected_after_edge` when the reference has an intron downstream; accept an unscored extension only when none is expected, otherwise leave it `provisional`. Test it on this same table before any code.
3. The provisional tier (30-199 aa) is **justified for 2-exon genes and not for multi-intron genes**; the total-length screen should use **per-order** ranges (global ranges let 15 undershoots through) once each order has >= 3 curated records.
4. **No overshoot** supports reporting extensions at all (with their class), but the acceptance bar "zero wrong promotions" must be read as "zero overshoot": undershoot is a **missed end** (a partial gene reported as extended), which the curator should decide how to count.

**Limits.** 34 genes in 11 orders, skewed toward well-studied taxa; rows are not independent (32 genes, 4 truncations x 2 radii); truncation is simulated by cutting the reference, so it shows whether the extension finds the right stop, **not** whether a real divergent tail would still align elsewhere in the protein; the reference used for each model was not recorded; genes with frameshift models were not special-cased; the window flanks are 3 kb (an extension never reached the window end in 189 of 191). The 85/80 and the length-range rules were chosen after seeing these genes, so the 123/145 figure is fitted, not held out.

### 11b. The reference-structure guard (tested 2026-10-21, job 30234694)

`rescue_test2.py` repeats the test and also records, per hit, the best-hit reference and where the alignment ended on it (`vulgar`), and the number of introns of the **deposited gene** lying after the model end (the oracle cause). The 242 models are identical to the first run. Reference intron positions (aa) come from the curated records' exons (`ref_introns.tsv`); only **34 of the 55 curated references carry an exon structure**. 189 truncated rows reached a stop.

**Oracle (confirms the cause) [M].** Introns after the edge in the deposited gene: **0 -> 128 of 132 rows exact or within 5 aa (97%)**; 1 -> 6 of 41; 2 -> 0 of 9; 3 -> 0 of 7. Undershoot (> 20 aa short) occurs in 28 of 41 (1 intron), 9 of 9, 7 of 7 and 4 of 132 (none). An intron after the alignment edge is the cause of the failures.

**Guard = the best-hit reference has >= 1 intron downstream of the alignment end** (144 rows whose reference has exon info):
| | undershoot | fine |
|---|---|---|
| flagged | 17 | 4 |
| not flagged | 15 | 108 |
Precision 0.81, **recall 0.53**; the not-flagged rows are 108 of 123 right (88%); 4 flagged-but-correct extensions would be held back (4 of 21 flagged).
Why recall is limited: the reference and the target often differ in intron number or position (across the 144 rows the reference under-predicts the number of introns after the edge in 19 rows (13 rows 0 vs 1, 4 rows 1 vs 3, 2 rows 2 vs 3) and over-predicts in only 2 rows (2 vs 1); 123 rows agree), so reference structure is a partial proxy. 45 of the 189 rows (24%) have a best reference **without** exon info and fail more often (16 undershoots, 22 correct).

**Combined with the length-range rule of section 11:**
| accepted set | rows | exact or within 5 aa | undershoot | overshoot > 5 aa |
|---|---|---|---|---|
| length rule alone | 145 | 123 (85%) | 15 | 0 |
| length rule **and** reference known with no intron after the edge | 111 | 101 (91%) | 10 | 0 |
| length rule but guard not clear (intron after, or reference without exon info) | 34 | 22 | 5 | 0 |

**Conclusions.**
1. The guard raises acceptance precision from 85% to 91% at the cost of holding back 22 correct extensions as provisional; it is **worth having as a flag, not as a sufficient gate**: about 10% of accepted extensions in multi-intron genes would still stop early.
2. **Record the exon structure of every reference in the roster** (cheap, no run-time cost): 21 of 55 references lack it today; with full coverage the 24% "unknown" rows become judgeable. For references from deposits without an annotated CDS, derive the structure from the pipeline's own alignment to the source locus and have a curator confirm it.
3. The remaining early stops (an intron the reference does not predict) cannot be found from the reference alone. They need a sequence-level signal: a **splice-aware extension** (search for a donor/acceptor pair that reopens the frame), which is the same job a gene-finder-class candidate generator does, so this result is further support for adding one in stage 4 and for the miniprot `--spsc` route of section 5b. A first deterministic version to test on this same table: from the alignment edge, enumerate GT/GC donors and AG acceptors (40-300 bp apart), keep spliced continuations that stay open for >= 100 aa to a stop and land inside the curated length range, and report the **best-scoring** one with the number of competing candidates (ambiguity) rather than a single answer.
4. Same limits as section 11 (34 genes, simulated truncation, rules fitted on these genes), plus: reference intron positions are in aa along the reference protein and the alignment end is taken from `vulgar`; the position of a phase-1 or phase-2 intron is rounded to 0.1 aa.

### 11c. A minimal splice-aware extension, tested (`splice_ext.py`, curator asked for a test, not a build)

From each alignment edge, enumerate `GT..AG` introns of 40-400 bp that keep the frame open to a stop at least 30 aa beyond the unspliced stop; choose, among candidates whose total protein length is in the curated 10th-90th percentile, the one closest to the curated median length. Deterministic; sequence and the curated length range only. 189 truncated rows.

- **Highly ambiguous.** Median **18 candidate introns per row** (90th percentile 69); **12** inside the length range (90th percentile 48); the **true end is reachable by a candidate in only 28 of 189 rows** (26 inside the length range): the true intron is often not a `GT..AG` of 40-400 bp the search can reach, or the donor is not at a position the open frame allows.
- **Applied everywhere it is harmful:** it would replace the correct unspliced answer with a wrong spliced one in 103 of the 134 rows the unspliced extension got right (exact or within 5 aa overall drops from 134 to 51 of 189).
- **As a fallback only (applied to the 44 extensions the length rule declined):** 36 have a candidate; the spliced choice is exact in 13 (36%) against 7 within 5 aa for the unspliced extension on the same rows; for the 29 genuine early stops with a candidate it is exact in 13 (45%). The in-range candidates for those rows number a median of 48, and the 3 rows with exactly one candidate were all wrong. Restricted to the rows the reference-intron guard flags: 12 exact of 15 (unspliced: 0 of 15), **but those 15 rows come from only 6 genes (4 with any exact row; one gene supplies 6 of the 12)**, so this is not 15 independent successes.
- **Reading.** A hand-written donor/acceptor enumeration is ambiguity-limited: it works on a minority of genes (one-intron genes whose intron the length prior can pick) and is wrong otherwise. This is a **negative result for building it ourselves**: if splice-aware extension is wanted, the right tool is one with trained splice and exon-length models (the miniprot `--spsc` route of section 5b with clade-trained splice scores, or a gene-finder candidate generator), not this enumeration. **Recommendation: do not implement it**; keep multi-intron cases `provisional` and decide on a gene-finder candidate in stage 4 on the truth-set residual.

## 12. Decision needed: how acceptance counts an early stop

*Historical: decided 2026-10-21 (option 2) and withdrawn 2026-10-22 (section 14). Kept because the reasoning and the numbers are the record of why.*

An extension that stops early is **never worse than not extending**: in 187 of the 189 truncated rows the extended model ends closer to the true end (2 unchanged, 0 farther); of the 48 early stops, 46 are closer (median distance to the true end 162 aa against 262 before). What it gets wrong is the **label**: a still-partial gene reported as extended or complete. Today's pipeline treats all 189 as complete (the ORF test passes vacuously), so any acceptance rule improves on the legacy.
**Curator decision 2026-10-21: option 2 (below).** Options for the acceptance bar (section 9):
1. **Strict:** any end error beyond 5 aa in either direction is a wrong promotion; the bar is zero. Not attainable now (best accepted set: 10 of 111 early stops); it would mean never promoting multi-intron genes.
2. **Overshoot zero, early stops capped and reported separately (recommended):** zero accepted extensions more than 5 aa beyond the deposited end (measured: 0, maximum +2 aa); early stops (more than 20 aa short) are counted as a **false-complete rate**, reported per exon-count class and per order, with a ceiling the curator sets (measured on these genes: 0 of 95 for one-intron genes, 10 of 111 overall with the reference guard). Gene models never labelled `complete` unless the guard is clear.
3. **Relative bar:** no accepted extension may end farther from the true end than the unextended model (measured: 0 of 189), and the false-complete rate must be lower than the legacy rate in every order (legacy: 100% of truncated alignments are called complete).

**Decided (curator, 2026-10-21): option 2.** The acceptance table of section 9 reads: zero accepted extensions more than 5 aa beyond the deposited end; early stops (more than 20 aa short) are reported as a false-complete rate per exon-count class and per order against a ceiling the curator sets; a gene is labelled `complete` only when the reference-structure guard is clear.

### 11d. Exon structures for the whole roster (curator approved 2026-10-21; a proposal table, db untouched)

21 of the 55 reference proteins had no exon structure (22 gene rows in 18 records; one of the 22 rows, *A. fumigatus* A1163 `MAT1-2-1`, has no protein in the reference set, which is why it is 21 and not 22). `results/2026-10-21_gene_model_extension/roster_exon_structures/roster_exon_structures_proposed.tsv` derives them from the NCBI protein records: **15 ok** (CDS length matches the database protein, accession matches), **5 partial_cds** (Pezizales deposits incomplete at an end; the CDS length still matches the database protein), **2 without a database protein to cross-check**. Exon counts: 1 exon x 1, 2 x 8, 3 x 10, 4 x 3. Four RefSeq proteins needed the genomic CDS feature because `coded_by` names an mRNA. **A curator must confirm before the metadata is edited.**
Effect on the guard of section 11b, re-run with all references (`analyze4.py`): the guard can now judge **189 of 189 truncated rows (was 144 of 189)**.
| | before (144 judged) | after (189 judged) |
|---|---|---|
| flagged (reference has an intron after the edge) | 21 | 42 |
| flagged and an early stop | 17 | 31 |
| precision / recall for early stops | 0.81 / 0.53 | **0.74 / 0.65** |
| not flagged, correct | 108 of 123 (88%) | 130 of 147 (88%) |
| length rule **and** guard clear: rows / correct / early stops / overshoot | 111 / 101 (91%) / 10 / 0 | **129 / 119 (92%) / 10 / 0** |
The full roster covers every row, raises recall, and the combined acceptance rule now accepts 129 extensions at 92% correct with no overshoot, holding back 16 (4 of them correct). The same limits apply (34 genes, simulated truncation, rules fitted on these genes).

### 12b. Options for the false-complete ceiling (measured on the 189 truncated rows with full roster structures; Wilson 95% intervals)

| rule that labels an extension `complete` | accepted | early stops | rate (95% CI) | genes with an early stop |
|---|---|---|---|---|
| E. every extension of >= 30 aa, no screen | 180 | 39 | 21.7% (16-28%) | 12 |
| A. length-range rule only | 145 | 15 | 10.3% (6-16%) | 4 |
| B. length rule + reference-intron guard clear | 129 | 10 | 7.8% (4-14%) | 4 |
| C. B and ORF extension >= 200 aa only | 42 | 1 | 2.4% (0-12%) | 1 |
| **D'. B and the best-hit reference has <= 2 exons (usable at run time)** | **93 (15 genes)** | **0** | **0.0% (0-4%, rows as independent)** | **0** |
| F. B and the best-hit reference has >= 3 exons | 26 | 10 | 38.5% (22-57%) | 4 |

(An earlier variant restricted to genes with <= 2 exons used the deposited gene's own exon count, which is not known at run time: that was an oracle and is not an option. The run-time proxy is the best-hit reference's exon count; reference 2 exons vs deposited 2 / 3 / 4 exons: 95 / 7 / 0; reference 3 exons: 10 / 49 / 8; reference 4 exons: 0 / 4 / 2.)
All 10 early stops among B's accepted set come from 4 genes (*Helotiales* 3, *Mycosphaerellales* 1), all with a 3-4-exon reference; by deposited exon count, 2-exon genes 0 of 95, 3-exon 8 of 32 (25%), 4-exon 2 of 2. By order: *Helotiales* 8 of 14, *Mycosphaerellales* 2 of 14, every other order 0 (*Pleosporales* 0 of 36, *Teloschistales* 0 of 24, *Hypocreales* 0 of 14).

Cutoff options for the curator (section 12 option 2):
1. **One global ceiling** on rule B's false-complete rate, for example <= 10% (measured 7.8%, CI 4-14%, so a 10% ceiling is not statistically distinguishable from the measurement; with n of 129 from 23 genes in 9 orders it cannot be tested more finely).
2. **Per-class ceilings, with the reference exon count as the class (recommended):** `complete` only when rule D' holds (reference has <= 2 exons); ceiling **5%** for that class (measured 0 of 93, upper CI 4%; the 93 rows come from only **15 genes in 7 orders**, so the interval treats rows as independent and is optimistic); every extension whose best-hit reference has >= 3 exons is `provisional` and never `complete` (measured 38% would be false-complete); the provisional class carries no ceiling but is reported.
3. **Zero tolerance:** a ceiling of 0% for `complete` with rule D' (measured 0 of 93, but a zero bound cannot be verified on 93 rows: it is the observed rate, not a guarantee; the CI upper bound is 4%).
4. **Relative bar:** the false-complete rate under the new rule must be lower than the legacy rate (100% of truncated alignments are called complete) in every order, with no ceiling number.
Caveats: n is small and clustered (one order and a few genes drive every early stop); the thresholds and rules were chosen after seeing these genes, so the figures are fitted, not held out; the real test of any ceiling is the same measurement on genes not used here (new curated records, or the held-out radius of the curated-record holdout).

**Decided (curator, 2026-10-21): per-class ceilings (option 2 of 12b).** `complete` is allowed only when the best-hit reference has <= 2 exons **and** the length rule passes **and** the reference-intron guard is clear; the false-complete ceiling for that class is **5%**, always reported with n, the number of distinct genes and the orders covered (measured so far: 0 of 93 rows from 15 genes in 7 orders). Every extension whose best-hit reference has >= 3 exons is `provisional` and is never labelled `complete`; that class has no ceiling and is reported. The ceiling is verified on genes not used to choose the rule (new curated records or the held-out radius of the curated-record holdout) before it is called met.
**Curation of the roster exon structures (curator, 2026-10-21):** a curator may be an Opus agent or the project curator depending on complexity. The 15 `ok` structures are confirmed by an independent Opus agent re-deriving them from NCBI; the 5 partial deposits and the 2 rows without a database protein are decided by the project curator on the agent's evidence. A separate draft PR adds only confirmed structures to the records.

**Roster PR (2026-10-22):** the curator-confirmed exon structures for the 22 references were merged as #119 (`roster-exon-structures`; the `record_version` bumps it missed are in #123): 17 verified rows plus the five *Tuber* rows (canonical GT-AG boundaries for the two `MAT1-2-1`, the deposit's for the three `MAT1-1-1`, per the curator). The guard of 11b/11d needs this PR merged to have full coverage at run time.

**Merge order (curator, 2026-10-22):** roster PR #119 merged before this design PR (#118; done): the reference-structure guard of 11b/11d needs the full-coverage roster. A Fable agent reviews #119 before merge.

### 13. Stage 0c: held-out validation of the frozen rules (pre-registered 2026-10-22, before any result)

**Held-out set.** Direct GenBank deposits of `MAT1-1-1` / `MAT1-2-1` from NCBI protein records (esearch of the gene and title fields over Fungi: 3,763 records), filtered to: a named `MAT1-1-1` or `MAT1-2-1`; not a RefSeq protein; source nucleotide a small direct deposit (old-style accession, 1-2 letters + 5-6 digits, not a genome contig); CDS with >= 2 exons and complete (no `<`/`>`); **genus not in the 55 curated reference set** (so no reference is a close relative by genus); >= 1,000 bp of sequence downstream of the CDS; at most 2 genes per genus, spread over exon counts, PubMed-cited records preferred. Result: **106 genes, 66 genera, 20 orders; 52 two-exon, 44 three-exon, 10 four-exon; 57 `MAT1-1-1`, 49 `MAT1-2-1`; 64 cite a PubMed article.** Selection script and tables are in `results/2026-10-21_gene_model_extension/heldout/`.
**Caveat on truth quality:** these are the depositors' annotated CDS structures, many from population or phylogeny studies of PCR-amplified MAT regions; they are not RNA-supported, and some may be predictions. They are independent of this pipeline and of the 55 references, which is what a held-out test needs, but a disagreement is ambiguous (the deposit may be wrong).
**Frozen rules (not re-tuned on these genes).** Method exactly as stage 0b: the real exonerate with production settings, all 55 curated references truncated to 30/50/70/100%, the in-frame extension, comparison with the deposited end; length rule (>= 200 aa ORF extension, or 30-199 aa with total length inside the curated 10th-90th percentile of 234-482 aa for `MAT1-1-1` and 176-395 aa for `MAT1-2-1`), reference-intron guard (best-hit reference has no intron after the alignment edge, using the full 56-row roster table), and class D' (best-hit reference has <= 2 exons).
**Pre-registered criteria.** (1) Overshoot: zero accepted extensions more than 5 aa beyond the deposited end. (2) The `complete` class D': false-complete rate (> 20 aa short) <= 5%, reported with n, the number of distinct genes and orders, and an interval computed per gene rather than per row. (3) Reported for transparency, not as pass/fail: the rate for classes with a >= 3-exon reference, by exon count and by order, and how the stage 0b rates (7.8% for rule B, 0 of 93 for D') compare. **If criterion 2 fails, the ceiling or the rule is revised in a new, dated section and the failure stays in this note.**

### 13b. Stage 0c results (2026-10-22): the frozen rules did NOT meet the pre-registered criteria

HPCC job 30258079; 106 held-out genes, 100 with a model, 424 runs (358 with a model; 66 `none`), **334 truncated alignments reaching a stop from 100 genes**; no strand mismatches. Rules applied exactly as frozen in section 13. **Tables: `results/2026-10-21_gene_model_extension/heldout/`.**

| rule | accepted rows | early stops (rows, 95% CI) | genes with an early stop (per gene, CI) | overshoot > 5 aa | orders |
|---|---|---|---|---|---|
| E. every extension >= 30 aa | 316 | 49 (15.5%, 12-20%) | 25 of 99 (25%, 18-35%) | 8 | 18 |
| A. length rule only | 272 | 33 (12.1%, 9-17%) | 15 of 84 (18%, 11-27%) | 8 | 17 |
| B. + reference-intron guard clear | 232 | 17 (7.3%, 5-11%) | 5 of 77 (6.5%, 3-14%) | **6** | 17 |
| D'. B + best-hit reference <= 2 exons | 172 | **14 (8.1%, 5-13%)** | 5 of 54 (9.3%, 4-20%) | **6** | 15 |
| F. B + best-hit reference >= 3 exons | 60 | 3 (5.0%, 2-14%) | 1 of 27 (3.7%, 1-18%) | 0 | 13 |

**Criteria (section 13).**
1. **Overshoot: FAIL.** 6 accepted rows overshoot by more than 5 aa (maximum +95 aa), from 2 genes: *Calonectria naviculata* `MAT1-1-1` (+95 aa, deposit 372 aa, extension total 416) and *Diplocarpon coronariae* `MAT1-1-1` (+71 aa, deposit 306 aa, total 395 or 352). Stage 0b had none.
2. **`complete` class D' <= 5% false-complete: FAIL.** 14 of 172 rows (8.1%), 5 of 54 genes (9.3%, CI 4-20%). The rows come from five genes: *Beauveria bassiana*, *Ophiocordyceps sinensis*, *Ustilaginoidea virens* (all *Hypocreales*, 12 of the 14 rows; deposited gene with **three** exons, best-hit reference *Fusarium* with **two**), *Ascochyta lentis* and a *Pseudogymnoascus* (1 row each).
3. **Reported:** the exon-count distinction **did not generalize**: class F (best-hit reference >= 3 exons) had 5.0% here, against 38.5% in stage 0b; class D' had 8.1% here, against 0% (0 of 93). The ordering the curator's per-class ceiling was built on was an artefact of 34 genes. The stage 0b rate for rule B (7.8%) is close to this one (7.3%), so rule B's overall false-early-stop rate is stable at about 7-8%; it is the *class split* that failed.

**The mechanism replicates, the proxy does not.**
- **Oracle (true introns in the deposited gene after the alignment edge):** 0 -> **246 of 252 rows exact or within 5 aa (97.6%)**, 6 early stops; >= 1 -> **2 of 82 (2.4%)** (61 early stops, 8 overshoots, 11 neither; 75 rows with one true intron after the edge and 7 with two). Stage 0b: 97% and 11%. **All 8 overshoots occur where one true intron follows the edge:** the unspliced extension reads through the intron to a later stop. Among rule B's accepted rows: 198 of 199 right when no true intron follows the edge, **0 of 33 right when one does.**
- **The reference as a proxy fails in this held-out set:** the best-hit reference's exon count equals the deposited gene's in 248 of ~330 rows; it differs in **77 rows**, of which 37 are reference 2 exons vs deposited 3 and 3 are 2 vs 4 (the reference lacks an intron the target has, mostly in *Hypocreales*), and 23 are reference 3 vs deposited 2. The reference-intron guard flags 48 rows (26 early stops, 2 overshoots, 20 right) but misses 34 early stops and 6 overshoots among the unflagged.
- **So the extension itself is sound whenever no intron follows the alignment edge (97.6% right), and unsafe otherwise; what has not been solved is predicting from the reference or the sequence whether an intron follows the edge.**

**Consequences.**
1. The per-class ceiling decision of section 12b (`complete` only for reference <= 2 exons, 5% ceiling) **cannot be met with the rules as frozen** (8.1%, and overshoot occurs). It must be revisited by the curator; this note makes no change to the recorded decision until then.
2. The reference's exon count is not a reliable class label across clades (*Hypocreales* `MAT1-1-1` carry an extra intron the *Fusarium* references lack); a rule tuned on 34 genes did not transfer to 100. This is the argument for the pre-registered held-out step, and it worked as intended.
3. Options (not chosen): (a) report every extension with its evidence and **never label it `complete`** (rescue use only; the false-early-stop rate of about 7-8% for rule B is disclosed per clade); (b) develop a better "intron after the edge" predictor (for example, use the intron positions of **all** references aligned to the region rather than the best hit only, per-clade reference structure, or the splice-site route of 5b), and validate it on a **third, untouched set**: 101 further direct-deposit genes with >= 1 kb downstream flank were not used here and are reserved (`heldout/candidates_direct.tsv` minus `selected.tsv`); (c) accept a higher ceiling (the curator's call; the measured rates are 7-8% overall, clade-dependent).
4. The 106 genes are now **development data**: any new rule examined on them needs the reserved set for its own validation.
Caveats: truth is the depositors' annotation (not RNA-supported; a deposit may itself be wrong, e.g. the two overshoot genes could be deposits shorter than the real gene, which this test cannot distinguish); rows repeat across four truncations (genes are the unit that counts, so per-gene figures are given); the selection prefers PubMed-cited records and at most two genes per genus.

## 14. Decision (curator, 2026-10-22): option 1, rescue only; supersedes the per-class ceiling of section 12b

After the held-out validation of section 13b the curator chose **option 1**: the extension is a **rescue** step and never earns the label `complete`.
**What this means operationally (to implement; nothing is implemented yet):**
1. The 3' extension is computed whenever its preconditions hold (>= 30 aa open frame to a stop; no contig edge; no neighbouring modelled CDS overlap; total length <= 1.25 x the longest curated record) and is **always reported with its evidence**: `ext_aa`, `ext_evidence` (`orf_only`, `orf+independent_stop`, `orf+tail_similarity`, `orf_weak`), whether the length rule passes, whether the reference-intron guard is clear, the best-hit reference's exon count, `donor_GT`, and `contig_edge`/`neighbour` flags. The model class is `extended_provisional`.
2. **`complete` is reserved for models whose own alignment already reaches the stop** (`closed_by_alignment`): the extension never promotes a model to `complete`. The 5% false-complete ceiling and the per-class `complete` rule (12b, 12 "Decided") are **withdrawn**, because the held-out rates (7.3% for rule B, 8.1% for class D', overshoot up to +95 aa) do not meet them.
3. **Disclosed error rates** accompany the feature in the documentation and the report schema, per clade where n allows: the false-early-stop and overshoot rates of section 11 (34 fitting genes: rule B 7.8%, 0 overshoot) and 13b (106 held-out genes: rule B 7.3% early stops, 6 accepted overshoot rows from 2 genes, class D' 8.1%); the oracle (97.6% right when no intron follows the alignment edge, 2.4% otherwise); and the statement that these are depositors' annotations, not RNA-supported truth.
4. **Downstream consumers** (section 4) may read `extended_provisional` as a lower-confidence input behind their own policy flags, each replayed as already planned; none treats it as `complete`. The homothallic completeness test (#117) may use an extended length at reduced confidence; it must not report `complete` from it.
5. **Shadow mode first** (section 3.9): `legacy` stays byte-identical; the acceptance table for the shadow stage becomes: no change to any legacy call; extension lengths reproduce the survey; the disclosed rates reproduce on a fresh re-run; no extension labelled `complete`. The false-complete ceiling no longer applies because nothing extended is labelled `complete`.
6. A better "intron after the edge" predictor remains an open research item. A feasibility check (2026-10-22, 130 genes, gene-held-out cross-validation) found cheap signals insufficient: reference-intron flag AUC 0.67; with length and splice-signal features 0.79; self-trained codon-usage coding potential 0.48 (no value); the lowest-risk 15% of rows by model risk was still 7.8% wrong (24.8% overall). The bounded experiment now under way tests whether a larger, leave-genus-out roster improves the reference guard (section 14b, when finished). Any predictor is validated on the 101 reserved genes (section 13b).

### 14b. Roster-growth experiment (2026-10-22): more references with exon structures improve the guard, modestly

Question: is the guard weak because the roster is thin? Method: the 55 curated references plus the 106 held-out deposited proteins (161 references, with their exon structures) against the same 140 genes (34 curated + 106 deposits), **leave-genus-out** (every reference of a gene's own genus withheld), the real exonerate with production settings, references truncated to 30/50/70/100%, rules frozen as before. Baseline: the same genes and truncations with the 55 references (genus radius for the 34; the 106 are already genus-disjoint). HPCC array 30261437. **386 truncated rows from 127 genes compared**; files in `results/2026-10-21_gene_model_extension/roster_growth/`.

| | 55 references | 161 references |
|---|---|---|
| best-hit reference has exon info | 372 of 386 | 382 of 386 |
| best-hit reference exon count equals the deposited gene's | 292 of 372 (78%) | 313 of 382 (82%) |
| guard ("reference has an intron after the edge") precision / recall for a true intron after the edge | 0.79 / 0.52 | 0.77 / **0.68** |
| AUC of the guard | 0.734 | **0.799** |
| rule A (length rule only): wrong (early or overshoot > 5 aa) | 14.6% (CI 11-19%) | 15.5% (12-20%) |
| **rule B (length rule + guard clear): accepted / early stops / overshoot / wrong** | 261 / 16 / 6 / **8.4%** (6-12%) | 246 / 8 / 6 / **5.7%** (3-9%) |
| rule B genes with any wrong row (per gene) | 7 of 95 (7.4%) | 5 of 101 (5.0%) |
| class D' (B and reference <= 2 exons) wrong | 9.6% (6-15%) | 5.9% (3-10%) |
| class F (B and reference >= 3 exons) wrong | 4.7% (2-13%) | 5.3% (2-13%) |
The best-hit reference came from the added deposits in 273 of 386 rows.

**Reading.** More references with recorded structures help: guard recall rises from 0.52 to 0.68, early stops among rule B's accepted rows halve (16 to 8), and the wrong rate falls from 8.4% to 5.7% (per-gene 7.4% to 5.0%). **Overshoot does not improve (6 rows, the same two genes)**, and precision is unchanged, so roster growth reduces missed introns but does not remove the through-intron overshoot. The reference's exon-count agreement is still 82%, so about a fifth of best-hit references disagree with the target's structure even with 161 references.
**Caveats.** The 106 genes are development data and the rules were frozen, but this is **not a validation**: a validation uses the 101 reserved genes with the same leave-genus-out roster; leave-genus-out does not remove same-family relatives (realistic for roster growth, but it flatters closely related clades); deposits are the depositors' annotations (some may be predictions); rows repeat across four truncations (per-gene figures given); the intervals treat rows as independent.
**Consequence for the plan.** Under the rescue-only decision (section 14) roster growth is the **cheapest lever that measurably lowers the disclosed rates**: record exon structures for every roster reference (PR #119 and the follow-up agent branch do this for the existing 55) and **add annotated, genus-diverse references** (for example the deposits used here, after curation). The reserved 101 genes can validate that effect; whether to spend them is the curator's call.

## 15. Validation of the roster-growth effect on the reserved set (pre-registered 2026-10-22, before any result)

**Decision (curator, 2026-10-22):** spend the 101 reserved genes now, before investing curation effort in a larger roster.
**Set.** The direct-deposit candidates of section 13 with >= 1,000 bp downstream flank that were **not** among the 106 selected genes (207 - 106 = 101 genes), fetched with +/-3 kb windows; same quality caveat as section 13 (depositors' annotations, not RNA-supported).
**Configurations (identical genes, truncations, rules and settings):** (a) baseline: the 55 curated references; (b) augmented: the 161 references of section 14b (55 curated + 106 deposits with exon structures). Both **leave-genus-out** (every reference of a gene's own genus withheld, including the deposits). The real exonerate with production settings, references truncated to 30/50/70/100%, rules frozen as in section 13 (length rule, reference-intron guard, class labels).
**Pre-registered criteria** (the effect "replicates" only if both hold):
1. Rule B (length rule + guard clear): the wrong rate (early stop > 20 aa short, or overshoot > 5 aa) is **lower** under the augmented roster than under the baseline on the same genes, **and the per-gene wrong rate is not higher**.
2. Guard recall for a true intron after the edge rises by **at least 0.10** (absolute) with precision not falling by more than 0.10.
Reported regardless: overshoot counts (must not worsen: the 2 held-out overshoot genes showed it does not improve), class D'/F, per-gene figures, intervals per gene.
**If either criterion fails**, the roster-growth lever is recorded as not validated and the curation investment is **not** recommended; the failure stays in this note.

### 15b. Result (2026-10-22): the roster-growth effect replicates on the reserved genes; both pre-registered criteria are met

HPCC array 30263502; 101 reserved genes (98 with a model; **288 truncated rows from 98 genes**, same genes and truncations in both configurations); leave-genus-out; rules frozen. Files: `results/2026-10-21_gene_model_extension/roster_growth_validation/`.

| | 55 references (baseline) | 161 references (augmented) |
|---|---|---|
| best-hit reference has exon info | 278 of 288 | 286 of 288 |
| best-hit reference exon count equals the deposited gene's | 177 of 278 (64%) | 233 of 286 (81%) |
| guard precision / recall for a true intron after the edge | 0.74 / **0.32** | 0.74 / **0.49** |
| guard AUC | 0.638 | 0.697 |
| rule A (length rule only): wrong | 16.2% (12-22%) | 14.3% (10-20%) |
| **rule B (length rule + guard clear): accepted / early / overshoot / wrong** | 170 / 17 / 5 / **12.9%** (9-19%) | 172 / **3** / 5 / **4.7%** (2-9%) |
| rule B genes with any wrong row (per gene) | 7 of 60 (11.7%) | 5 of 81 (6.2%) |
| class D' (B and reference <= 2 exons) wrong | 14.6% (9-22%) | 6.8% (3-13%) |
| class F (B and reference >= 3 exons) wrong | 8.5% (3-20%) | 1.4% (0-8%) |
The best-hit reference came from the added deposits in 239 of 288 rows.

**Criteria (section 15).** (1) Rule B's wrong rate is lower with the augmented roster (4.7% against 12.9%) and the per-gene rate is not higher (6.2% against 11.7%): **met**. (2) Guard recall rises by at least 0.10 (0.32 to 0.49, +0.17) with precision unchanged (0.74 to 0.74): **met**. **The effect replicates on genes never used before.** As in the development set, overshoot does not improve (5 accepted rows in both configurations): roster growth lowers the early stops (17 to 3) but not the reading-through-an-intron overshoot.
**Caveats.** The reserved set is more clonal than the development set (32 genera; *Phialocephala* 17 genes, *Stagonosporopsis* 8, *Didymella* 7, *Calonectria* 7), so the per-gene intervals are wide and the baseline is worse here than on the development genes (12.9% against 8.4%); leave-genus-out keeps same-family relatives in the roster (realistic for roster growth, flattering for closely related clades); the truth is the depositors' annotation, not RNA-supported; the rules were chosen on the first 34 genes and are unchanged.
**Consequence.** Growing the roster with genus-diverse annotated references that carry exon structures is validated as the cheapest lever that lowers the disclosed early-stop rate (about 3-fold fewer early stops on the reserved genes). The untouched set is now spent; any later predictor needs new deposits.

