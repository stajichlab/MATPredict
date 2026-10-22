# Gene-model selection and 3' extension: design, blind-spot audit and evaluation plan (2026-10-21)

Status: **proposal, nothing implemented.** This note is a design for review. It is written so that a reviewer with no
access to the conversation can check it. Statements marked **[M]** were measured, **[R]** were read in code, **[I]** are
inference or recollection not tested. Where I am unsure of a tool's behaviour I say so (section 8).

Origin: the homothallic-loci patch (PR #117, `docs/notes/2026-10-21_homothallic-patch-measurement.md`) found that the
196 aa *Pyrenophora tritici-repentis* `MAT1-1-1` model ends at the edge of the exonerate alignment, not at the gene's end.
A read-only review (Fable agent, 2026-10-21) then surveyed all 3,181 patched-run reports and prototyped a model chooser.
Its report is summarised here; the prototype, survey table and null test are in
`results/2026-10-21_gene_model_extension/`.

---

## 1. The problem, in one table

For each gene, the pipeline models the target with `exonerate --refine region` and `miniprot`. `polish.classify()` makes
**exonerate canonical whenever it produced a model** (`src/MATPredict/detect/polish.py:116-118`, **[R]**), whatever its length, stops or
agreement with miniprot. Exonerate (`protein2genome`, local) ends where its alignment to the *one best reference* ends **[R/M]**.

Survey of 3,181 reports, in-frame continuation after the canonical model's last exon **[M, `survey.tsv`]**:

| 3' state of the exonerate model | MAT1-1-1 (n=1,721) | MAT1-2-1 (n=1,781) |
|---|---|---|
| internal stop in the model | 142 | 350 |
| next codon is a stop (complete) | 182 | 444 |
| stop within 1-29 aa (edge slop only) | 281 | 82 |
| open 30-99 aa | 184 | 47 |
| open >= 100 aa | 924 | 855 |

So 64% of `MAT1-1-1` and 51% of `MAT1-2-1` canonical models have >= 30 aa of open frame past their end. The effect is clade-structured:
essentially all *Macrophomina*, *Botryosphaeria*, *Pyrenophora*, *Corynespora*, *Phyllosticta*, *Didymella*; about 60% of *Aureobasidium* (a further 63 have internal stops);
*Parastagonospora*, *Zymoseptoria* and *Alternaria* are essentially complete. Of `polished_disagree` rows, 76% (`MAT1-1-1`) and 66% (`MAT1-2-1`) have an open exonerate frame:
"disagree" is mostly "exonerate stopped early", not two competing structures.

Consequences already seen: the ORF-completeness test passes vacuously (it counts internal stops, and an edge-ended model has none); a length floor (200, then 190 aa)
tests aligned extent, not the gene; *P. tritici-repentis* passed the 190 aa floor by 6 aa.

Facts corrected during review (these were wrong in earlier notes) **[M]**:
- The *Botrytis* remnant's in-frame stop is at **residue 10**, not 209: the 219 aa model is `LRDVDARRI*SKYFALL...`; the 209 aa protein is the stop-free stretch after it, with no M, and it is 3'-complete. It is **5'-truncated**.
- Exonerate models never contain a stop codon (0 of 3,502 rows); miniprot alternates end with one in about 19%. "Terminal stop in the model" is not a usable test on the canonical; "the next codon is a stop" is.
- `_translate_model` (`pipeline.py:821-880`, **[R]**) returns the longest stop-free stretch over **all three frames**, not the canonical frame, in spite of its docstring ("stopping at the first stop"). A frame error or an intron in the CDS can still produce a long "ORF".

---

## 2. What "right" means: definitions

For one target gene at one locus, a **candidate model** is a spliced CDS (list of exons, strand, frame) with its translation.

| term | definition (deterministic, from sequence only unless stated) |
|---|---|
| `clean` | no internal in-frame stop (using the frameshift-corrected block translation when the model has `cds_blocks`, `pipeline.py:837-850`) |
| `3'-closed` | the codon after the CDS is a stop (or the CDS ends in one) |
| `3'-open` | not 3'-closed, and the frame continues >= `EXT_MIN_AA` before the next stop |
| `starts M` | first codon is a start |
| `n_term_truncated` | no M, and an in-frame stop lies <= `NTERM_STOP_WINDOW` aa upstream with no M between that stop and the model start |
| `complete` | clean, 3'-closed, (starts M or flagged `n_term_unknown`), aligned coverage of the best reference >= 80% |
| `partial` | clean, aligned coverage 50-80% |
| `remnant` | clean, identity >= 85% to the best reference and aligned coverage < 80% (the existing rule, `homothallic.py`), or `n_term_truncated` with coverage < 80% |
| `interrupted` | internal stops remain in every candidate; reported with the first-stop position and how many reference columns align *after* it |

Principle (the curator's rule, refined): **among candidates, prefer longer *clean* models that are 3'-closed and better explained by the reference; a longer model with
internal stops is never promoted; a clean model that is a truncated copy is flagged, not called complete.** Longer-without-stops is a good ranking default but not a correctness test:
miniprot can skip an intron (the 419 aa model with 6 stops), a clean model can be a remnant (*Botrytis*), and some stop-interrupted *Aureobasidium* models have the reference aligning across the stops
(missed intron or frameshift rather than a pseudogene, **[M]** 79-261 reference columns after the first stop, at 31-44% identity).

---

## 3. The algorithm: deterministic code, no model in the loop

Everything below is ordinary Python run inside `detect`; no LLM is involved at run time. Inputs are the genome sequence, the candidate exons, the genetic code and the curated reference proteins already in the database.
Every step is a pure function of those inputs, so the output is reproducible and unit-testable.

### 3.1 Where it runs

A new module `src/MATPredict/detect/model_selection.py`. The pipeline calls it once per polished gene, right after `polish.classify()` returns a `PolishOutcome`
(`pipeline.py` around the `_own_model` / classify calls, **[R]** lines ~2722-2726) and before `GeneEvidence` is built. It returns a `ModelAssessment`; in
**shadow mode it only attaches the assessment to the report** and changes nothing else (section 3.6).

```
select_model(outcome: PolishOutcome, contig_seq: str, genetic_code: int,
             references: Sequence[RefProtein], policy: ModelPolicy) -> ModelAssessment
```

### 3.2 Candidate construction (pure sequence logic)

1. **Base candidates:** the exonerate model and the miniprot model of the `PolishOutcome` (each may be absent).
2. **Frame.** A model's first exon does not always start on a codon boundary (the reason `_translate_model` tries three frames). The module takes the frame from the tool when it is
   available (exonerate and miniprot both write a GFF phase; `PolishModel` does not keep it today, so Phase 1 adds `phase` to `PolishModel`, **[R]** the dataclass is `polish.py:39-67`);
   until then it falls back to: the frame with **zero internal stops over the whole model**; ties broken by the longer stop-free stretch. If no frame is stop-free, the model is `interrupted` and is **not extended**.
3. **3' extension** of a clean model (`extend_3p`):
   ```
   cds   = spliced CDS in transcript order, from the chosen frame
   rest  = len(cds) % 3                       # leftover bases complete the next codon
   walk downstream on the gene's strand from the 3' end of the last exon,
        reading codons (the leftover `rest` bases first), until
          a stop codon of `genetic_code`      -> closed, ext_aa = number of codons read
          or EXT_MAX_BP (3,000 bp) / contig end -> unclosed, no extension claimed
   accept when closed and ext_aa >= EXT_MIN_AA (30)
   new last-exon end = old end + 3*ext_aa      # stop codon excluded, as exonerate does
   ```
   Handles both strands (reverse-complement the downstream window), partial codons at exon junctions, and non-standard genetic codes (the report's `genetic_code`).
   Contig edge: if the walk reaches the contig end before a stop, record `contig_edge` and make no claim.
4. **Strength of an extension.** 30-99 aa: kept as a candidate but preferred only if an independent model (miniprot, or helixer when present) ends within `EXT_AGREE_AA` (10) of the same position.
   >= 100 aa: accepted on the open frame alone (random-ORF arithmetic: P(>= 100 sense codons in a random frame) is roughly 1-3%, P(>= 200) well under 0.1%; **[I]** GC-rich genomes are more permissive). A reviewer should check this arithmetic with the GC of the affected genomes.
5. **Candidate set:** `exonerate`, `miniprot`, `exonerate+ext3`, `miniprot+ext3` (extension only for clean 3'-open models). **Introns are never invented by extension**: a continuation stays in the same exon, same frame. A longer spliced structure has to come from miniprot (380 clean miniprot alternates are longer by a splice jump that extension cannot reproduce, **[M]**).

### 3.3 Scoring and choosing

- **One reference per target.** From the curated proteins of that gene (all clades; there are only 26 `MAT1-1-1` records, so this is cheap), the best-scoring reference of the best clean candidate is used for *all* candidates, so scores are comparable.
  (Per-candidate references made scores incomparable in the prototype, **[M]**.)
- **Score.** Raw BLOSUM62 local-alignment score (affine gaps), the same function the helixer rule already uses (`helixer_rule.query_score`, `MIN_SCORE_GAIN = 10.0`), which has the right semantics for "an extension must add aligned residues, not unaligned ones". The prototype used a MAFFT identical-columns surrogate; production should use the existing function.
- **Choice:** (1) drop every candidate with internal stops; if none remain, keep the exonerate model, flag `interrupted`; (2) take the highest score; (3) scores within `TIE_TOL` (3) are ties and the **longer** wins.
- **Not used as a gate:** alignment of the added tail to the reference. The null test shows the tail aligns only marginally better than shuffled sequence in 5 of 9 runs (Ptr ~63 identical columns vs shuffled max 56-60; **[M]** `nulltest.py`). The extension rests on the open-frame length and on independent-model agreement.

### 3.4 Flags and completeness

Computed from the chosen candidate: `terminal_stop`, `internal_stops` (positions), `starts_M`, `n_term_truncated`, `ext3_aa`, `after_stop_ref_cols`, `aligned_coverage` (reference residues inside the local alignment span / reference length; **not** the existing length ratio), `identity`, and the class
`complete | partial | remnant | interrupted`. Expected length range per order: 10th-90th percentile of the curated records of that gene in the order (fall back to the all-clade range for orders with < 3 records). A flat 190 aa floor is wrong in both directions
(curated *Pezizales* `MAT1-1-1` are 140-171 aa; *P. tritici-repentis* needs ~400 aa).

### 3.5 5' end

Phase 1 reports only: `starts_M`, `up_to_stop` (aa to the preceding in-frame stop), `up_first_M`, and `n_term_truncated`. It does **not** move the 5' end.
Reason: only 21% of `MAT1-1-1` canonical models start with M (median 11 aa open upstream to the first M), and an upstream intron can place the real start beyond the in-frame stop **[M/I]**. A 5' extension step is a Phase 2 candidate once a truth set (section 7) can score it.

### 3.6 Staging: shadow mode first

`--gene-model-policy {legacy, shadow, choose_extend}`; default `legacy` (today's behaviour byte-for-byte). `shadow` computes the assessment and writes it into a new optional `GeneEvidence.model_assessment` field and the report, and **changes no call**. Because adding the field to the report changes
the report, a test asserts that `legacy` reports are unchanged. `choose_extend` makes the chosen model the canonical one. Each downstream consumer (section 4) is switched to the new signals by its **own** policy flag after its own replay; `choose_extend` alone does not change any vote.

### 3.7 What code gets written (for review)

- `model_selection.py`: `extend_3p`, `find_frame`, `build_candidates`, `score_candidate`, `choose`, `assess`, plus constants (`EXT_MIN_AA`, `EXT_MAX_BP`, `EXT_AGREE_AA`, `TIE_TOL`, `NTERM_STOP_WINDOW`, all provisional and in one place with the evidence they were fitted on, as `homothallic.py` does).
- `polish.py`: add `phase` to `PolishModel` (default `None`); `search.py` parses it from the exonerate/miniprot GFF.
- `pipeline.py`: call in shadow mode; `GeneEvidence` gains `model_assessment: dict | None = None`.
- `report`: serialise it; `tests/detect/test_model_selection.py`: both strands, multi-exon, partial codon at a junction, contig edge, non-standard code, frameshift model, the *Botrytis* case (clean remnant stays `remnant`), the miniprot no-intron case (stops rejected), the controls (a complete gene is not extended), tie-breaking, determinism.

---

## 4. Blind-spot audit: every test that treats the canonical model as the gene

Each test below reads the canonical model's translated length, ORF fraction, or the report's `coverage` (a **length ratio**, model aa / reference aa, capped at 100: `pipeline.py:573-587`, **[R]** — it *rewards* longer models regardless of alignment) or `identity`.
The strategy is the same for all: **do not edit the thresholds; give each test the right signal (the chosen model's `aligned_coverage` and flags), switch it behind its own policy flag, and replay it.**

| # | test | where | what it assumes | failure | strategy | replay that must pass before the switch |
|---|---|---|---|---|---|---|
| 1 | homothallic completeness | `homothallic.py:48-50,89` | ORF fraction and a flat aa floor describe completeness | passes vacuously for edge-ended models; excludes true genes (1,684 of 2,418 `short`/`orf_truncated` rows would reach >= 190 aa clean after extension, an upper bound, **[M]**); a flat floor threatens short true genes | use `complete` from the chosen model; drop the aa floor; per-order range | the 3,181-genome patched set: *Pyrenophora* and *Aureobasidium* recovered, *Botrytis* still excluded, all controls at 0 |
| 2 | homothallic remnant check | `homothallic.py:78-84` | report `coverage` is aligned coverage | after extension coverage rises, so a true remnant is judged on the wrong number; a 5'-only truncation is invisible to it | apply to the chosen model's `aligned_coverage`; add `n_term_truncated` as a second remnant signal | *Botrytis* (29 loci) and *Sclerotinia* (3) still removed |
| 3 | `full_length_models` | `pipeline.py:535-571`, `HOMOTHALLIC_MIN_MODEL_COVERAGE = 0.5` | length ratio >= 0.5 means full length | an edge model reads 32-60%, straddling the cutoff; median canonical `MAT1-1-1` is 189 aa vs a 359 aa median reference (~53%) **[M]** | `aligned_coverage` of the chosen model | audit sample of 50 loci changing class |
| 4 | partial-core vote trigger | `partial_core.py:46-48,112-121` (`PARTIAL_COVERAGE 80`, `FULL_COVERAGE 90`) | a partial winner is a real partial gene | an edge model reads partial and can fire the flip; **but extension will make some winners full and turn the rule off where the flip was right** (the flip review found 41 of 46 flips correct under today's behaviour) | **highest risk.** Replay the 46 flip-review loci and the 17 Botrytis/Sclerotinia truth loci first; the profile-margin cross-check stays; no switch until the flips that were right still fire or are shown to be unneeded | flip-review: 41 of 46 right must be preserved |
| 5 | confidence cap (#100) | `confidence_cap.py:96` | partial winner means a weak call | edge-truncated correct genes are capped | switch to chosen-model flags; the cap is off by default, so lowest urgency | the `confidence-cap` evaluation (24 high-to-medium changes on 94 controls) rerun |
| 6 | strong single gene rule (#97) | `strong_single_gene.py:57,79` (`min_coverage = 90`) | a complete lone gene covers >= 90% | a correct edge-truncated lone gene fails the bar; the *Cochliobolus* case (101 uncalled genomes, one gene in the cluster) is where a correct extended model matters | switch to `complete`; re-run its 12-locus check (decision: after the clean run) | the 12-locus sequence check |
| 7 | idiomorph resolution | `idiomorph.py:50-70,669-670` | winner/loser coverage from canonical models | the *Pyrenophora* `MAT1-1-3` fragment (18% coverage, 96 bp) takes the vote over `MAT1-2-1` in 59 loci | not in Phase 1. Candidate Phase 3 change: resolve on chosen-model scores; needs its own measurement | n/a for now |
| 8 | `_translate_model` | `pipeline.py:821-880` | longest stop-free stretch = the protein | the wrong frame can win; docstring says otherwise | add a strict translator (`protein, frame, internal_stops`) used only by the new module; leave `_translate_model` alone for classifier scoring (domain scoring tolerates fragments) until a replay shows no score change | profile bit-scores unchanged on the 3,181 set |
| 9 | helixer rule | `helixer_rule.py:2-25` | core genes keep their homology model | helixer cannot fix a core-gene edge | out of scope; revisit only with the held-out test (a decision already recorded) | n/a |
| 10 | `mat_gene_gate` | `mat_gene_gate.py:68,87` | domain score | unaffected (alpha domain lies in the aligned part) | no change; confirm by replay | bits unchanged |
| 11 | `two_idiomorphs` / `split_locus` / `_short_orf_genes` | `two_idiomorphs.py:58`, `split_locus.py:117`, `pipeline.py:650-700` | identity floors, reference lengths, 30 aa | not length-of-model blind | no change | n/a |
| 12 | `report/html.py`, `track_report.py`, `scoring.py`, `benchmark.py`, `holdout.py` | — | **not audited** | unknown | assign an audit task; grep for `coverage`, `protein_aa`, `MIN_`, `orf` | — |

Ordering: 1 and 2 first (PR #117's follow-up, the motivating case), then 3, 8 (no vote change), then 5-6 (off by default), and **4 last and only with its replay**, because it changes calls the curator has already judged right.

---

## 5. Do we need a better HMM gene model, or is extension + exonerate + splicing enough?

**Short answer: extension plus miniprot/exonerate splicing is enough for the 3' problem and should come first; a gene-finding HMM is warranted only for the residual classes it cannot touch, and we should measure those before building anything.**

What each tool can and cannot do:

- **The profile HMMs already shipped** (alpha-box, HMG-box, per-gene profiles used for classification) are *domain* models. They say a domain is present and where; they carry no exon/intron structure and no start or stop signals, so they cannot decide where a protein ends. Using them as the gene model would give the same edge-ended answer. Their right use here is a second opinion on completeness (domain envelope coverage) and as an evaluation instrument.
- **Splice-aware, protein-guided gene finding** (exonerate, miniprot: already used; GeneWise with a profile; **AUGUSTUS with protein-profile extension (PPX)**; **GeMoMa** with several references at once) do carry splice-site, phase and stop/start signals. They address what extension cannot: **splice jumps** (380 clean miniprot alternates longer by a jump), **interrupted models** (142 `MAT1-1-1` and 350 `MAT1-2-1` with internal stops, where the reference aligns across the stops in some *Aureobasidium* models) and **5' starts** (79% of `MAT1-1-1` models lack M).
- **Helixer** (deep learning, already available as helixer-lite for flank rules) is an independent opinion. It only extends flanking genes today; whether it should for core genes needs the held-out test already on the revisit list.

Where the evidence points **[M/I]**: the survey says 1,108 of 1,721 `MAT1-1-1` models are 3'-open (extension fixes them), 142 have internal stops, and about 380 across both genes need a splice jump; so extension plus miniprot should resolve the large majority of `MAT1-1-1` models, but **the residual (interrupted + N-terminal + splice) is unknown until there is a truth set**.

Decision rule I propose: after Phase 1 is scored on the truth set (section 7), if the residual errors are mostly **splice or N-terminal**, add one more *candidate generator* (AUGUSTUS-PPX with a clade species model, or a GeMoMa-style multi-reference consensus) to the same chooser, never as a replacement; if they are mostly **interrupted/stop** models, check frameshift repair and read evidence (section 7e) before adding a gene finder. We should **not** build our own HMM gene finder.

A cheap lever to test first: exonerate's `protein2genome:bestfit` model (global in the query: it tries to place the whole reference protein), alongside multiple references per target rather than the single best. Both may remove the local-alignment edge at its source. **[I]** I recall the `bestfit` model variant from the exonerate documentation; it must be confirmed against the installed version before relying on it (section 8). The risk is the opposite error: forcing the full query across a divergent C-terminus produces spurious alignment, which the chooser's stop filter and score would have to catch.

---

## 6. Is it an exonerate-vs-miniprot problem or a query problem?

Mostly **design**: exonerate is always canonical, and its model is then tested with a definition of completeness it satisfies vacuously. Amplified by exonerate's local-alignment edge **[R/M]** (`search.py:855-870`: the model is the alignment, with no extension to a stop). The reference is a secondary factor: *P. tritici-repentis* aligns to *Bipolaris* C5 at 62% identity over 193 columns and the alignment edge is at reference residue 199 of 383; a closer reference would extend further, but cannot be guaranteed for new clades.
Miniprot helps (it is longer and clean in about 1,000 rows per gene) but is also edge-limited (628 of 1,349 `MAT1-1-1` alternates still sit >= 30 aa before a stop) and can skip introns. Neither tool should be treated as the default winner; the chooser picks per locus.
There are only 26 curated `MAT1-1-1` proteins (140-711 aa, median 359), so many clades have no close reference, which is the argument for section 7(b).

---

## 7. What to train or build so `MAT1-1-1` is right more often

(a) **A truth set of annotated proteins, built first.** Without it none of the thresholds can be set honestly (the existing 85/80 remnant thresholds were fitted after seeing *Botrytis*). Source: genomes in the BFD set that carry a gene annotation, with the `MAT1-1-1`/`MAT1-2-1` protein of the annotated gene at the locus; the pipeline's own `tracks` comparison and GenBank input (#88) and the held-out benchmark (`docs/holdout-benchmark.md`, records withheld at radii record/species/genus/family/order) already provide the machinery. Metrics: exact-protein match, length ratio, 3' and 5' end offsets, internal-stop rate, per clade. Caveat **[I]**: annotated proteins are themselves predictions (often exonerate/Augustus/Funannotate output) and can be wrong at the same edges; prefer annotations backed by RNA-seq or curated records, and say which.

(b) **Grow the reference roster with verified full-length proteins.** After Phase 1, chosen models that are `complete`, 3'-closed, start with M and agree with an annotation or with miniprot, in clades with no close curated reference, are candidates for new tier-3 records (curated by a person, never promoted automatically). This addresses the cause (a divergent single reference) rather than the symptom.

(c) **A calibrated chooser instead of hand thresholds.** Features: internal stops, `ext3_aa`, score margin, miniprot/exonerate/helixer agreement, aligned coverage, `n_term_truncated`, GC. Label: concordance with the truth set. Start with the rule in section 3.3; fit a small interpretable model (logistic regression or a shallow tree) only if the truth set is large enough for **leave-one-genus-out** validation, and report held-out performance, not fit.

(d) **Clade splice-site models** (position weight matrices for 5'/3' splice sites and branch region, trained from annotated introns of that clade) to score miniprot's splice jumps and exonerate's 53 bp introns, and to test whether a missed intron explains the stop-interrupted *Aureobasidium* models.

(e) **Read evidence where it exists** (the `reads` module and the SRA-derived panels): spliced RNA-seq reads across a proposed intron settle the missed-intron and 3'-end questions directly; optional, per-genome, not a requirement.

(f) **Full-length clade profile HMMs** (hmmbuild on MAFFT alignments of the extended, verified proteins per order) as a *completeness* instrument: domain-envelope coverage of the chosen model, and a check that the C-terminal region is conserved at all (the *P. tritici-repentis*/*P. teres* extensions share `...NDAFRVGADEDATLPTFDGANNA`).

---

## 8. Open questions and what I have not verified

- Whether `protein2genome:bestfit` exists in the installed exonerate and what it does with a divergent C-terminus (section 5): not tested.
- Whether exonerate and miniprot GFF phases are reliably usable as the frame (section 3.2): not read. The fallback (stop-free frame) is untested on frameshift models.
- The random-ORF probabilities (section 3.2) are arithmetic, not measured on these genomes.
- Only `MAT1-1-1` was prototyped; `MAT1-2-1` (350 models with internal stops, 855 open >= 100 aa) is surveyed but no chooser was run on it.
- The prototype score is a MAFFT identical-columns surrogate, not the pipeline's BLOSUM62 score; score differences between `exonerate+ext3` and `miniprot(+ext3)` were often 1-10 points at 31-44% identity (not decisive); ties went to the longer model by design.
- The prototype's test set (35 genomes) includes the cases the thresholds were fitted on: it demonstrates, it does not validate.
- Extensions are open reading frames, not proven coding; no independent annotation or RNA was checked.
- The *P. tritici-repentis* "ends `...LFA*`" in an earlier note is not reproduced by the canonical translation (which has no terminal stop); that note's wording should be ignored.
- Items 9-12 in section 4: not audited.
- Run-to-run noise (one *Ogataea* start coordinate, one intermittent test) is not investigated and could mask small effects in a replay.

## 9. Plan and acceptance criteria

| stage | work | done when |
|---|---|---|
| 0 | truth set (7a); baseline concordance of today's canonical models by clade | a table of baseline exact-match and length-ratio rates, with the held-out split fixed in advance |
| 1 | `model_selection.py`, `phase` on `PolishModel`, shadow mode, tests | `legacy` reports byte-identical; shadow assessments on the 3,181 set reproduce the prototype's lengths on the 35 genomes (+/- the score change) |
| 2 | scored on the truth set | held-out clade results: fraction of models 3'-complete, exact-protein rate, and no increase in wrong promotions (a model promoted to longer that disagrees with the annotation) |
| 3 | consumer switches in the order of section 4, one PR each, each with its replay | each replay's stated criterion; no change to complete-gene calls |
| 4 | residual analysis; add a candidate generator only if the residual is splice/N-terminal (section 5) | a written decision |
| 5 | roster growth (7b), splice PWMs (7d), calibrated chooser (7c) | held-out improvement |

Principles carried over from earlier decisions: defaults stay unchanged until a replay is read; thresholds are provisional and recorded with the evidence they were fitted on; the curator merges PRs; a held-out test is needed before anything is called validated.
