# MATPredict: gene-model length, stops and the canonical-model choice (review, 2026-10-21)

Read-only review. Nothing in any repo was edited, committed or pushed. Scratch files: `/private/tmp/claude-501/-Users-jstajich-projects-MATPredict/519c0d9b-c67c-487f-bb14-2f7c2115d61e/scratchpad/work/`
(`modelkit.py` survey, `modelkit2.py` chooser prototype, `survey.tsv`, `core_annot.tsv`, `test2_out.tsv`, `nulltest.py`, `peek.py`). The scripts were copied to `~/fable_review/` on the HPCC and run with `srun -p short` only.
Repo paths below are relative to the worktree `/Users/jstajich/projects/MATPredict/.claude/worktrees/agent-a769420ca3477c529`.

Labels: **[M]** = measured by me in this session; **[R]** = read in code/notes; **[I]** = inference, not tested.

---

## 0. Corrections to the brief (things I re-checked and found different)

1. **Botrytis: the in-model stop is at residue 10, not 209.** [M, `peek.py`] The 219 aa canonical model of GCA_000227075.1 (and the other two genomes, and GCA_004379285.1) reads `LRDVDARRI*SKYFALL...`: one stop at index 9 (0-based), then 209 stop-free residues. So exonerate's model starts 9 aa *upstream* of an in-frame stop and the 209 aa stretch begins immediately after that stop with no M (`up_to_stop` 18, `up_M` 0). The miniprot model (214 aa, 2 exons, no stop) starts at exactly that post-stop `SKYFALL...` and ends at a stop (`...PLFADF`, 0 aa to the next stop). The remnant is therefore **5'-truncated (no start codon, frame opens right after a stop) and 3'-complete**, not "stop at 209". The existing note ("in-model stop at residue 209", `docs/notes/2026-10-21_homothallic-patch-measurement.md` ~line 145) should be corrected; the 209 is the *longest stop-free stretch*, which is the whole remnant. This matters because it gives a cheap *5'* discriminator between "alignment edge" and "real truncated copy".
2. **Exonerate models never contain a stop codon; miniprot models sometimes do.** [M] Over all 3,502 core MAT1-1-1/MAT1-2-1 evidence rows, `exo` ends with `*` in 0 (exonerate's CDS excludes the stop). Miniprot alternates end with a stop in 252/1349 (MAT1-1-1) and 406/1301 (MAT1-2-1). So "terminal stop present in the model" cannot be a test on the exonerate model; the usable test is "the codon immediately after the model end is a stop" (my `ext_aa == 0`).
3. **`protein_aa` in the homothallic gate is not the canonical-frame protein.** [R] `_translate_model` (`src/MATPredict/detect/pipeline.py:821-880`) returns the **longest stop-free stretch over all three frames** of the spliced CDS (lines 868-878), and the classifier input pools **both** tools' models (`pipeline.py:2849-2858`, loop over `(miniprot_model, exonerate_model)`). So a model with a frame error or an intron-in-CDS can still yield a long "ORF" from the wrong frame. This weakens `MIN_ORF_FRACTION` further (see section 5).
4. **The earlier note says the P. tritici-repentis model ends `...EDIIRVAQGK` then `LFA*`.** [R] My translation of the canonical model (196 aa) has no terminal stop, as `orf_check.py` also reports. The `LFA*` is presumably from the MAFFT input where a stop was appended or from a different model; I did not resolve this and I do not rely on it.

---

## 1. Is "longer without stops is better; longer with stops is wrong" right?

**As a ranking rule among candidate models for the same locus it is a good default, but it is not a sound definition of correctness.** It fails in both directions.

Where it is right [M]: P. tritici-repentis (196 -> 395-401 aa), P. teres (227 -> 369-372), Aureobasidium (124-201 -> 360-386), Macrophomina (144 -> 344-399). In every one, the canonical model has zero stops, the frame stays open for +142 to +237 aa, and the clean longer models land at 85-92 % of the length of the best curated reference (383-440 aa). The extension is also 3' only: the 5' side is not resolved (see 3b).

Exceptions, with what the data say:

| Exception | Consequence for the rule | Evidence |
|---|---|---|
| **A long open frame is not a gene.** A stop-free extension can run through intron + downstream sequence. | "Longer, no stops" over-extends if the real gene has an intron after the alignment edge (the in-frame stop then falls inside the intron or far downstream). | [I] Random ORF arithmetic: with ~3-4 % stop codons per codon in frame, P(open >= 100 codons) is about 1-3 %, P(>= 200) well under 0.1 %; GC-rich genomes are more permissive. Not measured empirically. A 5-30 aa extension is meaningless; 140-240 aa is strong evidence of coding. |
| **Added region need not align to the query.** | Requiring the added tail to align to the curated reference would reject genuine extensions in divergent clades. | [M] `nulltest.py`: aligning only the 3' extension (141-237 aa) to the best of 26 curated MAT1-1-1 gives 48-64 identical columns vs 34-52 (mean) / 43-61 (max) for 10 shuffles of the same extension. The signal is only marginally above shuffled in 5 of 9 runs (Ptr ~63 vs shuffled max 56-60; Macrophomina 59 vs max 50-58) and clear in P. teres (49 vs max 43-44) and one Aureobasidium (48 vs 44). So alignment support of the *added* region to distant references is weak evidence; do not use it as a gate. |
| **Introns that look like coding sequence.** | A stop-free stretch can include an unmodelled intron (short, no in-frame stop). Extending then looks "longer, clean" and is wrong in the middle. | [I] Only a spliced model (miniprot, exonerate with splice sites) catches this. 380 of 2,437 clean miniprot alternates are >30 aa longer than the exonerate model without exonerate having an open frame at its end, i.e. the longer model comes from a *splice jump* that an in-frame extension cannot reproduce [M, `core_annot` analysis]. |
| **Miniprot is not safe either: no-intron models.** | Ptr134 and Ptr902: miniprot skipped the 53 bp intron and gave 419 aa with 6 internal stops. 'Longer' would pick it unless stops are screened. | [M] `test2_out.tsv`: mini 419 aa, internal 6, score 57-67 vs 146-148 for the exonerate+extension 401 aa. |
| **Stops that are real and should be reported.** Premature-stop remnants/pseudogenes. | A stop-free *shorter* model (Botrytis 214 aa) beats the stop-containing one but is not a gene. | [M] Botrytis (see section 0.1): clean, 91-95 % identity, 60.6 % reference coverage, no M. The right call is "truncated copy", not "complete". |
| **Frameshifts / sequencing errors in a real gene.** | A real gene with one indel has a stop downstream; "has stops -> wrong" misclassifies it. | [R] The pipeline already handles exonerate frameshifts (`pipeline.py:837-850`, curator ruling 2026-10-04): translation from the aligned blocks, used only if it then has no internal stop. In my survey: 142/1721 MAT1-1-1 and 350/1781 MAT1-2-1 exonerate models have internal stops [M]; I did not test how many are repaired by the block translation (my script ignores `cds_blocks`). |
| **Stops with reference support across them.** | Some stop-interrupted models have the reference aligning well past the stop (Aureobasidium `polished_single` rows: 82-207 reference columns after the first stop), consistent with a missed intron/frameshift rather than a pseudogene. | [M] `test2_out.tsv`, `after_stop_ref_cols`. At 31-44 % identity this is suggestive only. |
| **N-terminal start uncertainty.** | Unlike the 3' end there is no stop to bound the 5' end; the first ATG is ambiguous, and an upstream intron can hide the real start. | [M] 79 % of exonerate MAT1-1-1 models do not start with M (353/1721 do); median 11 aa open upstream to the first M; 1 % have no upstream stop within 3 kb. [I] 5' extension is much less reliable than 3'. |
| Selenocysteine/readthrough | Irrelevant here. | |

**Replacement rule.** "Prefer the candidate that is (i) free of internal stops, (ii) closed at the 3' end by a real stop, (iii) best explained by the reference/profile *where it can be aligned*, and (iv) the longest among those within a score tolerance. A model that cannot satisfy (i)+(ii) is reported as `truncated`/`interrupted`, never silently promoted or demoted." Also: a clean model that starts right after a stop with no M and covers a minority of its best reference at high identity is a remnant (Botrytis), whatever its length.

---

## 2. Where does the problem come from?

Short answer: **mainly a design problem** (how the canonical is picked and then tested), amplified by **exonerate's local-alignment edges**; miniprot helps but is not sufficient alone; reference choice contributes in divergent clades.

### 2.1 Over all 3,181 reports [M, `survey.tsv`; 2,916 genomes have a modelled core MAT1-1-1/MAT1-2-1; 3,502 evidence rows]

Exonerate canonical model, in-frame continuation after its last exon (my `ext_aa`), MAT1-1-1 (n = 1,721) / MAT1-2-1 (n = 1,781):

| 3' state of exonerate model | MAT1-1-1 | MAT1-2-1 |
|---|---|---|
| internal stop(s) in the model | 142 | 350 |
| stop is the very next codon (complete) | 182 | 444 |
| stop within 1-29 aa (essentially complete, normal edge slop) | 281 | 82 |
| open 30-99 aa | 184 | 47 |
| open >= 100 aa | 924 | 855 |
| no stop within 3 kb | 8 | 3 |

So **1,108 of 1,721 (64 %) MAT1-1-1 and 902 of 1,781 (51 %) MAT1-2-1 exonerate canonical models have an open frame of >= 30 aa past their end** (2,010 rows, 1,807 genomes of 2,916). Median extension 200 aa for MAT1-1-1. 692 of the 1,108 MAT1-1-1 models are < 190 aa as modelled.

By genus (MAT1-1-1, `open` = >= 30 aa open, no internal stop): Macrophomina 250/251, Botryosphaeria 142/142, Aureobasidium 94/160 (another 63 have internal stops), Pyrenophora 76/76, Hortaea 61/63, Cladosporium 36/39, Corynespora 27/27, Phyllosticta 26/26, Didymella 21/21, Alternaria 19/159 (137 are within 1-29 aa: close to a reference, edge slop only), Parastagonospora 0/92 (92 complete), Zymoseptoria 0/41 (38 complete), Aspergillus 12/41 (28 complete). By order: Botryosphaeriales 462/464, Dothideales 104/170, Mycosphaerellales 140/291, Pleosporales 262/505. MAT1-2-1: Macrophomina 239/239, Aureobasidium 100/197, Venturia 61/63, Cladosporium 43/44, Hortaea 39/39, Phyllosticta 37/37. **The effect is clade-structured**, as expected if it follows reference similarity.

Where the miniprot alternate exists (`polished_disagree`, 1,349 + 1,301 rows):
- exonerate open >= 30 aa in 1,030/1,349 (MAT1-1-1) and 859/1,301 (MAT1-2-1) of these. Within `polished_agree` it is 13/203 and 21/247; within `polished_single` (exonerate only, no second opinion) 65/169 and 22/233. **"Disagree" is mostly "exonerate stopped early"**, not two competing gene structures.
- miniprot alternate is > 30 aa longer than the exonerate model *and* has no internal stop in 996 (MAT1-1-1) and 1,005 (MAT1-2-1) rows. Of the 996 MAT1-1-1 rows, 742 coincide with an exonerate frame open >= 100 aa; 380 (over both genes) longer-clean miniprot models are *not* explained by an open frame at the exonerate end, i.e. they use a splice jump.
- the exonerate model has an internal stop while miniprot is clean in 130 rows (40 + 90); the reverse (miniprot has stops, exonerate clean) in 147 rows (75 + 72).
- miniprot itself is edge-limited: 628 of 1,349 MAT1-1-1 alternates still sit >= 30 aa before a stop; 252 do end with a stop. (Caveat: miniprot with frameshifts may make my extension frame wrong in a few; only 7 rows have no stop within 3 kb.)

### 2.2 Exonerate vs miniprot vs references (this is where I would not blame one tool)

- Exonerate `protein2genome` is a local alignment of the best reference; it ends where the alignment to *that reference* ends (`src/MATPredict/detect/search.py:855-870`; the model returned is the alignment, with no extension to a stop or start). When the best reference is divergent at the C-terminus (P. tritici-repentis vs Bipolaris, 62 % identity over 193 columns, edge at reference residue 199 of 383), the model is the aligned part only. [R/M]
- The reference choice is per tool and per clade (`pipeline.py` per-gene best by raw tool score, `search.py:675-711`): `polish_with_*` select one alignment per gene per tool by raw alignment score; for Ptr the best reference differed between candidate models in my first run (`5016_c5` 383 aa vs `5022_unknown-1` 440 aa), which changes `coverage`. The chooser must score all candidates against *one* reference (done in `modelkit2.py`).
- `coverage` on the report is a **length ratio** (model aa / reference aa, capped at 100) (`pipeline.py:573-587`). It therefore *rewards* longer models regardless of alignment, and it penalises short-reference/long-gene pairs. It is well matched to "an extended model will raise coverage" but does not know whether the added residues are homologous.
- Reference length varies by clade [M, db]: curated MAT1-1-1: 26 records, 140-711 aa, median 359, 5 %/95 % 164/623. Pezizales has records of 140-171 aa (5 records, median 171); Mycosphaerellales 296-388. A fixed floor of 190 aa is therefore wrong for some real clades in the other direction too (it would exclude a complete short Pezizales gene).

**Verdict.** (1) Exonerate (as configured, local, no start/stop extension) systematically ends at the alignment edge: that is structural, not clade-specific bad luck. (2) The reference is a modest factor (divergent C-terminus is why the edge is early); using longer/closer references would help but cannot be guaranteed for new clades. (3) The *design* is the root: the canonical is **always exonerate when both tools produced a model** (`polish.py:116-118`), irrespective of length, stops or agreement; the downstream gates then test that model's "ORF completeness" with a definition (longest stop-free stretch relative to coding length, `homothallic.py:89`) that a model ending at an edge satisfies vacuously.

---

## 3. Proposed strategy

### 3a. Choosing among candidate models

Candidates per target: `exo`, `mini`, `exo+ext3`, `mini+ext3` (extension to the next in-frame stop; only when no stop already follows within a few aa), plus a helixer/ab-initio model when present (see "extend-only" rule below).

Procedure I prototyped (`modelkit2.py:choose`) and recommend:
1. Pick **one reference per target**: the curated protein with the best support for the best clean candidate (so candidates are compared on identical ground; per-candidate references made scores incomparable in my first run).
2. For each candidate compute `score = identical columns with that reference - 0.5 x residues with no reference counterpart` (my local pairwise MAFFT surrogate; in production use the existing helixer-rule BLOSUM62 local-alignment score `helixer_rule.py` `MIN_SCORE_GAIN = 10.0`, which already has exactly this "an extension must not replace aligned residues with unaligned ones" semantics).
3. **Hard filter**: zero internal stops (or the exonerate frameshift-corrected block translation is stop free). If no candidate passes, keep the exonerate model, flag `interrupted`, and report the first-stop position and the number of reference columns that align *after* the first stop (readthrough/frameshift/missed intron vs pseudogene). Do not promote a model from the wrong frame.
4. Among clean candidates take the highest score; scores within ~3 are a tie and the longer wins. (Longest-among-clean alone is not safe: the Botrytis 214 aa clean miniprot model is clean but is a remnant.)
5. Flag the winner: `remnant` (identity >= 85 % and reference coverage < 80 %, the existing rule in `homothallic.py:78-84`, kept) or `partial` (coverage < 80 %).

I would **not** require the added tail to align to the reference (weak signal, section 1), and I would **not** use "longest stop-free stretch in any of three frames" as a completeness measure.

### 3b. Extending or trimming the ends

- **3' end:** extend to the next in-frame stop when the codon after the model is not a stop, the open stretch is >= ~30 aa, a stop is found within ~3 kb, and the result is stop-free. For long extensions (>= 100 aa) the open frame itself is strong evidence; for 30-100 aa I would only accept it when miniprot (or a helixer model) agrees on the end within a few aa. [I] No first-principles guarantee; random ORF arithmetic only.
- **Introns:** do not try to splice during extension. Let miniprot (or exonerate) propose the spliced structure; the extension step only handles continuation in the *same frame, same exon*. If miniprot's splice-jump model is clean and longer, it is a candidate in its own right (380 such rows).
- **5' end:** extend to the first in-frame ATG after the preceding in-frame stop *only when the model does not start with M and the upstream segment to the stop is short* (survey median 11 aa to the first M; Ptr: stop 64 aa upstream, M 22 aa up; Botrytis: stop 18 aa up and no M => no extension possible). Treat 5' extension as optional/lower confidence; an upstream intron can place the real start further away than the in-frame stop. If there is no M between the preceding stop and the model start and the model begins within ~20 aa of the stop, label **`n_terminal_truncated`**. This was the right call on all three Botrytis/Helotiales remnants tested (`up_to_stop` 4-18, `up_M` 0).
- **Trimming:** when the frame has an internal stop near the model start (Botrytis, stop at residue 10), trim to after the stop (what miniprot did). Trimming at an internal stop elsewhere (e.g. residue 92 of 233 in GCA_004379285.1) would drop most of the gene; that case is `interrupted`, not trimmed.

### 3c. Stop screening

Report separately: `terminal_stop` (codon after the model, or the last codon of a miniprot CDS, is a stop), `internal_stops` (count, positions, stop-free prefix/suffix lengths), `n_terminal_start` (M present / no M but preceded by a stop), and `after_stop_ref_cols` (reference columns aligning after the first stop). Premature-stop remnants (Botrytis) are those with clean frame, 3' complete, **no M and stop upstream**, and reference coverage < 80 % at high identity. Frameshift-tolerant: keep the existing block-translation correction (`pipeline.py:837-850`), and use the miniprot frameshift-bearing model as an additional candidate rather than trusting only exonerate.

### 3d. Completeness definition downstream

Recommend: **reference coverage of the best candidate after extension/splicing, against one curated reference per target, with a per-clade expected length range** rather than a flat aa floor. Concretely: a gene is `complete` when the chosen model is stop-free, closed at the 3' end by a stop, starts with M or is flagged `n_terminal_unknown`, and covers >= 80 % (the `PARTIAL_COVERAGE`) of its best reference **by alignment** (not only by length ratio); `partial` between 50 and 80 %; `remnant` as defined. A terminal stop should be *required for "complete"* but **not** for the profile/vote tests (those only need the domain). Expected length range: curated MAT1-1-1 5-95 % is 164-623 aa overall, and clade medians are 171 (Pezizales), 353 (Helotiales), 383 (Pleosporales); I would compute a 10th-90th percentile range per order from the curated records and fall back to the global range only for orders with < 3 records. What I would not do: keep a fixed 190/200 aa floor as the sole gate.

### Recommended combination

(a) candidate set exo, mini, exo+ext3, mini+ext3 chosen by one-reference BLOSUM-type score with an internal-stop hard filter and ties to the longer; (b) 3' extension to the next in-frame stop, optional 5' extension to first ATG, `n_terminal_truncated` when stop-preceded and M-less; (c) separate flags for terminal stop, internal stops, N-term, post-stop reference support; (d) completeness = aligned coverage of the chosen model vs best reference with per-order length range. Keep the helixer rule extend-only (`helixer_rule.py:2-25`), and consider letting it apply to core genes only once a held-out test exists (memory: helixer only extends; held-out test pending).

**What I would not do:** (1) "longest model wins" with no stop filter (picks Ptr134's 419 aa/6 stops); (2) "no stops wins" with no completeness/remnant flag (picks Botrytis 214 aa clean as complete); (3) require the added tail to align to the reference; (4) treat exonerate and miniprot as symmetric votes by default (the pooled "disagree" is dominated by exonerate edge truncation); (5) extend through introns blindly; (6) tune thresholds on Aureobasidium/Botrytis alone (the existing 85/80 remnant thresholds were already fitted post hoc, per the note, and the new rules would be too).

---

## 4. Test of the proposal [M: actually run]

`modelkit2.py test` on the 17 genomes of `hl_align_genomes*.txt` plus 18 further genomes from the survey (Botrytis variants, Sclerotinia, Alternaria/Parastagonospora/Aspergillus/Zymoseptoria controls, Aureobasidium/Friedmanniomyces/Helotiales rows with internal stops, 2 Macrophomina). 35 genomes, 88 candidate rows, `test2_out.tsv`. References: the 26 curated MAT1-1-1 in `db/Ascomycota/*/*/proteins.faa`. Pairwise MAFFT `--localpair`; coverage = reference columns aligned / reference length (alignment-based, global-pairwise so inflates coverage a bit; identity is over aligned columns).

Chosen model per group (canonical aa -> chosen aa; internal stops; ref coverage; flag):

| Group (genomes) | Canonical exo | Chosen | Internal stops | Ref coverage (ref aa) | Notes |
|---|---|---|---|---|---|
| P. tritici-repentis (4) | 196 aa, cov 0.445-0.50 | exo+ext3 401, 401, mini 344, exo+ext3 395 | 0 | 0.88, 0.88, 0.86, 0.87 (440/383 aa) | miniprot 419 aa/6 stops rejected in Ptr134 and Ptr902 (score 57-67 vs 146-148). Genome 3 chose miniprot 344 (spliced, score 161) over exo+ext3 401 (145) on the 383-aa reference. Genome 4's exo+ext3 and mini are identical (395). |
| P. teres MAT1-1 (2) | 227, cov 0.585 | mini 372 | 0 | 0.914 (383) | exo+ext3 369 within 0 score points; tie -> longer |
| Aureobasidium (4) | 124-201, cov 0.32-0.52 | 386 (exo+ext3), 332 (mini), 371 (mini+ext3), 362 (exo+ext3) | 0 | 0.918, 0.820, 0.881, 0.869 | all four clean; every canonical rose from 32-52 % to 82-92 % coverage |
| Macrophomina (2) | 144, cov 0.357 | mini+ext3 399 | 0 | 0.851 (389) | |
| Botrytis remnants (3 + 1 + 1) | 219, 1 stop (at residue 10), up_M 0 | mini 214 (clean) | 0 | 0.606-0.62 at ident 0.91-0.95 | flagged `REMNANT`; 74 aa Botrytis -> 79 aa, cov 0.22, `REMNANT`. Note: the stop is N-terminal, no M, `up_to_stop` 18 |
| Botrytis full-length (1) | 348, cov 0.986 | mini 353 | 0 | 1.000 | not penalised |
| Botrytis 243 aa/5 stops (1) | 243, stops at multiple positions | kept exo | 5 | 0.688 | `interrupted+partial`; after_stop_ref_cols 237 |
| Sclerotinia (2) | 219, cov 0.567 | mini+ext3 281 | 0 | 0.728 | `partial(cov<80)`: stays flagged |
| Controls: Alternaria (2), Parastagonospora (2), Aspergillus (2), Zymoseptoria (2) | 375, 344, 368/371, 296 | mini 378 / exo 344 / exo 368/371 / exo 296 | 0 | 0.97-1.00 | **no over-extension**: complete models stayed, no candidate with a spurious +ext was chosen; extension only offered when >= 20 aa open |
| Aureobasidium/Friedmanniomyces with stops (5) | 104-147, 2-3 stops | `exo+ext3` kept as `interrupted` (4), mini clean 273 (1) | 2-3 (4) / 0 (1) | 0.30-0.70 | `interrupted+partial`; after_stop_ref_cols 79-261: reference aligns across the stops, so these are unresolved (missed intron/frameshift vs pseudogene), not "wrong because they have stops" |
| Helotiales with stops (2) | 232 (10 stops), 233 (1 stop at residue 92) | mini 242 / 237 | 1 / 1 | 0.72 / 0.67 | `interrupted`; mini also has a stop |

Summary of the 12 Pyrenophora/Aureobasidium/teres-type genomes with an open-frame canonical (4 Ptr, 2 teres, 4 Aureobasidium, 2 Macrophomina): the chosen model has 0 internal stops and 0.82-0.92 reference coverage in all 12 (Ptr 0.86-0.88, teres 0.914, Aureobasidium 0.82-0.92, Macrophomina 0.851); canonical coverage was 0.32-0.585.

Caveats on the test [M/I]:
- 35 genomes, one gene, one reference set; Botrytis and Pyrenophora are the same cases the thresholds were fitted on, so this is a demonstration, not a held-out test. The alignment scorer's identity is 31-40 % for Aureobasidium/Pyrenophora extensions (distant references), so the score differences between `exo+ext3` and `mini(+ext3)` (often 1-10 points) are not decisive. Ties are broken toward the longer clean candidate, which is by design.
- My coverage is a pairwise MAFFT surrogate, not the pipeline's length-ratio `coverage` or the helixer-rule BLOSUM62 score; I did not use pyhmmer profiles.
- The `nulltest.py` result (section 1) says the tail is only weakly supported by alignment: the extension decision rests mainly on the open frame length and on miniprot agreement.
- Not run: any test with helixer; no MAT1-2-1 chooser run; no frameshift (`cds_blocks`) repair; no reproduction of the pipeline's own polishing.

---

## 5. Downstream tests that share the blind spot (not changed)

All of these treat the canonical model's translated length, ORF fraction, `coverage` or `identity` as if it described the gene:

1. `src/MATPredict/detect/homothallic.py:48-50, 89`: `MIN_ORF_FRACTION = 0.9`, `MIN_PROTEIN_AA = 190`; `confirmed_class` (line 87-95). Longest stop-free stretch over coding length; passes vacuously for an edge-ended model, and the 190 aa floor is a threshold on aligned extent. [R/M] In `why_not.tsv` 2,418 core rows fail on `short` / `orf_truncated`; for 1,684 of them (MAT1-1-1 and MAT1-2-1) the exonerate model plus in-frame extension would be >= 190 aa with no internal stop (Macrophomina 250+239, Botryosphaeria 142, Aureobasidium 101+125, Pyrenophora 64, Venturia 61, Friedmanniomyces 40+58, Hortaea 55, Cladosporium 33+43, Phyllosticta 26+37). That is an upper bound on the effect (it is not a claim these genomes are homothallic). It is also asymmetrical: a flat 190 aa also threatens complete short genes of clades with short references (Pezizales 140-171 aa).
2. `homothallic.py:78-84`: `REMNANT_MIN_IDENTITY 85 / REMNANT_MAX_COVERAGE 80` use the report `coverage` (length ratio of the canonical). After extension coverage rises, so the remnant check must be applied to the *chosen* model; Botrytis stays a remnant (cov 60.6-62 %, ident 91-95 %), but remnants with a 5' truncation only are detectable only by `n_terminal_truncated`.
3. `src/MATPredict/detect/pipeline.py:535-571` `full_length_models`, `HOMOTHALLIC_MIN_MODEL_COVERAGE = 0.5`: canonical aa // reference aa. An alignment-edge model reads 32-60 % here; the 0.5 threshold falls inside the range. [M/I] Survey: median canonical MAT1-1-1 model is 189 aa against a median curated reference of 359 aa (~53 %); I did not tabulate the report `coverage` column itself.
4. `src/MATPredict/detect/partial_core.py:46-48` `PARTIAL_COVERAGE = 80`, `FULL_COVERAGE = 90`, and the trigger at lines 112-121 (`w.coverage >= partial`): a canonical edge model reads as partial (e.g. Ptr 51 %), which can trigger a vote flip toward the other idiomorph. The existing profile cross-check limits harm but the trigger is the same blind spot.
5. `src/MATPredict/detect/confidence_cap.py:96` (`best.coverage < partial` and no full-coverage gene) and `strong_single_gene.py:57,79` (`min_coverage = 90`): both use the report `coverage`; a correct but edge-truncated gene fails the 90 % bar; Botrytis-like clean remnants at 92 % identity are handled only by `identity`.
6. `src/MATPredict/detect/idiomorph.py:50-70, 669-670`: `winner_coverage` / `loser_coverage` are carried into the idiomorph resolution evidence; `mat_gene_gate.py:68` (`DEFAULT_MAT_GENE_MIN_SCORE = 100.0`) scores the canonical translation (longest stop-free stretch in any of three frames), unaffected by the choice if the alpha domain is inside the aligned part, which it usually is. `pipeline.py:2849-2858` pools both tools' proteins for the classifier, so the vote is more robust to this than the length tests.
7. `src/MATPredict/detect/two_idiomorphs.py:58` `MIN_PROTEIN_AA = 30` (loose; not affected in practice) and `split_locus.py:117`, `mat_gene_gate.py:87` (identity floors, which are unaffected by length). `pipeline.py:650-700` `_short_orf_genes` (60 aa floor) uses curated reference lengths, not model lengths, so it is not blind in this way.
8. `src/MATPredict/detect/helixer_rule.py:2-25`: core genes always keep their homology model, so helixer cannot fix a core-gene edge; extending that rule to core genes needs the held-out test noted in memory.

Not checked: `report/html.py`, `track_report.py`, `scoring.py`, and anything in `benchmark.py`/`holdout.py`. [I] Likely candidates for the same assumption, not verified.

---

## 6. Uncertainties and what I did not verify

- Survey assumes genetic code from the report (1 or 4 only handled specially). Exons are used as reported; exonerate `cds_blocks`/frameshift correction is not reproduced. 5' upstream is read in frame without splicing.
- The extension counts are *open frames*, not proven coding. The null test shows weak alignment support for the tail; I believe the long extensions (>= 100 aa) are real (conserved C-terminus between P. tritici-repentis and P. teres, and 342-399 aa lengths matching 383-389 aa curated records), but I did not verify them against any independent annotation or RNA.
- Botrytis conclusion (5'-truncated, M-less, stop-preceded) rests on 4 genomes (GCA_000227075.1, GCA_015147945.1, GCA_015147965.1, GCA_004379285.1) plus consistent survey rows (`up_M` 0, `up_to_stop` 4-18).
- The homothallic status of the 62 P. tritici-repentis loci is out of scope here; the extension only changes the completeness test, not the biology.
