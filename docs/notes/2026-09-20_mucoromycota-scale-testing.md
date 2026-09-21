# Mucoromycota MAT-locus detection at scale — 283 BFD genomes

Run 2026-09-20 on branch `mucoromycota-scale-testing`, based on
`detect-scoring-idiomorph` @ `505437d` (v0.5.0). Baseline before any change:
**463 tests passing**, reconfirmed in this worktree.

All numbers below are measured from real runs. Inferences are marked
`[INFERRED]`.

## Branch-point correction

The task asked for `git reset --hard detect-targeting`, then to verify the
result contains PR #2 (`34117e9`) or v0.5.0 (`505437d`). It cannot contain
either: `detect-targeting` is `e129280`, which is an **ancestor** of
`505437d`. Resetting there would have discarded the entire scoring/idiomorph
branch this work depends on. This branch is therefore based on
`detect-scoring-idiomorph` @ `505437d`, which satisfies the stated
verification. `34117e9` exists in the repo but is not reachable from
`505437d`; it is a GitHub-side merge of the same branch.

## What was run

| experiment | genomes | mode | wall time |
|---|---|---|---|
| 1. Paired mode comparison | 23 ground-truth Mucoromycota | annotated AND genome-only | 13 min + 104 min |
| 2. BFD sweep | 283 BFD Mucoromycota | genome-only | 13 SLURM jobs, 8-45 min each |
| 3. BFD paired subset | 145 of the 283 | annotated (NCBI proteomes) | 1 job, 17 min |

Zero crashes. 283/283 and 145/145 reports written; every batch log reads
`ran N genomes, 0 failed`.

Outputs (not committed): `/bigdata/stajichlab/jstajich/mucoro_scale_out/`
(`bfd_runs/`, `annotated_runs/`, `proteomes/`); analysis scripts and CSVs in
`/scratch/jstajich/28365455/mucoro_scale/`.

---

## Finding 1 — genome-only loses NO accuracy on the ground-truth set

This was the gating question. The prediction in the task brief was that
genome-only mode would degrade the idiomorph call to ~19/23, because
`idiomorph._rank` prefers a proteome-found hit and every genome-only hit is
`tblastn_genome`, so the tiebreak degenerates to identity.

**It does not degrade at all.** Both arms run on v0.5.0, same 23 genomes,
same `--phylum Mucoromycota`:

| metric | annotated (`--proteins`) | genome-only |
|---|---|---|
| locus on truth scaffold | 23/23 | **23/23** |
| idiomorph correct | 23/23 | **23/23** |
| `locus_class = mat_locus` | 23/23 | 23/23 |
| confidence | 18 high / 5 medium | 18 high / 5 medium |
| gene recovery | algA 23, rnhA 23, tptA 23, glrA 10, sexP 16, sexM 7 | **identical** |
| total loci reported | 23 | 24 |
| mean wall time | 33.6 s | **271.9 s (8.1×)** |
| max wall time | 51 s | 381 s |

Methods confirmed: every genome-only locus used only `tblastn_genome` and
`exonerate_refine`. No diamond path ran.

**Why the prediction failed.** The "identity scores 19/23" measurement was
identity applied to the *annotated* run's identities. Genome-only does not
just remove the tiebreak, it also changes the numbers being compared: sexM
and sexP are then both scored by tblastn against the genome, which is a
symmetric comparison. Three *Cunninghamella* Minus genomes make this visible
— in the annotated arm their idiomorph margins are **negative** (−3.546,
−1.721, −1.721), i.e. identity alone would have called them Plus and only the
search-path tiebreak saved them; in the genome-only arm the same three have
**positive** margins (+4.533, +1.633, +3.184) and are called correctly
without any tiebreak. [INFERRED] mechanism; the margins are measured.

**Cost side.** 8.1× runtime, spent almost entirely in polishing: genome-only
has no proteome to pre-filter candidates, so far more (cluster, family) pairs
reach exonerate. And one extra locus appeared: *Absidia cuneospora*
`scaffold_57:56138-122856`, a **66.7 kb** span with three genes at 28.4–39.0%
identity, reported as `mat_locus`, `strict`, `medium`. That is a false
positive the annotated arm did not produce.

**Conclusion: BFD-scale genome-only testing is worth doing.** That is what the
rest of this note reports.

---

## Finding 2 — annotated Mucoromycota proteomes ARE obtainable, for half of BFD

Measured, not assumed:

* 283 BFD Mucoromycota; only **14 are RefSeq (`GCF_`)**, 269 are `GCA_`.
* Querying NCBI `datasets summary` for all 283: **146 carry an
  `annotation_info` block** (14 GCF + 132 GCA). 137 have no annotation.
* `datasets download genome accession <acc> --include protein,gff3` yields
  `protein.faa` + `genomic.gff`. Joining CDS features by `protein_id` and
  taking min-start/max-end per protein produces the
  `contig:start-end:strand` defline `search.PROTEOME_DEFLINE_FORMAT`
  requires. 145 of 146 converted; 1 failed
  (`GCA_900175165.2_FCH_5_7`, no protein/gff in the package).
* **The BFD assembly FASTA uses NCBI's own sequence accessions as contig
  names** (verified on `GCF_001638985.1`: `NW_017265134.1` etc. match the
  GFF3 seqids exactly). So no renaming or remapping is needed — the BFD
  assembly and the NCBI-derived proteome can be handed to `detect` together.

Script: `/scratch/jstajich/28365455/mucoro_scale/fetch_proteomes.py`. This is
the path to make `homothallic_candidate` reachable on BFD assemblies at all;
it requires both genes proteome-supported and so can never fire genome-only.
**Zero `homothallic_candidate` calls in the 283-genome genome-only sweep**,
exactly as predicted.

---

## Finding 3 — the relaxed pass is the scale failure mode

283 genomes, genome-only. 687 loci total.

| | n | median loci/genome | mean | max |
|---|---|---|---|---|
| best call from **strict** pass | 232 | **1** | 1.44 | 6 |
| best call from **relaxed** pass | **50** | **8** | 7.04 | **13** |

A haploid fungus has one MAT locus. Strict behaves: 170 of 283 genomes (60.1%)
report exactly one locus. Relaxed does not: those 50 genomes average seven
loci each and contribute most of the 43 genomes reporting ≥5 loci.

Per-locus identity confirms the split:

| pass | n loci | median best identity | median span |
|---|---|---|---|
| strict | 335 | **77.6%** | 16,067 bp |
| relaxed | 352 | **38.1%** | 5,538 bp |

Worked example, *Rhizomucor miehei* CAU432 (`GCA_000611695.1`): **13 loci,
9 of them `mat_locus`**, every single one `relaxed`/`medium`, best identities
21.7–46.9%, spans to 53 kb, and the idiomorph calls contradict each other
across the genome (3 Plus, 4 Minus, 6 undetermined). There is no true locus
in the output; the relaxed pass has converted "nothing confident found" into
a dozen confident-looking calls.

**It is not random which genomes fall to relaxed.** By family:

| family | n | relaxed | rate |
|---|---|---|---|
| Endogonaceae | 4 | 3 | 75% |
| Lichtheimiaceae | 30 | 18 | **60%** |
| Rhizopodaceae | 106 | 22 | 21% |
| Cunninghamellaceae | 22 | 4 | 18% |
| Mucoraceae | 66 | 1 | **2%** |
| Umbelopsidaceae | 14 | 0 | 0% |

Lichtheimiaceae has **no curated reference record**; Mucoraceae is where the
curated set is densest. Within Rhizopodaceae the split tracks species:
*R. arrhizus* 7/32 relaxed and *R. microsporus* 5/10, but *R. delemar* 0/8 and
*R. stolonifer* 0/2. [INFERRED] the relaxed rate is a proxy for reference
distance, not for genome quality.

**Recommendation.** The relaxed pass should be capped — emit at most the best
N relaxed loci per (genome, family), or refuse to emit `mat_locus` from the
relaxed pass at all and let those calls carry a weaker class. Reporting nine
mutually contradictory `mat_locus` calls for one genome is worse than
reporting none, because a downstream consumer cannot tell which to believe.
This could not be seen on 23 genomes: all 23 reported exactly one locus.

---

## Finding 4 — `mat_locus` is emitted with no MAT gene in it (fixed, TDD)

`classify_locus`'s docstring says `mat_locus` is "a core gene with at least
one flanking gene". The code's final statement was `return LOCUS_CLASS_MAT`
as an unguarded catch-all, and `if not live: return LOCUS_CLASS_MAT` at the
top did the same for an empty cluster.

Measured on the sweep: **11 of 534 `mat_locus` calls contain ZERO `core_MAT`
gene, and all 11 were admitted by the STRICT pass**, so nothing downstream
flags them. Examples:

| species | locus | genes |
|---|---|---|
| *Umbelopsis ramanniana* AG | `NW_026252103.1:214691-218060` | rnhA\|algA\|btbA |
| *Cokeromyces recurvatus* | `NW_026251440.1:50778-71656` | tptA\|algA\|btbA |
| *Rhizomucor pusillus* | `FWWN02000524.1:48009-116058` | rnhA\|algA\|btbA |
| *Actinomucor elegans* | `BCHK01000001.1:185393-231177` | tptA\|rnhA\|btbA |

Fixed test-first in this branch. Two failing tests were added to
`tests/detect/test_locus_class.py` and watched fail on `mat_locus`, then
`classify_locus` gained a `flanking_gene_only` class — the exact mirror of
`idiomorph_gene_only`, filed rather than discarded per the curator's
2026-09-20 ruling. **463 → 465 tests, all passing, no regressions.**

**The class name is provisional and needs a curator ruling.** These clusters
are not junk: the flanking neighbourhood is precisely where a missed core
gene would sit, so they are leads for the short-ORF problem.

This fix changes only classification. It was made AFTER every run in this note
completed, so no number here is affected by it; re-running would move 11 calls
out of `mat_locus`.

---

## Finding 5 — genome-only and annotated disagree on unselected genomes

The 23 ground-truth genomes agreed 23/23 (Finding 1). On the 145 paired BFD
assemblies — the same assembly run both ways — they do not:

| agreement | n/145 |
|---|---|
| same primary contig | 100 (69.0%) |
| same idiomorph call | 110 (75.9%) |
| same locus count | 78 (53.8%) |
| same gene set | 56 (38.6%) |

**22 of 145 (15.2%) are flat Plus↔Minus contradictions** (neither side
`undetermined`), and **11 of those 22 have BOTH sides admitted by the strict
pass**. Idiomorph is a binary biological fact, so in each of those 11 at least
one mode is wrong and neither is hedging.

The direction is strongly asymmetric: of the 22 contradictions, **18 are
genome-only=Minus vs annotated=Plus**, only 4 the other way. [INFERRED] this
is the search-path tiebreak plus reference imbalance: in annotated mode a
proteome-found sexP outranks a tblastn-rescued sexM regardless of identity,
pushing toward Plus; genome-only removes that and the 9-Minus/6-Plus curated
set pulls the other way. The handoff already warns "reference balance is not
neutral" — this measures it on 145 genomes.

Mode also changes how often a genome reaches strict at all, and not in the
direction one might guess: **genome-only 117/145 strict vs annotated 98/145**.
tblastn finds genes the submitted annotation missed, which is the known
annotation-gap problem for small MAT genes.

**The 23-genome corpus cannot detect any of this.** It was built from a
cblaster analysis, so by construction it contains only genomes whose locus is
intact and co-located on one scaffold. It is a recall test, not a
representative sample.

---

## Finding 6 — locus spans are not bounded

687 loci: median span 13,565 bp, p90 57,085 bp, **max 154,355 bp**. 94 loci
span >50 kb and 8 span >100 kb; all 94 of the >50 kb ones are class
`mat_locus`. At the other end, 94 loci span <1 kb — including
*Rhizomucor miehei* `KK100165.1:81475-82681`, which reports **sexP, sexM and
glrA inside 1,206 bp**. Three genes cannot fit there; those are overlapping
HSPs of one HMG region plus a spurious glrA that the 0.5 overlap-collapse bar
did not merge.

For reference, the 23 ground-truth loci span 6,795–13,089 bp. A per-locus
plausible-span bound (curation data, like `max_homothallic_separation_bp`)
would remove both tails. **Not implemented** — the bound is a curator
decision, not something to set from this sweep.

---

## Finding 7 — behaviour outside the curated set (all 15 records are Mucorales)

| order | n | pass | loci/genome | best identity (min/median/max) |
|---|---|---|---|---|
| Mucorales | 264 | — | — | 28.3 / **85.8** / 100.0 |
| Umbelopsidales | 14 | **14/14 strict** | 3.07 | 34.0 / **63.6** / 71.4 |
| Endogonales | 5 | 1 strict, 3 relaxed, 1 no call | 1.40 | 27.5 / **45.2** / 56.2 |

**Umbelopsidales degrades gracefully.** All 14 reach strict at a sensible
63.6% median identity. Idiomorph splits 6 Plus / 6 Minus / 2 undetermined, so
there is no obvious one-sided bias. Caveat: 3.07 loci per genome is well above
the one a haploid genome should have.

**Endogonales does not.** *Bifiguratus adelaidae* (`GCA_002261195.1`) is the
clearest failure: **4 strict `mat_locus` calls on 4 different contigs**, two
Plus and two Minus, identities 22.2–56.2%, spans to 78.6 kb — and the
first of them carries **no core MAT gene at all** (Finding 4). The three
*Jimgerdemannia* genomes are the largest in the set (52–77 MB gzipped) and
produce the sweep's only zero-locus genome (`GCA_003990745.1`, 76.6 MB) plus
two relaxed `idiomorph_gene_only` calls at 27.5–35.6% identity.

Endogonales should be excluded from Mucoromycota routing, or given its own
curated references, before any of its calls is believed.

---

## Operational findings

* **`scripts/run_detection_batch.slurm` needs its batch JSON on SHARED
  storage.** The first 13-job submission failed instantly, all with
  `FileNotFoundError` on the batch JSON, because it had been written to
  `$SCRATCH` per the project's own HPCC guidance. `$SCRATCH` is node-local and
  the job lands on a different node. The same applies to `--out-dir` and to
  any proteome. Worth a line in the slurm script's header, which currently
  documents `$SCRATCH` only for decompression.
* **Job sizing.** `estimate_seconds_from_gzipped_size` is calibrated at
  96.4 s per 12 MB gzipped. Real genome-only runtime on this corpus is about
  **2.4× that**. Passing `target_seconds_per_job=2250` produced 13 batches of
  9–27 genomes that ran 8–45 min — within the 1–1.5 h target. The estimator
  should be recalibrated, and it should know whether a proteome is supplied:
  annotated mode averaged **53.3 s/genome** (max 121 s) against genome-only's
  ~235 s.
* **No resource pathologies.** No OOM, no timeout, no crash across 428 genome
  runs. The largest genome (76.6 MB gzipped) completed.

---

## Recommendations, in priority order

1. **Cap or reclassify relaxed-pass output** (Finding 3). Nine contradictory
   `mat_locus` calls in one genome is the single worst behaviour at scale.
2. **Curator ruling on `flanking_gene_only`** (Finding 4) — name, and whether
   these clusters should be surfaced as short-ORF leads.
3. **Build the annotated BFD proteome set into the pipeline** (Finding 2).
   146 of 283 are available and the contig names already match. This is the
   only route to `homothallic_candidate` on BFD assemblies.
4. **Get ground truth for unselected genomes** (Finding 5). The 15.2% flat
   contradiction rate cannot be resolved without it; the existing 23-genome
   set is structurally unable to expose it.
5. **Add a per-locus plausible-span bound** (Finding 6), as curation data.
6. **Exclude or separately reference Endogonales** (Finding 7).
7. **Recalibrate `estimate_seconds_from_gzipped_size`** and document the
   shared-storage requirement in the slurm header (Operational).

## What this note does NOT establish

There is **no ground truth for the 283 BFD genomes**. Every statement about
them is either internal consistency (locus counts, spans, class contracts) or
disagreement between two modes. Where the two modes contradict each other,
this note does not claim which is right. The only genomes with an external
answer are the 23, and on those both modes are 23/23.

---

## Addendum — the 11 `flanking_gene_only` candidates, catalogued

Curator ruling 2026-09-20: keep the name `flanking_gene_only`; these are
biologically interesting and should be catalogued; we need the region SIZE and
what OTHER genes sit inside it; `idiomorph_gene_only` and `flanking_gene_only`
both stay as classes, for later synteny / dotplot / clinker comparison to
separate real evolution from assembly artefacts.

| species | region | span | MAT genes hit | NCBI proteome |
|---|---|---|---|---|
| *Rhizomucor pusillus* | `FWWN02000524.1:48009-116058` | 68,049 | rnhA\|algA\|btbA | no |
| *Umbelopsis ramanniana* | `CDSBDH010000019.1:45571-94613` | 49,042 | rnhA\|algA\|btbA | yes (24 genes inside) |
| *Actinomucor elegans* JCM 22485 | `BCHK01000001.1:185393-231177` | 45,784 | tptA\|rnhA\|btbA | no |
| *Actinomucor elegans* ASM2602732 | `JAMSLZ010000005.1:1106792-1152570` | 45,778 | tptA\|rnhA\|btbA | no |
| *Umbelopsis vinacea* | `CDSBDG010000020.1:593788-633492` | 39,704 | rnhA\|glrA\|btbA | yes (13 genes inside) |
| *Radiomyces spectabilis* | `NW_026251926.1:1569571-1599194` | 29,623 | tptA\|glrA\|btbA | yes (9 genes inside) |
| *Umbelopsis* sp. AD052 | `JAIXMS010000011.1:566926-594549` | 27,623 | tptA\|rnhA\|algA\|btbA | yes (15 genes inside) |
| *Cokeromyces recurvatus* B5483 | `JNEH01002377.1:14730-35757` | 21,027 | tptA\|algA\|btbA | no |
| *Cokeromyces recurvatus* | `NW_026251440.1:50778-71656` | 20,878 | tptA\|algA\|btbA | yes (11 genes inside) |
| *Bifiguratus adelaidae* | `MVBO01000034.1:16948-29301` | 12,353 | tptA\|glrA\|btbA | yes (3 genes inside) |
| *Umbelopsis ramanniana* AG | `NW_026252103.1:214691-218060` | 3,369 | rnhA\|algA\|btbA | yes (3 genes inside) |

All 11 are `strict`/`medium`, best identities 23.2–56.2%. **7 of 11 already
have an NCBI proteome** built in this work, so their gene content is available
now; the other 4 would need a targeted Augustus run.

**Two observations worth following up.**

1. **`btbA` is in 11 of 11.** Across all 687 loci `btbA` appears in only 133
   (19.4%), so its presence in every flanking-only candidate is not chance.
   [INFERRED] `btbA` is the flanking gene most able to anchor a cluster on its
   own. Note the opposite-sign fact: loci containing `btbA` have a *higher*
   median best identity (97.6%) than those without (40.5%), so `btbA` is not
   simply a promiscuous low-quality hit — it is doing both jobs.

2. **Two of them reproduce across independent assemblies.**
   *Actinomucor elegans* appears twice from unrelated assemblies with the same
   gene set and spans agreeing to **6 bp** (45,784 / 45,778).
   *Cokeromyces recurvatus* likewise, spans 21,027 / 20,878 with the same
   `tptA|algA|btbA`. A shared assembly artefact would not reproduce this way.
   These two species are the best starting points for the synteny comparison.

**Empty-cluster branch: measured, never reached.** Of 687 reported loci, zero
had empty `gene_evidence` and zero were fragmented/multi-segment. The
`if not live` guard is defensive only, so routing it to `flanking_gene_only`
changed no observed output.

**Still open:** the spans run 3.4 kb to 68 kb against ground-truth loci of
6.8–13.1 kb, so some of these regions are probably over-wide clusters rather
than eroded loci. Deciding the span bound (Finding 6) before curating these
would change which of the 11 survive.

---

## Addendum 2 — the btbA conflation, and a correction to Finding 4

### Correction: the first fix did not work

The fix committed in `2f21259` **did not reclassify any of the 11 cases it was
written for**, and the commit message overstated it. Verified by calling
`classify_locus` with the real `db/Mucoromycota/order.yml` family against the
six real gene sets: 6 of 6 still returned `mat_locus`.

Cause: `classify_locus` split its genes on `present_in_idiomorphs`
(idiomorph restriction), not on `role`. **`btbA` is `flanking_variable` yet
carries `present_in_idiomorphs: ["Plus"]`**, so it landed in `restricted`
beside sexP/sexM. The new `if unrestricted and not restricted` branch
therefore could never fire for a cluster containing `btbA` — and all 11
candidates contain `btbA`. The test that passed used a fixture family with no
`btbA` in it, so it never exercised the real shape.

### The same conflation caused a second, larger bug

The homothallic loop iterated `restricted` and paired genes whose idiomorph
sets are disjoint. `btbA`(Plus) paired with a lone `sexM`(Minus) satisfies
that, so a flanking gene could stand in for a second core gene.

Measured on the 145-genome annotated BFD run: **17 of 18
`homothallic_candidate` calls have only ONE core gene** (sexM) and rest on
`btbA`. Every hit was `diamond_proteome`, so the both-genes-proteome-supported
gate did not stop them. Example: `GCA_011763815.1`,
`JAANIU010000625.1:2001-15277`, genes `sexM|rnhA|glrA|btbA`.

Only **1 of 18** is a genuine two-core-gene call: *Radiomyces spectabilis*,
`NW_026251940.1:262111-274785`, sexP + sexM both from the proteome.

The genome-only sweep produced 0 homothallic calls, so this is invisible
without a proteome — which is why the 283-genome sweep alone did not find it.

### Fix

`classify_locus` now splits on `role == "core_MAT"`, not on idiomorph
restriction. Both the homothallic test and the flanking-only test are
role-based. Three regression tests added against a fixture that models the
REAL family including `btbA`. Re-verified: 6 of 6 real gene sets now return
`flanking_gene_only`, and the *Radiomyces* two-core-gene case still returns
`homothallic_candidate`. **474 tests passing.**

### Where btbA sits, and whether it is required

Gene order, read from the curated records. `btbA` appears in exactly **one**
of the 15: `64495_cbs346-36_MAT_Plus` (*Phycomyces blakesleeanus* CBS 346-36):

```
tptA(+) -> btbA(+) -> sexP(+) -> rnhA(+)
```

Against the consensus Plus architecture from the other records
(`algA -> tptA -> sexP -> rnhA -> glrA`), `btbA` sits **between tptA and the
core gene**, i.e. inboard of tptA, immediately outboard of sexP:

```
algA — tptA — [btbA] — sexP/sexM — rnhA — glrA
```

Is it required? **No.** Measured across the 687 sweep loci:

| | n | btbA present |
|---|---|---|
| Plus, all loci | 220 | 57 (**26%**) |
| Minus, all loci | 207 | **0 (0%)** |
| Plus, strict `mat_locus` | 134 | 39 (29%) |
| Minus, strict `mat_locus` | 127 | **0 (0%)** |

So `btbA` is **perfectly specific but weakly sensitive**: it is never seen in
a Minus call, and seen in barely a quarter of Plus calls. Presence argues
Plus; absence argues nothing. Caveat: `btbA` has only ONE reference protein
in the whole curated set, so the 26% is a floor, not a property of the gene.

**Consequence that needs a ruling.** `expected_genes_for_idiomorph` includes
every gene whose `present_in_idiomorphs` matches, so a Plus locus is scored
against **6** expected genes and a Minus locus against **5**. Since 74% of
Plus loci lack `btbA`, most Plus calls carry a `fraction_found` penalty that
no Minus call can incur. [INFERRED] this is an asymmetric bias against Plus
detection, and it is suggestive that the handoff's measured Plus
`fraction_found` is exactly 0.833 = 5/6.

The data supports keeping `present_in_idiomorphs: ["Plus"]` on `btbA` (0/207
in Minus). What it does not support is `btbA` counting as an EXPECTED gene
for Plus. That distinction — idiomorph-specific but optional — does not exist
in the schema today. **Not changed; needs a curator ruling.**

---

## Addendum 3 — span flag and Augustus (both curator-ruled, both implemented)

**Span: flagged at 200 kb, never dropped.** Curator ruling 2026-09-20 —
large loci are real and consistent with unpublished findings by a former
graduate student. Implemented as `max_plausible_locus_span_bp` on the family
(curation data, like `max_cluster_gap_bp`), a public
`pipeline.span_exceeds_plausible_bound()` predicate, and a
`span_exceeds_plausible_bound` field emitted on **every** locus in the report
so "not flagged" is distinguishable from "output predates the field".

Measured consequence: **0 of 687 loci exceed 200 kb** (widest observed
154,355 bp). This is a guard-rail against a runaway cluster, not a filter on
present output. The 120 kb first discussed would have flagged 3, one of them a
*Blakeslea trispora* call with six genes at 87.2% identity sitting 473 bp over
the line.

**Augustus added to `pixi.toml`** (bioconda, resolves to 3.5.0, verified in
the environment). For ab initio prediction over a sliced locus window on the
137 of 283 BFD Mucoromycota with no NCBI annotation.

Augustus rather than Helixer, for packaging reasons and not accuracy.
`helixerlite` (PyPI 25.5.27, nextgenusfs) is the better annotator but cannot
join this environment: it requires `tensorflow>=2.6.2`, `tensorflow-addons`
(archived/EOL) and `keras<3.0.0`, and ships wheels for cp39/cp310/cp311 only
while this workspace resolves to Python 3.14. It is not in bioconda, and no
Helixer image exists in `/bigdata/stajichlab/shared/singularity`. Adopting it
would require a second environment or a container runtime that this
standalone tool does not otherwise need. Augustus installs like every other
binary `detect` already shells out to.
