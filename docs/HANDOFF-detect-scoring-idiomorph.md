# Detect scoring + idiomorph calling — handoff (2026-09-20)

Branch `detect-scoring-idiomorph`, based on `detect-targeting` (PR #1).
**463 tests passing. NOT MERGED.**

Read this before touching idiomorph calling, the scoring denominator, or the
Mucoromycota reference set.

## What this branch does

It implements the fixes the `detect-targeting` validation raised, and then the
fixes that implementing *those* raised. Measured state on the 23 ground-truth
Mucoromycota genomes:

| | before | after |
|---|---|---|
| loci on the truth scaffold | 23/23 | 23/23 |
| **idiomorph correct** | **0/23** (all `undetermined`) | **23/23** |
| Plus `fraction_found` | 0.571 | **0.833** |
| Minus `fraction_found` | 0.571 | 0.800 |
| genes unsearchable | `algA`, `glrA` | **none** |
| confidence | — | 18 high / 5 medium |

And on the 44-genus discovery sweep, strict `mat_locus` calls went from
**6 genera to 14**, with `algA`/`glrA` appearing in 12 of them.

## The thing that would have bitten you

**The handoff's own recommended fix would have broken Plus detection.** It said
to resolve the sexM/sexP overlap, then "confirm `fraction_found` rises for
every genome". Every one of the 23 found exactly `{tptA, sexP, sexM, rnhA}`, so
resolving the pair alone moves Plus from 4/7 = 0.571 to **3/6 = 0.500** — the
rejection boundary, surviving only because the floor test is `<` and not `<=`.

The cause was that `algA` and `glrA` sat in the roster with **no reference
protein anywhere in `db/`**, so they could never be found and yet counted
against every genome. The denominator fix was not an improvement, it was a
prerequisite.

## Decisions, and what measured them

| decision | value | basis |
|---|---|---|
| idiomorph tiebreak | **search path**, identity as fallback | identity scored 19/23 after new references were added; search path 22/23, and 16/16 for both on the smaller set |
| overlap collapse bar | **0.5** | raised to 0.8 on the 23-genome corpus, then reverted: the sweep showed 12 candidates OVERLAPPING yet reported as two genes |
| `min_idiomorph_margin` | 5.0, per locus | Plus margins separate by 44-64 points, Minus by 2.3-13.0 |
| `max_homothallic_separation_bp` | 20 kb, per locus | the 50 kb clustering gap grouped a real sexM with a sexP fragment **47 kb** away |
| homothallic call | **both genes proteome-supported** | 54 of 75 candidates had NEITHER gene annotated, median separation 349 bp |
| min alignment length | 90 bp | a 27 bp "sexM" (nine codons) was being reported as a locus |
| `min_identity` | still **None** | Minus identities run 25.9-43.5%; any global cutoff destroys Minus detection |

Every threshold is PROVISIONAL and documented at its definition with the
measurement behind it. The `idiomorph_resolution` rows in the evidence
diagnostics are the corpus to revise them from.

## Classification, not filtering

The curator ruled on 2026-09-20 that sub-threshold calls are filed, not thrown
away. `locus_class` says WHAT was found; `detection_pass` says HOW it was
admitted. They are orthogonal.

* `mat_locus` — a core gene with at least one flank. The ordinary case.
* `homothallic_candidate` — both idiomorphs, same contig, within the
  separation bar, **both proteome-supported**. This is real biology: *Syzygites
  megalocarpus* encodes both HMG transcription factors, each flanked by its own
  intact gene with the other flank pseudogenised (Idnurm 2011).
* `idiomorph_gene_only` — idiomorph core genes with no flank at all. **Kept
  deliberately**: a lone sexM or sexP is training material for a per-idiomorph
  HMM, a search strategy this project intends to build.

Sweep result: 14 strict `mat_locus` / 28 relaxed `mat_locus` / 84
`idiomorph_gene_only` / **2 `homothallic_candidate`**.

**The one to look at: *Protomycocladus faisalabadensis* NRRL 22826**,
`scaffold_11:49110-58292` — sexP and sexM 2,571 bp apart, both from the
annotated proteome, with rnhA in the same cluster. A genuine Mucorales genus
showing the documented homothallic shape. *Spiromyces aspiralis* is the other,
but it is Zoopagomycota and was force-routed to Mucoromycota by the sweep
script, so its reference set does not really apply.

## Reference database: 6 -> 15 records

All tier-1 published deposits. Assembly-derived loci for *Cunninghamella* and
*Chaetocladium* were considered and rejected: they would be tier-2
homology-inferred, which the database design spec excludes in this phase, and a
literature search confirmed **no MAT-locus deposit exists** for those genera
(nor for *Actinomucor* or 11 others checked).

| record | accession | why |
|---|---|---|
| Mooraboolomyces wintlei | `OR965930.1` | the ONLY source of `algA` and `glrA` proteins |
| *Absidia urquhartii* ×2 | `PP971768/9` | first Cunninghamellaceae; first Minus outside Mucor/Phycomyces/Rhizopus |
| *Mucor mucedo* | `JN587498.1` | the alginate lyase the algL/algA question is about |
| *Rhizopus azygosporus* | `MG967659.1` | second `glrA` |
| *Blakeslea trispora* | `HG939558.1` | sexM only; its tptA (38 aa) and rnhA (82 aa) are fragments |
| *Syzygites megalocarpus* ×2 | `JN112239/40` | **first homothallic reference**; third `glrA` |
| *Parasitella parasitica* | `KY081664.1` | Parasitellaceae Minus |

Now 9 Minus / 6 Plus, 2 homothallic. Every gene validates at 100% identity AND
100% coverage.

**Reference balance is not neutral.** Adding two Plus references flipped four
*Cunninghamella* Minus genomes to Plus. Check the per-idiomorph count before
and after any ingest.

## Open

* **`Syzygites` sp. MES_3091 is not resolved.** A homothallic genus whose loci
  still come out `idiomorph_gene_only`. The published *S. megalocarpus* genes
  are now references, and `scaffold_132`'s sexM matches one at 100% identity,
  but no locus there clears the strict bar.
* **Non-Mucoromycota genera in the sweep.** The script forces
  `--phylum Mucoromycota` on all 44, including *Alternaria* (Ascomycota) and
  *Basidiobolus*/*Conidiobolus* (Entomophthoromycota). Their results are
  meaningless by construction. Re-run those under correct routing before
  drawing conclusions about them.
* **Resolution runs pre- AND post-polish.** Necessary, because polishing moves
  coordinates. But the pre-polish pass compares pre-polish identities, so a
  call could in principle be made on numbers the report does not show.
* **`min_identity`** stays None, to be calibrated on the diagnostics corpus.
* **Schulz et al. 2016** (`Endocytobiosis Cell Res` 27(4):39-57) has no DOI or
  PMID and is not online; the curator supplied the PDF. Its Table 1 lists 24
  species with protein IDs for all five genes — the richest remaining source of
  curation leads. See `docs/notes/2026-09-20_algA-glrA-mucorales-literature.md`.

## Reproducing

Ground truth (23 genomes) and the sweep (44 genera) both run from
`$SCRATCH/zygo/fasta/`. `testset/Zygo/` is **untracked in git**, so it does not
exist in a worktree — read it from the main checkout.

```
matpredict detect --genome <g>.fna --proteins <g>.faa --phylum Mucoromycota \
  --out-dir <out> --evidence-diagnostics <out>/evidence_diagnostics.jsonl
```

Never run a real experiment against a tree that is being edited: `detect` reads
both `src/` and `db/`, so a curation ingest mid-sweep makes different genomes
see different reference sets.
