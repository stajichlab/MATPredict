# Leave-one-out recall, measured

The first held-out generalisation measurement for this pipeline. 17 curated
records x 5 holdout radii, every run lineage-routed, scored against
coordinates NCBI supplied rather than against our own queries.

Run 2026-09-23. Inputs: `testset/coordinate_truth/`, SLURM 29040880-84.
Raw results: `results/2026-09-23_holdout/{record,species,genus,family,order}.yaml`.

## Why a holdout was needed at all

A curated record is the QUERY the pipeline searches with. Running detection on
the genome that record came from asks "did we find the thing we told it to look
for" -- self-consistency, not recall. `detect/benchmark.py` has always said so
and reported `sensitivity=None` rather than pretend otherwise. Withholding the
record makes the same genome a genuine held-out test.

## Results

| radius | hit | genuine miss | no_ref (family) | no_ref (idiomorph) | recall |
|---|---:|---:|---:|---:|---:|
| record | 13 | 1 | 1 | 2 | **93%** (13/14) |
| species | 12 | 2 | 1 | 2 | **86%** (12/14) |
| genus | 12 | 2 | 3 | 0 | **86%** (12/14) |
| family | 9 | 0 | 8 | 0 | **100%** (9/9) |
| order | 9 | 0 | 8 | 0 | **100%** (9/9) |

**Do not quote a pooled number.** Raw recall reads 76% -> 53% across the radii,
which looks like collapse under evolutionary distance. It is not: the falling
denominator is the holdout emptying the reference set, not the pipeline
failing. Every locus that still had a usable reference was found at family and
order radius.

### Two kinds of "no reference", both excluded from recall

* **family** -- the holdout removed every record of the family. A hole in the
  curated database, not a detection failure. Scoring it as a miss blames the
  wrong thing and would be fixed by curation, not by code.
* **idiomorph** -- records remained, but none for the idiomorph being sought.
  Idiomorphs are NON-HOMOLOGOUS by definition, so a reference for the other one
  cannot find this one. `mat1` holds exactly two records, *S. pombe* M and
  *S. octosporus* P; withholding either leaves only the opposite idiomorph and
  the target is structurally unfindable. Counting these as misses was this
  analysis's own first error, and it understated recall at record radius by
  12 points.

## By subphylum, at ORDER radius (nothing from the clade remains)

| subphylum | recovered | no_reference |
|---|---|---:|
| Pezizomycotina | **7/7** | 0 |
| Taphrinomycotina | 2/2 | 2 |
| Saccharomycotina | 0/0 | 6 |

**Pezizomycotina generalises.** *Aspergillus nidulans*, *A. fumigatus* (both
idiomorph records), *Botrytis cinerea*, *Sclerotinia sclerotiorum*,
*Coccidioides immitis* and *Neurospora crassa* are all recovered with their
entire order withheld. That is the question the whole exercise was built to
answer.

**Saccharomycotina cannot be tested at that radius**, because every family
holdout empties the family: `MATsc` has 6 records and all are
Saccharomycetaceae; `MTL` has exactly ONE, so *Debaryomyces hansenii* is
`no_reference` even at RECORD radius. That is a precise measurement of a
curation gap -- see [[matpredict-saccharomycotina-curation-queue]].

## The only genuine misses in the whole matrix

Two records, both *Saccharomycetaceae* `a` cassettes:

| record | radii failed | usable reference that existed |
|---|---|---|
| `28985_nrrl-y-1140_MATsc_MATa` (*K. lactis*) | record, species, genus | 5 MATsc records incl. S288C HMRa |
| `4932_s288c_MATsc_HMRa` (*S. cerevisiae*) | species, genus | 4 MATsc records |

This independently reproduces something the curation queue already suspected:
on held-out *K. marxianus* the current database finds a2 only on the SILENT
alpha cassette. The `a` idiomorph of Saccharomycetaceae is the weak spot, and
it is a reference-coverage problem rather than a detection-logic one.

## What the benchmark found that nothing else had

A routing defect, surfaced as an unexplained miss with a reference still
present (fixed in 7b870f9). `route` returned at the first matching tier, so a
SPECIES-scoped family shadowed a GENUS-scoped one:

    mat2, mat3   scope [4896]   silent cassettes   matched exactly, returned
    mat1         scope [4895]   the ACTIVE locus   never reached

Detection looked like it worked on *S. pombe* -- it reported the silent
cassettes every time -- while never searching for the locus that determines
mating type. No panel, test or review had surfaced this; every run produced a
plausible report.

## Limits, stated plainly

* **n=17**, and the five radii are not independent samples.
* The 100% at family/order rests on **9 records**, 7 of them Pezizomycotina.
* `genus` recall equalling `species` is partly a denominator artefact: two
  records moved from `miss` to `no_reference` between those radii.
* Every truth coordinate comes from NCBI's annotation of the assembly, which
  is independent of our references but not infallible -- one *Sclerotinia*
  gene model disagrees with miniprot by 280 bp at the 5' end.
* Scoring is overlap of any reported locus with the true span on the true
  contig. It does not check that the idiomorph CALL is right, only that the
  locus was found in the right place.

## What would raise the numbers

1. Curate a second family for Saccharomycotina outside Saccharomycetaceae, so
   a family holdout there leaves something to search with.
2. Add an `a`-idiomorph reference closer to *Kluyveromyces* than S288C.
3. Add a second `MTL` record, so *Debaryomyces* is testable at all.

None of these is a code change. The pipeline's measured failure mode at this
sample size is reference coverage, not detection logic.

## Correction, 2026-09-24

The tables above were scored on reports written BEFORE the routing fix
(7b870f9) they describe: the runner reused any report on disk. It also counted
any overlapping locus as a hit without checking the idiomorph, and the two
`no_reference` columns were assigned by hand. All three are fixed
(`holdout.score_holdout`, `source_commit` stamps), and the benchmark was re-run.
See `notes/2026-09-24_gap-lineage-pilots-and-holdout-rerun.md`.

What changes: there is **no genuine miss at any radius**. *S. pombe* mat1-M is
found (idiomorph M) and withheld by the modelled-gene bar, not structurally
unfindable. *K. lactis* MATa and S288C HMRa -- the "only two genuine misses" --
are bar losses on thin evidence. S288C HMRa at record radius is found with the
WRONG idiomorph (MATalpha). The Pezizomycotina order-radius result, 7/7 with
correct idiomorphs, stands.
