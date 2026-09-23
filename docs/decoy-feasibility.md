# Decoy HMG-box / alpha-box proteins: measured, and not built

**Verdict: not worth building right now.** The modelled-gene bar removed the
false positives a decoy test was meant to catch, and on the only ground-truth
set available there are none left to demote. Recorded here so the idea is not
re-proposed from scratch, and so the conditions that would make it worth
revisiting are written down.

Curator's request, 2026-09-22:

> I do wonder if we should have a decoy HMG box and alphabox proteins which are
> NOT the MAT loci and if the region has a better hit to those than the alpha or
> HMG the region is demoted?

## The idea is sound, and it is not the one already rejected

Two different cross-match problems exist, and only one of them is a decoy's job:

* **Genome-wide paralog noise** -- scattered HMG-box and alpha-box genes
  elsewhere in the genome. A decoy set addresses this directly.
* **Within-locus cross-match** -- `MAT1-1-3` vs `MAT1-2-1`, both HMG-box,
  E 2.5e-17 to each other, opposite idiomorphs. An earlier experiment showed a
  `MATa1` cutoff with **0 of 400 false positives against random genomic decoys**
  while 89 of 157 "MATa" calls still landed on an ORF whose true best match was
  `MATA2`. Decoys cannot see that, because it arises INSIDE the locus. That
  needs reciprocal-best-hit against the MAT set.

The earlier negative result is therefore not evidence against this proposal.

## The decoy pool would be cheap and plentiful

ZygoLife carries Pfam annotations for **812 genomes**
(`annotate_misc/annotations.pfam.txt`). Sampling 40 of them:

| Pfam | Name | Total | Median per genome | Max |
|---|---|---:|---:|---:|
| PF00505 | HMG_box | 989 | 24 | 35 |
| PF00046 | Homeodomain | 923 | 27 | 48 |
| PF04769 | MATalpha_HMGbox | 0 | 0 | 0 |
| PF08800 | MATA_HMG | 0 | 0 | 0 |

So ~24 non-MAT HMG-box proteins per genome are available to build from, and
none of them are annotated as MAT-specific. Decoys would add one query set to
the existing tblastn call and never polish, so the runtime cost is small.

## But the paralog tail is thinner than assumed

Diamond blastp, all 15 curated `sexM`/`sexP` proteins against the full
8,740,739-protein ZygoLife set:

* 2,313 HSPs over **1,095 distinct proteins across 812 genomes -- about 1.3
  per genome**
* bitscores: max 612, p90 152, median 64, min 45

At the PROTEIN level `sexM`/`sexP` are already specific. The pipeline's false
positives were never protein-level: they were genomic HSP fragments of 150-350
bp at 8-28% query coverage, which is a different failure and the one the
modelled-gene bar removes.

## What is actually left to demote: nothing

Ground truth: the 23-genome Zygo set, curated Plus/Minus per organism
(`testset/Zygo/{plus,minus}.abspres.csv`). Run with `--phylum Mucoromycota`.

|  | before the bar | after the bar |
|---|---:|---:|
| loci reported | 31 | 22 |
| correct idiomorph | 17 | **22** of 23 |
| wrong idiomorph | 0 | 0 |
| extra loci beyond the one real locus | 8 | **0** |

**Zero false positives.** Every reporting genome returns exactly one locus,
with textbook Mucoromycota architecture -- `tptA, sexP, rnhA, algA, glrA` --
and the correct idiomorph. A decoy test demotes a region whose best hit is a
decoy; there is no surviving region to demote.

The bar also FIXED five calls rather than merely tidying output. Six genomes
moved from `undetermined` to the correct idiomorph (*Actinomucor* NRRL A-23671,
*Circinomucor circinelloides* NRRL 22899, *Mucor alternans* NRRL A-15142,
*Mucor* sp. NRRL A-14906, A-21230, A-25970), because the spurious second locus
that was splitting the idiomorph vote is gone.

## The cost, stated plainly

**One genome of 23 was lost.** *Phycomyces blakesleeanus* NRRL 1554 is Plus in
the truth set. Its best cluster found `sexP` and `rnhA` but carried **0
modelled genes** -- neither exonerate nor miniprot could build a model for
either -- so the bar withheld it. The other two *Phycomyces* genomes in the set
are called `high` with 4 modelled genes, so this is assembly-specific, not a
species-level failure.

Net on ground truth: 17 correct -> 22 correct, one new miss. A decoy set would
not have recovered NRRL 1554: its problem is that nothing modelled, not that
something modelled wrongly.

## When to revisit

Build decoys if any of these becomes true:

1. A ground-truth set larger or more diverse than these 23 genomes shows
   surviving false positives. The Zygo set is Mucoromycota-only; Pezizomycotina
   has no comparable truth set, which is the real gap.
2. The modelled-gene bar is lowered or removed, at which point the population
   it suppresses returns and decoys become the next line of defence.
3. Recall pressure forces polish to be skipped for some candidates (the "cut 3"
   identity threshold), since a locus admitted without polishing cannot be
   judged by the bar at all and would need a different discriminator.

The reciprocal-best-hit check against the MAT set remains separately worth
doing, for the within-locus cross-match the decoys never addressed.
