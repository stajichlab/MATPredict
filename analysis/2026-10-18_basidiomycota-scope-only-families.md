# Basidiomycota scope-only families: redHD, redPR, rustHD, wallMAT on every in-scope BFD genome (2026-10-18)

Code: current `origin/main` at 60d7651 with default `detect` parameters and lineage routing (`detect --taxid`, one run
per genome); `scripts/family_performance_table.py` (new, this branch) and `scripts/run_asco_holdout.py` (unchanged).
Data and tables: `results/2026-10-18_basidio_scope_only/`. Curator: J. Stajich.

Descriptive only. Nothing here has a truth for the 433 BFD genomes except the 12 genomes that carry a curated record;
"called" means the pipeline returned a locus, not that the locus is right. What is measured is marked as such; judgement
is marked as judgement.

## Short answer

- **redPR (Sporidiobolales): works as far as it can be tested.** Called in 244 of 251 genomes, all `high`, one A1 or A2
  label per locus. Leave-own-record-out: 5/5 correct at record, species and genus radius (the 5 records that have a
  genome assembly). No reference exists at family or order radius, so nothing beyond the order-level scope is tested.
- **redHD (Sporidiobolales): the locus is found, the allele is not callable, and recovery was only testable inside
  *Rhodotorula*.** Called in 246 of 251 genomes. All 252 loci have `HD1` and `HD2` on one contig, all 252 divergently
  transcribed, intergenic 0 to 1,295 bp (median 490). Leave-own-record-out: 5/5 locus found at record and species
  radius. All five curated redHD records are *Rhodotorula*, so the genus radius leaves nothing: recovery in other genera
  is not measured by hold-out. The 32 *Sporobolomyces* and *Rhodosporidiobolus* genomes (called 32 of 34) rest on
  application with the arrangement check only.
- **rustHD (Pucciniales): works for *Puccinia* and *Phakopsora*, weaker or absent elsewhere.** Called in 107 of 131
  genomes (82%): *Puccinia* 42/51, *Phakopsora* 50/51 (one taxid, so not independent), *Melampsora* 5/8, *Hemileia* 0/5,
  *Uromyces* 1/4. Leave-own-record-out: 4/4 locus found at record, species, genus and family radius; no reference at
  order radius. Alleles are not callable by design (pattern vocabulary).
- **wallMAT (Wallemiales): 51 of 51 called, but the family has two records (one per version), so no held-out test of the
  version call is possible.**
- **Rhodotorula question:** in BFD, by NCBI taxonomy, every genome of *Rhodotorula*, *Rhodosporidiobolus* and
  *Sporobolomyces* (250 genomes) is inside Sporidiobolales. **None is outside**, so scope-only routing misses none of
  them by name. Related genera that were once *Rhodotorula* (*Cystobasidium*, 16 genomes) sit in other orders and are
  correctly not searched for redHD/redPR; whether they carry a red-yeast-type locus is unknown (no record).

## 1. What ran

Genomes are BFD's `input_clean_genomes/<ASM_ACCESSION>.fa.gz`, selected by lineage (taxonkit on the BFD NCBI taxids,
run under Slurm) containing the family's `taxonomic_scope` taxid: Sporidiobolales 231213 (251 genomes), Pucciniales 5258
(131), Wallemiales 431958 (52; one has no FASTA in `input_clean_genomes`, so 51 ran). The BFD CSV `ORDER` column agrees
with the lineage for all of them (no disagreement either way). One `detect --genome --taxid` per genome, default
parameters. In all 433 reports `routing_mode` is `lineage`, `routing_error` is null, `families_attempted` is exactly the
scope's families (Sporidiobolales: redPR+redHD; Pucciniales: rustHD; Wallemiales: wallMAT), and no family outside the
scope was searched. Taxonomy for routing came from NCBI E-utilities (`taxonomy_source` in every report), not the local
snapshot. Per genome wall time: median 84 s, 90th percentile 1,096 s, maximum 2,136 s (2 CPUs, 8 GB; the longest are
Pucciniales assemblies; 436 timing lines for 433 genomes). Four Pucciniales genomes had no report when their array
finished (most likely the 2 h `short` limit; I did not read the logs) and were run again in a final task.

## 2. Per-family table (tool output, whole run)

Columns are genome counts unless `loci`. `conf_*` counts loci; `best_*` counts called genomes by their best locus.

| family | in scope | called | loci | conf high / med / low | undetermined allele | partial | multi-locus genomes | both alleles in one genome |
|---|---|---|---|---|---|---|---|---|
| redHD | 251 | 246 | 252 | 252 / 0 / 0 | 246 | 0 | 5 | 0 |
| redPR | 251 | 244 | 248 | 248 / 0 / 0 | 0 | 0 | 4 | 4 |
| rustHD | 131 | 107 | 180 | 180 / 0 / 0 | 107 | 0 | 70 | 0 |
| wallMAT | 51 | 51 | 51 | 50 / 1 / 0 | 0 | 0 | 0 | 0 |

Allele or idiomorph labels (loci): redPR A1 = 59, A2 = 189; wallMAT v1 = 33, v2 = 18; redHD and rustHD: `undetermined`
for every locus, because their vocabulary is a pattern (alleles are not callable by homology; the labels in the records are
catalogue names).

Two things the table can mislead on:
- **`high` is not a discriminator for redHD and rustHD.** Every call is `high`. Both families have an all-`core_MAT`
  roster (two core genes, no flanking gene), so `classify_locus` files every call as `idiomorph_gene_only`
  (252 of 252 and 180 of 180) even with both genes present. This is the class's design, not a weak call.
- **`both alleles in one genome` for redPR (4 genomes)** is a count of genomes whose loci carry both A1 and A2:
  *R. pacifica* (GCA_049863195.1), *R. toruloides* (GCA_963855895.1, GCA_016808315.1) and a *Rhodosporidiobolus* genome
  (GCA_041531625.1). Each has two redPR loci on different contigs. I have not looked at whether these are mixed
  strains, diploids or contamination.

### By genus (called / genomes)

| family | genus | called / genomes | BFD taxids |
|---|---|---|---|
| redHD | *Rhodotorula* | 213 / 216 | 30 |
| redHD | *Sporobolomyces* | 19 / 21 | 10 |
| redHD | *Rhodosporidiobolus* | 13 / 13 | 10 |
| redHD | genus not named in NCBI | 1 / 1 | 1 |
| redPR | *Rhodotorula* | 211 / 216 | 30 |
| redPR | *Sporobolomyces* | 20 / 21 | 10 |
| redPR | *Rhodosporidiobolus* | 12 / 13 | 10 |
| rustHD | *Phakopsora* | 50 / 51 | 1 |
| rustHD | *Puccinia* | 42 / 51 | 25 |
| rustHD | *Melampsora* | 5 / 8 | 8 |
| rustHD | *Cronartium* | 5 / 6 | 6 |
| rustHD | *Austropuccinia* | 4 / 5 | 2 |
| rustHD | *Uromyces* | 1 / 4 | 2 |
| rustHD | *Hemileia* | 0 / 5 | 2 |
| rustHD | *Endocronartium* | 0 / 1 | 1 |
| wallMAT | *Wallemia* | 51 / 51 | 6 |

The 51 *Phakopsora* genomes are one taxid (170000); that row counts one species many times.

### Genomes not called

All uncalled genomes fail the same bar: the best cluster has fewer than 2 modelled genes (`modelled_gene_bar`). redHD 5,
redPR 7, rustHD 24. In none was an assembly gap flagged at the locus. In the 5 uncalled redHD genomes the withheld loci
include hits for both `HD1` and `HD2` (the union over that genome's withheld loci; I did not check that they sit in one
cluster); in 22 of the 24 rustHD genomes they include both `bE` and `bW`. So an uncalled genome here is often one where
the search saw the genes but the gene-model step did not model two of them; it is not evidence the locus is absent.

## 3. redHD: is a hit a divergent HD1/HD2 pair or two unrelated homeodomain genes?

Measured on every redHD locus (252 loci, 246 genomes), using the best-scoring `HD1` and `HD2` gene models on the locus's
contig:

| check | result |
|---|---|
| both `HD1` and `HD2` on the locus's contig | 252 / 252 |
| orientation | divergent 252, convergent 0, tandem 0 |
| intergenic bp | min 0 (one overlapping pair), 5th percentile 349, median 490, 95th percentile 870, max 1,295 |
| inside the records' 358-463 bp range | 54 |
| 1,000 bp or less | 247 |
| locus span | median 8,443 bp, 2,054 to 24,429 |
| best-hit identity to the curated proteins, median | `HD1` 65.5% (min 32.3%), `HD2` 82.9% (min 34.8%) |
| redHD and redPR locus on the same contig | 0 of 240 genomes that have both |

Judgement: a divergent pair under 1.3 kb in 252 of 252 loci is not what two unrelated homeodomain genes would give, so I
count these as HD loci. The 5 genomes with two redHD loci (3 loci in *R. pacifica*, 2 in each of four others) are not
resolved: they may be a second HD pair, a duplicated assembly region, or a second haplotype. The 54 within the record
range is a minority; the records span only five *Rhodotorula* and 358-463 bp, so I read the wider spread (up to
1,295 bp) as variation the records did not sample, not as a fault, but that is a reading.

## 4. Leave-own-record-out

Method: `scripts/run_asco_holdout.py run`, unchanged, on a manifest built with `build_asco_holdout_manifest.py`
extended with two radii the runner already accepts, `record` (drop only the record itself, the "same genome" test) and
`species`. For each curated record with a genome assembly in the library, `detect --taxid` runs on that record's own
genome with the radius's records excluded (`detect.holdout.records_to_withhold`). A task runs only if the radius leaves at
least one record of the same family. Pattern families (redHD, rustHD) are scored on whether a locus of the family is
found, because allele labels are not callable; enum families (redPR, wallMAT) on whether the record's own label is among
the called labels. 45 tasks ran. Not run: 35 for no reference left, 40 for no genome assembly.

| family | records (with assembly) | record | species | genus | family radius | order |
|---|---|---|---|---|---|---|
| redPR | 13 (5) | 5/5 | 5/5 | 5/5 | no reference | no reference |
| redHD | 5 (5) | 5/5 found | 5/5 found | no reference | no reference | no reference |
| rustHD | 4 (4) | 4/4 found | 4/4 found | 4/4 found | 4/4 found | no reference |
| wallMAT | 2 (2) | no reference for own version | same | no reference | no reference | no reference |

What the cells mean:
- **redHD, record radius** leaves the other four redHD records, all *Rhodotorula*, including the same species for two of
  the five. The locus was at the record's own coordinates in 5/5 (`found_at_truth`). **Genus radius leaves nothing**:
  all five records are *Rhodotorula*, and withholding the genus withholds all of them.
- **redPR genus radius** leaves 5 records (2 to 3 of the record's own allele) from *Sporobolomyces* genomes. 5/5 are
  correct; one is `medium` (R. toruloides NBRC 0880), the others `high`. Only 5 of 13 redPR records have a genome
  assembly in the library, so 8 cannot be held-out genomes.
- **rustHD family radius** is a real cross-family test: the Puccinia records are found with only the *Melampsora*
  record left, and the *Melampsora* record is found with only the three *Puccinia* records left. The locus was at the
  record's own span in the two records where coordinates are comparable; for the other two no contig is comparable.
  Scored as "locus found": 4/4 says nothing about alleles.
- **wallMAT:** with its own record withheld the locus is still found, but it carries the other version's label
  (the v1 genome is called v2 with `high`; the v2 genome is called v1 with `medium`). The runner records this as
  `no_reference_idiomorph`, a database gap, not a miss. Which version is which mating type is not known (placeholders).

The 12 record genomes are also in the whole-run counts above and were called there with the record's own label where
the family has labels: redHD 5/5 (locus), redPR 5/5 (A1 or A2 as in the record, `high`), rustHD 4/4 (locus), wallMAT 2/2
(v1 `high`, v2 `medium`).

## 5. The Rhodotorula question: where do these genera sit in BFD?

NCBI lineage from taxonkit (run under Slurm, local snapshot) for the BFD taxids, and a name search of the BFD species
names (`Rhodotorula`, `Rhodosporidium`, `Rhodosporidiobolus`, `Sporobolomyces`, `Sporidiobolus`):

| genus (NCBI) | order | genomes |
|---|---|---|
| *Rhodotorula* | Sporidiobolales | 216 |
| *Sporobolomyces* | Sporidiobolales | 21 |
| *Rhodosporidiobolus* | Sporidiobolales | 13 |
| (genus not named, family Sporidiobolaceae) | Sporidiobolales | 1 |
| *Cystobasidium* | Cystobasidiales | 16 |
| *Leucosporidium* | Leucosporidiales | 8 |
| *Erythrobasidium* | Erythrobasidiales | 4 |
| *Kondoa* | Agaricostilbales | 2 |
| *Cystobasidiopsis*, *Sterigmatomyces* | Agaricostilbales | 1 each |
| *Cystastrum* and one unnamed Cystobasidiaceae | Cystobasidiales | 1 each |
| *Curvibasidium*, *Symmetrospora* | no order in the NCBI lineage | 1 each |

**Genomes of *Rhodotorula*, *Rhodosporidiobolus*, *Sporobolomyces* (or the old names *Rhodosporidium*,
*Sporidiobolus*) outside Sporidiobolales: 0.** All 250 named ones, plus the one unnamed Sporidiobolaceae genome, are in
scope and were searched. Scope-only routing therefore leaves out no genome of these genera by name, in this BFD
snapshot.

What it does leave out, by design: the other Pucciniomycotina orders in the table, which carry names that used to be
*Rhodotorula* (*Cystobasidium*) or are close relatives. Whether they have a red-yeast P/R or HD locus is not known and
the family has no record from them. *Curvibasidium* and *Symmetrospora* have no order in NCBI, so no order-scoped family
reaches them. 147 BFD genomes have no taxonkit lineage; none of them is a Sporidiobolales, Pucciniales or Wallemiales by
the BFD CSV or by name.

## 6. Limits

- No truth for 421 of the 433 genomes. Counts of "called" can include a wrong locus; the arrangement check for redHD
  is the only independent check made here.
- redHD held-out recovery exists only at record and species radius, inside *Rhodotorula*. Whether it works in
  *Sporobolomyces* or *Rhodosporidiobolus* is not measured; 32 of 34 called there is the only evidence.
- rustHD alleles and redHD alleles are not callable; "5/5" and "4/4" mean the locus was found.
- rustHD loci per called genome: 70 of 107 genomes have two or more (180 loci). Dikaryotic assemblies carry two haplotypes,
  so this is expected, but I did not check haplotype against locus. The intergenic gaps are 536 to 1,797 bp (median
  1,739), 47 of 180 at 1,000 bp or less against 538-640 bp in the records; I did not break them down by genus.
- redPR labels: A1 59 against A2 189 loci. I have no explanation and no truth; I did not test whether A1 receptors are
  missed or relabelled.
- *Phakopsora* (51 genomes, one taxid) and *Puccinia* (51 genomes, 25 taxids) are the bulk of rustHD; the rate for
  Pucciniales is dominated by them.
- wallMAT 51/51 is 51 genomes of one genus, six taxids, with two records; a high call rate here does not test
  specificity or the v1/v2 call.
- Reports are written only when a genome finishes, so no partial report was used.
- Taxonomy for routing came from NCBI E-utilities during the run; the lineage table above came from the local taxonkit
  snapshot. I did not compare the two lineage by lineage; every report did attempt exactly its scope's families.

## 7. Open questions

Questions 2 and 4 are answered in sections 9 and 10 (follow-up).

1. A curated redHD record outside *Rhodotorula* (for example from *Sporobolomyces* or *Rhodosporidiobolus*) would make
   a genus-radius hold-out possible. Without one, redHD cannot be judged beyond the genus of its records.
2. Is the redPR A1:A2 skew (59:189) real, a reference imbalance, or a miss on A1?
3. The 5 redHD and 4 redPR multi-locus genomes: second haplotype, duplication or contamination.
4. Hemileia 0/5 and Uromyces 1/4: the bE and bW genes appear among withheld loci in 22 of 24 uncalled rust genomes.
   Does the modelled-gene bar (2 modelled genes) cost recall here, or are these true absences of the HD pair?
5. Should the `idiomorph_gene_only` class be renamed or special-cased for all-core pattern families, since it reads as
   a weak call and is the only class these families can produce?
6. The clean full-BFD run: use `scripts/family_performance_table.py` on that run directory, with a holdout manifest that
   includes the `record` and `species` radii (extended as here), and add the genomes with no taxonomy lineage to the
   routing log.

## 8. The per-family table tool

`scripts/family_performance_table.py REPORTS_DIR [--out T.tsv] [--md T.md] [--family Phylum:locus]
[--holdout-rows ROW.json|DIR ...] [--holdout-manifest M.tsv]`. It finds `detection_report.yaml` recursively (one
directory per genome, batches of those, or `<genome>/<mode>/`), counts a genome as in scope for a family when the
report's `families_attempted` lists it, and writes the table above plus held-out correct/tested by radius, with the
manifest's `no_assembly` and `no_reference_not_run`. For redHD (HD1/HD2) and rustHD (bE/bW) it adds the orientation and
gap of the pair. Nine tests, including one that pins its held-out verdicts to `run_asco_holdout.summarize`. Note that
`run_asco_holdout.py summarize` itself raises on the `record` and `species` radii (its order list has only genus,
family, order); the table tool handles them.

Reproduce: lists from `lin_ids.tsv` (taxonkit) and the BFD CSVs; one Slurm array of `python -m MATPredict detect
--genome G --taxid T --out-dir runs/<family>/<taxid>_<accession>` per chunk (2 CPUs, 8 GB, `short`); hold-out tasks
`run_asco_holdout.py run --manifest manifest_basidio_5radii.tsv --record R --radius X --library input_clean_genomes`;
then `family_performance_table.py runs/<family> --holdout-rows holdout/rows --holdout-manifest manifest_basidio_5radii.tsv`.

## 9. Follow-up: the redPR A1 59 : A2 189 skew (answers open question 2)

Data: `results/2026-10-18_basidio_scope_only/redPR_skew/` (per-locus call tables for the baseline and three
re-runs, `analysis_baseline.txt`, `balance_comparison.txt`); scripts in `followup_scripts/`. No change to `src/` or `db/`.

**Short answer.** It is not a reference imbalance and not a small-margin artefact; every call is decided by a wide
margin. The skew is mostly one species: *R. mucilaginosa* has 125 of the 251 genomes, from one taxid in 26
BioProjects, and 117 of its 123 loci are A2. Without it the loci are A1 53 : A2 72. Whether the excess of A2 in
*R. mucilaginosa* is biology (an A2-heavy sampled population or a clonal species) or sampling I cannot say: there is no
independent truth for it.

**References (measured).** The database holds 13 redPR records: A1 6, A2 7. *Rhodotorula* 8 (A1 4, A2 4: *toruloides*
A1 and A2, *glutinis* A1, *mucilaginosa* A1, *kratochvilovae* A1 and A2, *sphaerocarpa* A2, *graminis* A2),
*Sporobolomyces* 5 (A1 2, A2 3: *pararoseus* A1 and A2, *salmonicolor* A1 and A2, *longiusculus* A2). *Rhodosporidiobolus*:
none. Only 5 have a genome assembly in the library (A1 3, A2 2; all *Rhodotorula*). So the reference set is balanced
by allele (6 : 7) and, within *Rhodotorula*, exactly (4 : 4).

**Calls by genus (loci).**

| genus | genomes | A1 | A2 |
|---|---|---|---|
| *Rhodotorula* | 216 | 42 | 172 |
| *Sporobolomyces* | 21 | 14 | 6 |
| *Rhodosporidiobolus* | 13 | 3 | 10 |
| not named | 1 | 0 | 1 |

*Rhodotorula* without *R. mucilaginosa*: A1 36 : A2 55. By species (44 species with a call): A1 only 10, A2 only 22,
both alleles among that species' genomes 12. By taxid (48): A1 only 12, A2 only 24, both 12. Counting each
(species, BioProject) as one unit (131 units): A2 only 77, A1 only 43, both 11. *R. mucilaginosa*: 125 genomes, 26
BioProjects, A2 117, A1 6, no call 2. The largest BioProject (PRJNA723004, 51 genomes) is 49 A2, 2 A1. The 6 A1
*R. mucilaginosa* genomes are 99.7 to 100% identical to the *R. mucilaginosa* A1 record (JY1105); the A2 calls in that
species are about 74% identical to the nearest A2 record (a *R. sphaerocarpa* record), because there is no
*R. mucilaginosa* A2 record.

**Margins (measured, 248 loci).** For each locus I took the best `STE3a1` and `STE3a2` hit.

| | winner | loser |
|---|---|---|
| A1 calls (59): bitscore, median (min to max) | 435 (144 to 633) | 89 (0 to 122) |
| A2 calls (189): bitscore, median (min to max) | 299 (109 to 479) | 92 (30 to 117) |
| A1 calls: identity, median | 99 (min 63) | 34 (max 40) |
| A2 calls: identity, median | 74 (min 65) | 35 (max 51) |

Bitscore margin (winner minus loser): minimum 66, 0 calls under 50 bits, 3 under 100, 41 under 200. No call has a
winning identity under 60%. The losing allele's hit is a short cross-hit at about 35% identity; A1 and A2 receptors
are about 35% identical to each other, so the two alleles do not blur. No call is decided by a small margin.

**Does the call just follow the nearest record?** The label is by construction the allele of the best-hit
reference, so the real question is whether one allele wins only because it has the only near record. 133 of the 248
calls have the *R. sphaerocarpa* A2 record as nearest (all the *R. mucilaginosa* A2 calls, 74% identity), and 233 of 248
have a nearest record of their own genus (there is no *Rhodosporidiobolus* record). The re-runs below test this directly.

**Balance test (measured).** `detect --exclude-records` on the same db, no copy:

| run | references kept | genomes | result vs baseline |
|---|---|---|---|
| E1 | one A1 and one A2 per genus (*toruloides* cbs-14 and nbrc-0880, *salmonicolor* cbs-483, *longiusculus* cbs-9654; 4 of 13) | all 251 | A1 59 : A2 189, **identical call in 251 of 251 genomes**; 0 lost, 0 gained, 0 flips |
| E3 | A2 scarce: all 6 A1, 2 A2 | 80 (random, seed 1) | A1 20 : A2 56, identical 80 of 80 |
| E4 | A1 scarce: all 7 A2, 2 A1 | the same 80 | A1 20 : A2 56, identical 80 of 80 |

Cutting either allele's references to 2 changes nothing, so the skew is not driven by how many references each allele has.
Two limits: every variant keeps references of both alleles, so it does not test an allele wholly absent from the
references; and a detector that systematically misses a divergent A1 lineage cannot be excluded here. Only 7 genomes have
no call, and their withheld loci carry hits for both `STE3a1` and `STE3a2`, so I cannot say which allele each would have
been.

**Genomes carrying both alleles (4, as in section 2).** *R. pacifica* (GCA_049863195.1), *R. toruloides*
(GCA_963855895.1, GCA_016808315.1), *Rhodosporidiobolus fluvialis* (GCA_041531625.1). Each has an A2 and an A1 locus on
different contigs; in the two *R. toruloides* genomes the loci are 98 to 100% identical to the *R. toruloides* A1
(cbs-14) and A2 (nbrc-0880) records. What this means is not determined by this data. Possible causes are a diploid or
heterokaryotic strain, a hybrid, or a mixed or contaminated assembly; an A1/A2 pair in one genome is what a sexual
diploid would show, but so would contamination. The A1 locus of GCA_016808315.1 starts at position 1,000 of
contig JABGON010000578.1 and ends at 3,065, which fits a small fragment, but I did not check contig lengths.

**Independent truth.** None for the 251 genomes except the five record genomes (all five were called with the record's
own allele in the baseline, and recovered at record, species and genus radius in section 4). The allele names are the
curated catalogue's (Coelho et al. 2008 PMID 18408057; 2010 PMID 20700437; 2011 PMID 21880139 for the *Sporobolomyces*
deposits; the *Rhodotorula* records come from the not yet peer-reviewed Liu et al. bioRxiv 10.1101/2025.09.11.675505). I did
not find, and do not cite, allele frequencies for wild populations of *R. toruloides*, *Rhodosporidiobolus* or
*Sporobolomyces*, so there is no expected A1 : A2 ratio to compare with.

**Judgement.** The skew reflects the BFD sample (one species supplies half the genomes and is 95% A2), not the
detector: the calls are wide-margin, the 4-reference run reproduces them exactly, and species-level counts are far more
even (A1 only 10, A2 only 22, both 12). For the clean full-BFD run, report allele counts per species or per BioProject as
well as per genome; the per-genome 59 : 189 should not be read as an allele frequency.

## 10. Follow-up: why rustHD is uncalled in *Hemileia*, *Endocronartium*, *Uromyces* and *Melampsora* (answers open question 4)

Data: `results/2026-10-18_basidio_scope_only/rustHD_uncalled/` (`summary.tsv`, per-genome JSON in `rc_json.tar.gz`,
`relax_*.tsv`, `hd_pairs.txt`, `n50_bins.txt`). Analysis only: nothing in `src/` or `db/` was changed.

**Short answer.** Two causes, and they separate by genus. (a) Assembly fragmentation: 20 of the 24 uncalled genomes
have a BFD N50 under 10 kb, including all *Uromyces* and *Endocronartium* and most *Melampsora* and *Puccinia*.
(b) *Hemileia*: the locus is present, but the *Puccinia* and *Melampsora* references align only the homeodomain, so the
modelled-gene bar of 2 is not met even in the two chromosome-scale assemblies. Lowering the bar to 1 recovers 7 of the 24
uncalled genomes (all 5 *Hemileia*, 1 *Endocronartium*, 1 *Melampsora*) and adds one probable false positive.

**Contiguity (measured).** Calls by BFD N50 over the 131 genomes: under 10 kb 14 of 34, 10 to 100 kb 5 of 5, 0.1 to 1 Mb
13 of 15, 1 Mb or more 75 of 77. So 93 of 97 genomes with N50 of 10 kb or more are called (96%) against 14 of 34 (41%)
below. Median N50 of the 42 genomes I re-examined (computed on the cleaned FASTA): uncalled 4.9 kb (24 genomes), called
4.8 Mb (18 genomes, a stratified sample). By genus (median N50, contigs): *Hemileia* 11 kb, 33,176 contigs (5 genomes; two
are 794 and 31 contigs with N50 36 and 35 Mb); *Endocronartium* 4.2 kb, 10,136; *Uromyces* 2.5 kb, 4,387 (4 genomes, 2 to
174 Mb; the 174-Mb one is called); *Melampsora* 7.4 kb, 7,507; the 9 uncalled *Puccinia* genomes are 4 to 77 Mb (three
*P. striiformis* isolates at 41 to 47 Mb) against 83 to 173 Mb for the called controls and 111 Mb for the complete
*P. striiformis*. The same effect, much weaker, holds for the other families (called with N50 under 10 kb: redHD 1 of 4,
redPR 3 of 4, against 18 of 18 and 17 of 18 at 10 to 100 kb).

**Independent check (42 genomes: the 24 uncalled, the 18 *Hemileia*, *Endocronartium*, *Uromyces* and *Melampsora*
genomes, and 12 called controls).** Per genome: tblastn at 1e-3 and miniprot of the 8 rustHD record proteins, and Pfam
Homeodomain and Homeobox_KN on six-frame ORFs around the hits.

| class | meaning | target-genus genomes (18) | other uncalled (12) | called controls (12) |
|---|---|---|---|---|
| A | both genes modelled in full (miniprot coverage 0.8 or more) on one contig within 20 kb | 4 | 0 | 6 |
| B | a bE/bW pair on one contig, but a fragment (partial model or short contig) | 3 | 0 | 6 |
| C | both genes hit, but not on one contig | 6 | 3 | 0 |
| D | one gene only | 3 | 5 | 0 |
| E | no hit | 2 | 4 | 0 |

Of the 24 uncalled genomes, 23 are class C, D or E (what a split or incomplete assembly gives) and one, *Endocronartium*
(class B), is a divergent bE/bW pair about 700 bp apart on a 2.8-kb contig. Per genus, uncalled: *Hemileia* 5 (all C),
*Uromyces* 3 (D, D, E; genomes of 2 to 16 Mb, N50 about 2.5 kb), *Melampsora* 3 (E, C, D; N50 4 to 5 kb), *Puccinia* 9
(C 3, D 3, E 3), *Cronartium* 1 and *Austropuccinia* 1 (D), *Phakopsora* 1 (E: a 0.9-Mb, 2-contig assembly). Frameshifts
were not the problem: 2 of the miniprot models in the 42 genomes carry one.

***Hemileia* is different.** Two *H. vastatrix* assemblies are contiguous (GCA_030280995.1, 794 contigs, N50 36 Mb;
GCA_054948715.1, 31 contigs, N50 35 Mb) and still uncalled. In all five *Hemileia* genomes the references align only the
homeodomain: tblastn gives short hits of about 42 to 45 bits covering 16 to 17% of the bE query, miniprot gives bE models
of 13 to 35% of the query at 28 to 45% identity in four genomes and none in the fifth, and bW models of 32 to 81% at 24 to 30%
identity (in the two chromosome-scale genomes the bW "model" spans 32 kb, so it is not a gene). A genome-wide Pfam scan of the
two contiguous genomes finds, in each, a divergent HD pair on one contig at the tblastn site: JASXTE010000003.1
(Homeodomain on `-` at 13,705,138, Homeobox_KN on `+` at 13,707,238, 1,957 bp between the domains) and CM142961.1
(Homeobox_KN on `-` at 37,098,429, Homeodomain on `+` at 37,100,502, 1,957 bp). The same spacing in two independent
assemblies, in a divergent arrangement, is what a bE/bW-type pair looks like. The first genome has one more divergent
pair on a small contig (JASXTE010000732.1, 1,050 bp), perhaps a fragment of the same locus. Judgement: *Hemileia*
carries the locus and the pipeline misses it only because the HD proteins are too different from the *Puccinia* and
*Melampsora* references to be modelled. That these two Pfam domains are the *Hemileia* bE and bW cannot be shown without a
*Hemileia* record. The genome-wide scan is insensitive: on *P. striiformis*, with a known locus, it recovered only one of the
two domains, so a missing pair in it means little.

**Relaxing the bar (analysis only: `run_pipeline`'s own `min_polished_genes` set to 1 through a wrapper, no code change).**
All 24 uncalled genomes and 18 called ones were re-run.

| | result |
|---|---|
| uncalled genomes that gain a call | 7 of 24: *Hemileia* 5 of 5, *Endocronartium* 1, *Melampsora* 1 (so rustHD 107 to 114 of 131, 82% to 87%) |
| called genomes that gain extra loci | 0 of 18 re-run (the other 89 called genomes were not re-run) |
| new loci | 7, all `idiomorph_gene_only`, 1 polished gene, confidence medium 6, high 1 |
| arrangement of the 7 | divergent 6 (intergenic 696 to 2,265 bp), tandem 1 (gap 0) |
| the *Hemileia* loci | 5, divergent, 2,045 to 2,265 bp; the two chromosome-scale ones coincide with the Pfam pairs above |

The tandem one is *Melampsora* GCA_002157015.1 (contig KZ132129.1, 1,420 bp, `bE` and `bW` overlapping, confidence
`high`). It is the genome where tblastn and miniprot find a single full-length gene (class D), so `bE` and `bW` are
hitting the same gene: a false positive. An arrangement filter of the kind used for redHD (both genes on one contig,
non-overlapping, divergent, gap about 0.5 to 3 kb) keeps 6 of 7 and rejects it. The six accepted include all five *Hemileia*
genomes, one species. The 17 genomes still uncalled at bar 1 have no pair (classes C, D, E); the bar does not help them.

**What changes the clean-run counts.** (i) rustHD called would be 107 of 131 as it stands, 114 of 131 with the bar at 1
(7 genomes), or 113 with an arrangement filter. (ii) Any Pucciniales count should be stratified by N50: below 10 kb a no-call
is mostly "assembly", not "absent". (iii) A *Hemileia* record would make the locus callable without lowering the bar.
(iv) redPR and redHD are unaffected: the redPR tests above changed no call. Not measured: extra loci a relaxed bar adds in the
89 called genomes not re-run, and a bar of 0.

**Open question added.** Is the bar of 2 modelled genes worth keeping for all-core pattern families? For rustHD it costs 7 of
131 genomes (5.3%) in this run, at about one false positive in 7; that is a curator decision, not a fix I would apply unasked.
