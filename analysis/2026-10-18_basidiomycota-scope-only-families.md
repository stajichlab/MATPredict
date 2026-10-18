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
