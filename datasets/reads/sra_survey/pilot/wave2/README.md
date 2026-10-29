# Wave 2 (2026-10-28): reference-tier A, A1 and B species not in wave 1, 82 planned, 55 typed

`wave_driver.py` planned the tier A, A1 and B species that wave 1 had skipped (85 candidates, 82 planned after the yeast naming support and borrowing described below), built a panel for each from the BFD assemblies, picked up to 50 runs per species (natural-isolate studies, 0.6 to 12 Gb, 4 per study at most; up to 20 at the test stage) and typed them on the HPCC (`exfab` and `stajichlab`, 1 CPU and 6 GB per task). Every run also wrote a k=21 read sketch (`SKETCH=1`) for clone detection. *Botrytis cinerea* was then typed on all 333 eligible samples (see below).

Files: `gates.tsv` (test stage), `gates_final.tsv` (final table), `results_all_final.tsv` (one row per run), `species.tsv` (the species and tiers), `clone_pairs/<slug>.pairs.tsv` (pairwise sketch identity for the 32 species with 10 or more sketched runs), `botrytis_B05.10_MAT_locus.md`.

## Result
- Runs: 1,468; 1,467 finished, 1 unavailable at NCBI; 1,911 GB streamed (and discarded) from 469 studies. Calls: MAT1-1 357, MAT1-2 331, a 153, alpha 151, A 51, B 5, both 278, none 82, low_depth 59.
- Gate (final data): PASS 13, PROVISIONAL 16, PROVISIONAL_HOMOTHALLIC 7, HOLD 49 (the pages say confirmed / provisional / withheld). At the test stage the counts were 11 / 18 / 7 / 49; the full runs moved *Neurospora cerealis* and *Bipolaris oryzae* from provisional to confirmed (*Botrytis cinerea* was confirmed at both stages).
- Yeast naming: the panels name the idiomorphs a / alpha (and P / M, A / B), the names are carried through the summary, the gate and the pages. In a diploid yeast with one copy of each idiomorph a `both` call is the expected answer, so those species get the homothallic gate status (both expected) and the page says so; the single-call shares of such species are not allele frequencies.
- Borrowing: a species that had one idiomorph in the assemblies was given loci from congeners (`panel_from_assemblies.py --borrow "Genus "`): borrowed loci join only the species' own flank lineages, so a congener cannot bury the species' own loci.

## *Botrytis cinerea* (all strains)
333 samples of 11 studies (survey: 449 runs, 353 samples, 17 studies; 328 pooled strains after filters; mutant, cross and UV-mutagenesis studies excluded; PRJNA1295911 supplies 276 of them). Calls: MAT1-1 151, MAT1-2 157, both 24, low_depth 1; against `detect` 45 of 46 agree (SRR11527078, strain Rf2, is the one disagreement). Clones: 333 runs fall into 197 clones (average linkage at 0.998; the number is 172 at 0.997). By study the MAT1-1 share of single calls is 48% (PRJNA624742, 32 runs) and 51% (PRJNA1295911, 276 runs), so the pooled 49% is not one study's value; the reports show both views. The B05.10 runs are discussed in `botrytis_B05.10_MAT_locus.md` (the reads lack the MAT1-1 core).

## Held species (49), by reason
- **Panel build failed (16):** the panel could not be built (mostly yeasts where the assemblies show one idiomorph and borrowing from congeners did not supply the other; the reason was not examined species by species): *Candida dubliniensis*, *Candida viswanathii*, *Debaryomyces hansenii*, *Kluyveromyces lactis*, *Lachancea fermentati*, *Lachancea thermotolerans*, *Lodderomyces orthopsilosis*, *Schizosaccharomyces japonicus*, *Zygosaccharomyces bailii*, *Zygosaccharomyces favi*, *Zygosaccharomyces mellis*, *Zygosaccharomyces parabailii*, *Zygosaccharomyces pseudobailii*, *Zygosaccharomyces rouxii*, *Zygosaccharomyces sapae*, *Zygosaccharomyces siamensis*.
- **No idiomorph-specific k-mers in the largest lineage (3, recorded as `panel_L1_error` in `gates_final.tsv`; the driver now writes `panel_no_unique_kmers`):** *Ceratocystis fimbriata*, *Schizosaccharomyces pombe* and *Sclerotinia sclerotiorum* each have a lineage 1 (the one with most loci) in which one idiomorph has no k-mer absent from the other (`qc_panel.json`). The cause was confirmed for *S. pombe*: the lineage 1 representative locus (36 kb, GCA_006386395.1) is the same locus for P and M and its two cores are 14,111 bp of identical sequence, so the region holds both idiomorphs' sequence (the *mat1* cassette with the silent *mat2-P* and *mat3-M* cassettes) and every strain would read as both. *S. sclerotiorum* is homothallic with both MAT genes in one locus (`panel_homothallic.py` was written for it in the pilot); for *C. fimbriata* (unidirectional switching) the cause was not checked. Typing the presence of an idiomorph cannot work for these; the hold is the right answer, and the work for *S. pombe* would be a panel cut from the *mat1* cassette alone, which still reads as both because the silent cassettes are in every cell.
- **Too few loci (5):** *Aspergillus westerdijkiae*, *Candida ipomoeae*, *Candidozyma pseudohaemuli*, *Fusarium sambucinum*, *Morchella sextelata*.
- **Too few runs (6):** fewer than 4 runs after the run filters (`min_runs` 4 in `wave_driver.py`): *Aspergillus minisclerotigenes*, *Candida boidinii*, *Fusarium virguliforme*, *Ogataea angusta*, *Ogataea parapolymorpha*, *Pneumocystis jirovecii*.
- **Typing did not pass (19):**

| Species | Tier | Gate numbers |
|---|---|---|
| *Aspergillus hortae* | B | runs 5 done 5 typed 4 both 0 none/low 1 truth_typed 3 agree 3 |
| *Aspergillus sojae* | B | runs 11 done 11 typed 7 both 0 none/low 4 truth_typed 3 agree 3 |
| *Aspergillus tamarii* | B | runs 11 done 11 typed 8 both 0 none/low 3 truth_typed 3 agree 3 |
| *Bipolaris gigantea* | B | runs 4 done 4 typed 3 both 1 none/low 0 truth_typed 3 agree 2 |
| *Bipolaris maydis* | A | runs 10 done 10 typed 8 both 0 none/low 2 truth_typed 2 agree 2 |
| *Candida albicans* | A | runs 20 done 20 typed 7 both 8 none/low 5 truth_typed 5 agree 2 |
| *Candidozyma vulturna* | B | runs 11 done 11 typed 9 both 0 none/low 2 truth_typed 0 agree 0 |
| *Fusarium boothii* | B | runs 8 done 8 typed 0 both 0 none/low 8 truth_typed 0 agree 0 |
| *Fusarium keratoplasticum* | B | runs 7 done 7 typed 3 both 0 none/low 4 truth_typed 0 agree 0 |
| *Fusarium meridionale* | B | runs 5 done 5 typed 0 both 0 none/low 5 truth_typed 0 agree 0 |
| *Geotrichum candidum* | A | runs 12 done 12 typed 6 both 1 none/low 5 truth_typed 3 agree 3 |
| *Lodderomyces elongisporus* | B | runs 9 done 9 typed 0 both 0 none/low 9 truth_typed 0 agree 0 |
| *Morchella eximia* | B | runs 6 done 6 typed 2 both 2 none/low 2 truth_typed 1 agree 1 |
| *Ogataea haglerorum* | B | runs 5 done 5 typed 0 both 0 none/low 5 truth_typed 0 agree 0 |
| *Saccharomyces eubayanus* | B | runs 20 done 20 typed 2 both 0 none/low 18 truth_typed 0 agree 0 |
| *Saccharomyces mikatae* | B | runs 7 done 7 typed 0 both 0 none/low 7 truth_typed 0 agree 0 |
| *Saccharomyces paradoxus* | B | runs 20 done 20 typed 0 both 0 none/low 20 truth_typed 0 agree 0 |
| *Saccharomyces pastorianus* | B | runs 20 done 20 typed 3 both 0 none/low 17 truth_typed 0 agree 0 |
| *Saccharomyces uvarum* | B | runs 20 done 20 typed 14 both 0 none/low 6 truth_typed 0 agree 0 |

The *Saccharomyces* sensu stricto species (*S. eubayanus*, *S. paradoxus*, *S. uvarum*, *S. pastorianus*, *S. mikatae*) have a median flank depth below 2x, so most runs got no call; the cause was not examined (hybrid strains and a panel from few loci are possibilities). *Candida albicans* is called `both` in 40% of runs and 2 of 5 comparisons agree.

## Clone groups
The 32 species with 10 or more sketched runs have clone groups on their reports (`clone_pairs/`, rendered by `scripts/build_reads_site.py`: average linkage at 0.998, counts per clone only where the clone number is stable between 0.998 and 0.997). Most species are nearly all singletons (*Candidozyma auris* 50 runs, 50 clones; *Zymoseptoria tritici* 49 clones). Clonal species: *Komagataella phaffii* (50 runs, 10 clones, one of 41), *Lodderomyces parapsilosis* (2 clones, one of 49), *Aspergillus nidulans* (3 clones, one of 48), *Yarrowia lipolytica* (one of 25). The numbers printed by the cluster job on the HPCC came from an older copy of `clone_clusters.py` (219 clones for *B. cinerea*) and are not used; the pages use the repository version on the saved pairs.

## Not checked
Whether the single-call shares of yeast species with few strains are meaningful; why the *Saccharomyces* sensu stricto and *C. albicans* panels fail (the work list for a curated a / alpha panel); why *C. fimbriata* has no unique k-mers in its largest lineage; clone groups for species with fewer than 10 sketched runs. Counts are of sampled strains, not population frequencies.
