# B mating-type locus clustering in Agaricomycetes: what the receptor-array data can and cannot say about calling
Status: open (assessment only; no pipeline change, nothing merged, no PR)

## Question
The curator's point: pheromone receptors (STE3-like) cluster in B loci, so B-locus clustering should be assessed and used in calling. Here: (1) what curated B loci look like and which receptor array holds them; (2) which array properties distinguish the B array from other receptor arrays; (3) whether a "cassette" (receptor plus two or more precursor ORFs within a few kb) is a usable unit over the 1,287 quality-passing genomes; (4) options.

## What is known (literature, as given; not re-derived)
- Tetrapolar Agaricomycetes: the B locus holds pheromone receptors (STE3-like) and pheromone precursors (Casselton and Olesnicky 1998, MMBR 62:55). *S. commune* has Balpha and Bbeta (separated by up to about 3.5 map units), each with one receptor and several (three) precursors. *C. cinerea* has one locus of about 17 kb with three receptor-plus-two-pheromone cassettes. *U. maydis* a1/a2: one pheromone and one receptor.
- Nomenclature trap for this repository: the Ustilaginales `bLocus` records in `db/` are bE/bW (HD genes), not receptors; the Ustilago receptor locus is `aLocus`, and Sporidiobolales receptor records are `redPR`. Outgroup B-type contrasts below therefore mean the receptor-bearing `a`/P-R loci, and none was scanned here (the array study covers Agaricomycetes only).
- In the pipeline, `PR` (Basidiomycota:PR) is the generic receptor family; `locus_merge.py` merges PR, Balpha and Bbeta calls at one span into one B call (group "B").

## Data and method
Everything is reused; no scan was rerun on the full set. Study data: `results/2026-10-06_agaricomycetes_pr_arrays/` (`inventory/ste3_loci_table.tsv.gz`: 12,518 loci in 9,436 arrays of the 1,287 qpass genomes, 705 species; 50-kb single-linkage arrays), per-genome strict-CAAX ORF positions (`T`, `C[VI][IV][AVMG]` before a stop, Met 20-130 codons upstream, six-frame, no annotation) and tblastn precursor-homology hits (`Hx`, E <= 1, `pheromones_curated.faa`) from `scan_out_full.tar.xz` (the existing `scan_genome.py`/`array_scan.py` output; I did not reimplement the scan). Labels: `labelled_loci.tsv` from `origin/pr-receptor-ml-investigation` (copied here, with review `origin/pr-receptor-ml-review`). Curated records: `db/Basidiomycota/{Agaricales,Polyporales,Russulales}`.
- New HPCC step (Slurm job 29530344): tblastn of `pheromones_curated.faa` with the genome's own-species records removed, on the 6 genomes that carry a curated B record (`heldout_hx.py`, `run_heldout.sh`, results in `hx_heldout/`). Needed because `pheromones_curated.faa` holds the *S. commune* and *C. cinerea* precursors, so the main scan's homology count in those two genomes is a self-hit.
- Cassette definitions (positions as distinct ORFs; homology hits within 300 bp merged; "within 5 kb" is from either end of a receptor locus; `b_cluster.py`):
  A = at least 2 distinct strict-CAAX ORFs; B = at least 2 precursor candidates (CAAX or homology) and at least one that is both; C = at least 2 CAAX ORFs that each carry homology.
- Chance: the receptor locus placed uniformly at random on its own contig (40 draws), counting only CAAX ORFs that are not within 10 kb of any receptor locus. Same logic as the array study; it does not model local gene density.
- Code and outputs: `results/2026-10-06_b_locus_clustering/` (`b_curated.py`, `b_cluster.py`, `b_summarise.py`, `b_extra.py`; `curated_agaricomycete_B_records.tsv`, `b_vs_other_arrays_6genomes.tsv`, `b_loci_cassette.tsv.gz`, `b_arrays_features.tsv.gz`, `cassette_by_order.tsv`, `summary_output.txt`, `extra_output.txt`).

## Circularity ledger (read before the results)
| Item | Why it is circular, or not |
|---|---|
| Strict-CAAX ORF counts and the CAAX flag | The pipeline admits receptor + CAAX as a PR call (944 of 1,039 calls are `unverified` on that basis). Any "B array has CAAX ORFs" statement over a set of CAAX-selected calls is true by construction. |
| Labels of *Trametes*, *Grifola*, *Heterobasidion*, *Russula* (grade C) | The `labelled_loci.tsv` evidence column reads "genome-derived record, chosen from CAAX positional evidence (tier 2)". Their records were built around CAAX-bearing receptor clusters. CAAX counts in these 4 arrays carry no evidential weight. |
| Independent labels | Only *S. commune* (bar3/bbr2, transformation tests, PMID 7489716) and *C. cinerea* (B43, mapped, PMID 9539426) are genetically mapped. n = 2 genomes. Their precursor sets come from the literature, but those precursors are in `pheromones_curated.faa`, so unfiltered homology in them is a self-hit (see held-out numbers). |
| Flank genes | The flank set (STE20, PAN6 ...) was chosen knowing MAT linkage in Tremellales and Sporidiobolales; no Agaricomycete B-locus flank gene is established (the curated flank set in the db for Agaricomycetes, MIP1 and beta_fg, flanks the A locus). The curated Agaricomycete B records list zero flank genes. STE20 is used below only as a semi-independent co-location check, not as a label. |
| Called status (`pipeline_call`) | Derived from CAAX admission; not a label. Not used to define truth anywhere below. |

## Results

### 1. Curated B loci (`curated_agaricomycete_B_records.tsv`)
| Record | Species (order) | Span kb | Receptors | Precursors | Layout | Flank genes |
|---|---|---|---|---|---|---|
| Balpha 5334 | *S. commune* (Agaricales) | 6.9 | 1 | 2 (literature: 3) | one sublocus | 0 |
| Bbeta 5334 | *S. commune* | 17.2 | 1 | 7 | one sublocus | 0 |
| PR_B43 5346 | *C. cinerea* | 20.7 | 4 listed (literature: 3 groups) | 4 listed (literature: 2 per group) | receptor gaps median 2.3 kb | 0 |
| PR_B1 5325 | *Trametes versicolor* (Polyporales) | 24.5 | 3 | 9 | median receptor gap 6.1 kb | 0 |
| PR_B1 5627 | *Grifola frondosa* | 34.9 | 5 | 3 | gap 4.9 kb | 0 |
| PR_B1 2830151 | *Russula nobilis* (Russulales) | 14.6 | 3 | 1 | gap 2.9 kb | 0 |
| PR_B1 984962 | *Heterobasidion irregulare* | 23.5 | 5 | 3 | gap 3.8 kb | 0 |
- Per-genome curated B record: 6 genomes, 7 records (*S. commune* has two). Curated spans 7-35 kb; the combined *S. commune* Bbeta+Balpha region is 31.9 kb, with 7.7 kb between the last Bbeta gene and the first Balpha gene (physical separation of the two subloci is small compared with the map distance). *C. cinerea* 20.7 kb against the known 17 kb: the record lists 4 receptor genes where the literature says 3 groups (two adjacent receptor annotations, 1819643-1821277 and 1817324-1819437, merge into one inventory locus and may be one gene split in two; not resolved here).
- Counts are curated and partial ("completeness: partial" in every record), so they are descriptive. Curated precursors per receptor range 0.3 (*Russula*) to 7 (*S. commune* Bbeta) (Russula 1 precursor for 3 receptors; Trametes 9 for 3). The textbook "receptor plus two pheromones" is the *C. cinerea* pattern, not a rule in the other orders: 3 of the 7 records have fewer than 1 precursor per receptor.
- No record lists a flank gene; there is nothing curated to compare flank presence against.

### 2. Which array holds the B receptors, and how many paralogs sit with them (`b_vs_other_arrays_6genomes.tsv`)
| Genome | B array loci (of which curated B) | Span kb | Other STE3 loci in the same array | Curated receptors vs inventory loci |
|---|---|---|---|---|
| *S. commune* (independent) | 4 (2: Bbeta bbr2, Balpha bar3) | 88.5 | 2 (brl1 at 190 kb and brl3 at 268-279 kb flank the whole B region; excluded from both records) | 2 receptors, 2 loci |
| *C. cinerea* (independent) | 4 (3) | 77.4 | 1 (an unlabelled locus of 65% identity to its query, about 45 kb from the B receptors) | 4 listed, 3 loci |
| *T. versicolor* | 5 (2) | 54.9 | 3 | 3 receptors, 2 loci in the table |
| *G. frondosa* | 3 (3) | 30.0 | 0 | 5 receptors, 3 loci |
| *H. irregulare* | 3 (3) | 25.4 | 0 | 5 receptors, 3 loci (the first locus, 19 kb, merges 4 receptors) |
| *R. nobilis* | 1 (1) | 14.1 | 0 | 3 receptors, 1 locus (14 kb) |
- In 6 of 6 genomes all labelled B receptors sit in one array (the B array is unique per genome). In 3 of 6 (*Grifola*, *Heterobasidion*, *Russula*) the array holds only B receptors; in the other 3 it also holds 1-3 receptor loci outside the curated record (*S. commune* 2, *C. cinerea* 1, *Trametes* 3).
- Counting loci undercounts receptors: the miniprot merge joins tandem same-strand receptors whose alignments overlap, so the 16 curated receptors of *Trametes*, *Grifola*, *Russula* and *Heterobasidion* appear as 9 loci (22 curated receptors in all 6 genomes appear as 14 loci). A "receptor count in the array" from this table is a lower bound; *R. nobilis*'s three tandem receptors appear as one singleton locus, so "array size >= 2" would miss that B locus entirely.
- Span: the B arrays span 14-88 kb (median 54). The 17-kb *C. cinerea* locus and the 32-kb *S. commune* B region are each smaller than the array containing them (77 and 88 kb), because the array rule (gap <= 50 kb, single linkage) extends over adjacent paralogs. Array span is not B-locus span.
- Sublocus: in *S. commune* the Bbeta and Balpha receptors are two separate inventory loci 26 kb apart (gene-to-gene) in one array, each with its own CAAX/homology cassette (Bbeta locus: 2 CAAX within 5 kb, Balpha locus: 3). The array structure shows two cassettes but cannot label them Balpha and Bbeta; that needs a curated reference (the pipeline's Balpha/Bbeta families do this).

### 3. Which array properties distinguish the B array (descriptive; n = 6 B arrays, 25 other arrays in the same 6 genomes)
| Property | B arrays (n=6) | Other arrays in those genomes (n=25) |
|---|---|---|
| Loci per array | 1, 3, 3, 4, 4, 5 | 24 singletons, 1 of size 3 (*Russula*) |
| Span kb | 14-88 | 1.2-58 (the size-3 array 58) |
| Strict-CAAX ORFs within 10 kb of the array | 1, 4, 4, 4, 7, 9 (median 4) | 0 in 24, 1 in 1 |
| Cassette A present (>= 2 CAAX within 5 kb of a receptor) | 4 of 6 | 0 of 25 |
| Held-out precursor homology clusters within 10 kb | 0 to 2 (*S. commune* 0, *C. cinerea* 1, *Trametes* 2, *Heterobasidion* 1, *Grifola* 0, *Russula* 0) | not tabulated; the main-scan counts were self-hits in 2 genomes (*S. commune* 11, *C. cinerea* 6) |
| Flank (STE20 top hit within 250 kb; semi-independent, see ledger) | 5 of 6 (9.5-81 kb; *Trametes* none) | 0 of 25 |
| Contains a pipeline PR call (this table; *S. commune* is reported through the Balpha/Bbeta families, not checked here) | 5 of 6 | 0 of 25 (2 inside a withheld cluster) |
- Exact test of cassette A, B arrays against other arrays: 4/6 vs 0/25, Fisher p = 0.0005; on the 2 independently labelled genomes: 2/2 vs 0/9, p = 0.018 (the smallest p achievable). Loci and arrays inside one genome and genomes of one species are not independent, and the 4 grade-C arrays are circular for every CAAX-derived column. The honest content: in the only two independently mapped genomes the B array is the only array with a cassette, and it has two (*S. commune*) or one dominant (*C. cinerea*) cassettes; the other four genomes cannot be used to test it.
- Held-out homology is not informative outside the two curated Agaricomycete species: with the own-species curated precursors removed, *S. commune* drops from 11 to 0 homology clusters in the array and *C. cinerea* from 6 to 1; *Grifola* and *Russula* have 0. Precursor homology across orders is too diverged (E <= 1 tblastn of 31 curated sequences, mostly from other classes) to count precursors; the usable precursor signal is the CAAX ORF, which is circular for labelling.
- Array size and span do not distinguish: in the full set, 700 arrays have >= 3 loci (in 595 genomes), 311 have >= 4, so 3-5-locus arrays are about as common as genomes. 3 of the 6 B arrays have size <= 3; one is a singleton. Where the B array sits in the distribution: *S. commune* is in the top 25 arrays by span and top 13 by CAAX ORFs; *Heterobasidion* is at the 757th largest span and has 1 CAAX ORF (its record has 1 strict-CAAX precursor out of 3 by curation: the curated precursors are mostly relaxed CAAX or CAAX-like motifs not caught by the strict rule).
- Pipeline calls: the pipeline-called span for *Trametes* is 62 kb (1,556,459-1,618,028), covering the curated B region plus three other paralogous loci; the other calls roughly match the curated region (*Heterobasidion* exactly, *Grifola* 44 against 35 kb). So existing calls already merge a B locus with adjacent paralogs in at least 1 of 5.

### 4. Cassette over the 1,287 genomes (`cassette_by_order.tsv`, `summary_output.txt`, `extra_output.txt`)
Receptor loci (12,518): cassette A 1,101 (8.8%), B 626, C 199. Chance by contig-uniform placement: 22 expected for A (2%), in every order below 10 expected. Genomes with at least one cassette-A locus: 613 (47.6%); B 370; C 134.
| Order (genomes) | Loci | Cassette A (chance) | Cassette B | Cassette C | Genomes with A | PR-called loci with A / uncalled with A |
|---|---|---|---|---|---|---|
| Agaricales (676) | 7,032 | 708 (7.8) | 463 | 167 | 391 | 486 / 222 |
| Polyporales (175) | 1,473 | 230 (9.0) | 117 | 28 | 126 | 220 / 10 |
| Boletales (182) | 1,450 | 80 (2.6) | 22 | 2 | 45 | 54 / 26 |
| Russulales (76) | 993 | 41 (1.2) | 15 | 2 | 24 | 26 / 15 |
| Auriculariales (14) | 85 | 18 (0.2) | 0 | 0 | 11 | 18 / 0 |
| Cantharellales (77) | 674 | 2 (0.6) | 1 | 0 | 2 | 2 / 0 |
| Hymenochaetales (47) | 425 | 2 (0.2) | 1 | 0 | 2 | 2 / 0 |
| 12 other orders (40) | 386 | 20 (0.3) | 7 | 0 | 12 | 20 / 0 |
- The cassette is far above chance (1,101 against 22), and it is concentrated: 5 orders carry 98% of cassette-A loci; Cantharellales (99 calls) and Hymenochaetales (13 calls) have almost none. The CAAX-based cassette therefore cannot confirm the receptor calls in those two orders (73 of the 99 Cantharellales calls were already unsupported by array or homology). Either those orders lack CAAX-type precursors or the receptors are not B-locus receptors; the data do not decide it.
- Cassette is a locus-level unit, not an array unit: cassette-A loci per genome: 0 in 674 genomes, 1 array in 525, 2 in 79, 3 or more in 9. In genomes that have one, 86% (525 of 613) have exactly one cassette array, as expected for a single B locus per genome with multiple cassettes inside it; arrays with >= 2 cassette-A loci: 315 (of 1,864 arrays of size >= 2). This is consistent with the B locus being a cassette cluster but is not a test of it.
- By array size, share of arrays with any cassette class (A, B or C): singletons 2.2% (164 of 7,572), size 2 15.6%, size 3 49.1%, size 4 65.0%, size >= 5 64.1%. The 164 singleton arrays with a cassette are 25.6% STE20 co-located against 4.1% for singletons without.
- Semi-independent check (STE20 top hit within 250 kb, arrays in genomes with an STE20 hit; the flank was chosen from non-Agaricomycete MAT linkage, so this is supporting, not a label): size >= 2 arrays 48.7% with a cassette (n = 573) against 19.6% without (n = 1,283); size >= 3 with a cassette 56.1%. Among size >= 2 flagged arrays, 49.6% with a cassette (n = 548) against 34.9% flagged without (n = 456) and 11.6% for unflagged (n = 852). So the cassette separates within the CAAX flag, but the effect is smaller than flagged against unflagged.
- Uncalled cassettes: 273 loci (200 arrays, 189 genomes) have a cassette A and are not in a PR call; the contig-uniform null expects 12.6 loci among all uncalled loci, so the uncalled cassettes are essentially all above chance (96%). 138 arrays have >= 2 loci; 60 (57 genomes) are class C (each CAAX ORF with homology); 141 arrays are class C or size >= 2 and A. 125 of the 189 genomes have no PR call anywhere (class C: 52 of 57; C or size >= 2: 102 of 137). Agaricales carry 169 of the 205 arrays with any class of cassette and no call.
- Not shown, and not claimable: that a cassette is a B locus rather than any CAAX-bearing STE3 cluster, or that its CAAX ORFs are pheromones. It does support that the 273 uncalled loci are not random CAAX coincidences.

## What cannot be concluded
- Whether any array property predicts the B locus: labelled n = 6 genomes (2 independent). Four of six labels were chosen from CAAX positions.
- Array span or size as a B-locus boundary: arrays include paralogs and span 2-5 times the true locus; tandem loci are merged by the inventory step, so loci undercount receptors.
- Precursor counts per cassette: the strict-CAAX motif does not match all curated precursors (curated *Heterobasidion*: 1 strict of 3; *Trametes*: 7 of 9), and cross-order homology is weak (held-out *S. commune* 0 of 11).
- Flank genes: none established in Agaricomycete B loci. The tested flank (STE20) is not B-specific here.
- Sublocus labels (Balpha/Bbeta) from array data: the *S. commune* array shows two cassettes, but a single species cannot define a rule.
- Anything for Cantharellales, Hymenochaetales and the 12 minor orders (cassettes ~0; 14, 29 and 31 species).
- Outgroup receptor loci were not scanned.

## Options (estimates on the 1,287-genome scanned set; no code change here)
| # | Option | Purpose | Expected change | Risk | Validation needed |
|---|---|---|---|---|---|
| 1 | Report cassette fields on each PR call and array: receptor loci count (note: lower bound), CAAX ORFs within 5 kb, homology clusters, `cassette` A/B/C | Show which calls have a precursor cluster behind them; no label change | 0 calls changed; 828 of 2,311 called loci (36%) carry cassette A, 1,483 do not (Agaricales 486 of 1,237; Cantharellales 2 of 119) | Low. Needs the 5-kb window fixed, and a note that CAAX-based cassettes are circular for admission | Unit tests on a few records; check *S. commune*/*C. cinerea* against the curated cassettes |
| 2 | Cassette records for B-locus calls: a record per array carrying cassette positions (and sublocus count as the number of cassette loci) | One B-locus record per array with its cassette structure instead of a flat span; makes the *S. commune* two-cassette pattern visible without claiming Balpha/Bbeta | 1,039 calls to about 1,028 array records (per previous note); about 315 arrays with >= 2 cassette loci get a multi-cassette flag | Low-medium: sublocus naming must stay with the curated Balpha/Bbeta families; "two cassettes" is not "two subloci" | Hand check of the 315 arrays; compare against Balpha/Bbeta calls in Agaricales |
| 3 | B-locus-array confirmation tier (`array_support: cassette`): cassette A or C in the array of a call, shown beside `unverified`, never lifting it | Separate calls with a precursor cluster from those with receptor + one CAAX ORF | At locus level 828 of 2,311 called loci have cassette A (calls were not tabulated separately from loci); the called loci without one include nearly all in Cantharellales and Hymenochaetales | Medium: circular if used to confirm (CAAX admitted them); gives no information in 2 orders | Independent labels (>= 10 genera with mapped B receptors) before it is allowed to change a label |
| 4 | Admit uncalled cassette arrays as `array_candidate`, not as PR calls | Recover the receptors that the call rules withhold | 200 arrays in 189 genomes, 125 genomes with no PR call at all; restricted to class C or size >= 2: 141 arrays in 137 genomes (102 first records); chance about 5% (12.6 expected of 273 uncalled cassette loci, under a contig-uniform null that ignores gene density) | Medium: Agaricales 82% (169 of 205); same CAAX dependence as the call it supplements; no independent label | Label a sample from independent evidence (mapped receptors in *Pleurotus*, *Flammulina*, *Lentinula*); hold out by order; STE20 co-location by tier (already 49%) |
| 5 | Use cassette presence to pick the B receptor copy or array in a genome with several arrays | Choose among paralog arrays | A unique cassette array in 525 of 613 genomes; 88 genomes have two or more | High: n = 2 independent labels; the B array also holds paralogs | Not supportable yet; needs labels from >= 10 genera |
| 6 | Make sublocus labels (Balpha/Bbeta-like) from two cassettes in one array | Report B-locus architecture | 79+9 genomes have two or more cassette arrays; 315 arrays have >= 2 cassette loci | High: one curated species; Casselton 1998's subloci are a *S. commune* property | Do not do |

## Recommendation
1. Report, do not call. Options 1 and 2 (cassette fields, one record per array with cassette positions) cost no calls, make the B-locus structure visible, and use only the existing scan outputs. Keep the cassette separate from `confidence` and from `unverified`.
2. If option 4 is wanted, restrict it to class C or size >= 2 arrays in the five orders where it exceeds chance (141 arrays, about 4-5% chance expected), outside the PR call count, and only after labelling a few mapped-B genomes independently. It does not add independent evidence: it is the CAAX signal re-read at array scale.
3. Do not use array span, array size, flank genes or held-out precursor homology to define or confirm a B locus; do not use cassette to change any label before there are labelled B loci that were not chosen from CAAX positions.
4. Required before anything stronger: labelled B loci from at least 10 genera that are not selected by CAAX (mapped receptors; for example *Pleurotus*, *Flammulina*, *Lentinula*), a receptor-locus model that does not merge tandem receptors (loci undercount receptors, 9 loci for 16 curated), and a Cantharellales/Hymenochaetales look at their precursor motifs (cassette near 0 although 158 loci are called).

## Reproduction
`python b_curated.py; python b_cluster.py <extracted scan_out dir>; python b_summarise.py; python b_extra.py` in `results/2026-10-06_b_locus_clustering/` (pandas, numpy, pyyaml, scipy). The held-out step: `sbatch run_heldout.sh` on HPCC from a copy of `scan_genome.py`, `pheromones_curated.faa` and `ste3_all.faa`; job 29530344.
