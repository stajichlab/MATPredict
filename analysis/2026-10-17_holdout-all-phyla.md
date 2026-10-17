# Held-out test across the Ascomycota subphyla, Basidiomycota and Mucoromycota (2026-10-17)

Code: `scripts/run_asco_holdout.py` and `scripts/build_asco_holdout_manifest.py` (main at 2026-10-17, with the
profile builder inlined). Manifest: `testset/holdout/holdout_manifest.tsv` (166 curated records, 88 with a genome
assembly). Each task withholds one record's taxon at the genus, family or order radius
(`detect.holdout.records_to_withhold`), runs `detect --taxid` on that record's own genome with those records
excluded, and scores against the record's own mating type. 210 tasks ran; combinations that leave no reference of the
record's family (60) were not run. Rows and summary: `results/2026-10-17_holdout_all_phyla/`.

Descriptive: the sets are small and no significance is claimed. A record's own mating type is the truth.

## Results (pipeline call; `correct` = the record's idiomorph, either one for a homothallic record)

| family | subphylum | genus | family radius | order |
|---|---|---|---|---|
| Ascomycota `MAT` (6 heterothallic, 3 homothallic) | Pezizomycotina | 9/9 | 9/9 | 9/9 |
| Ascomycota `MTL` | Saccharomycotina | 6/6 | 5/6 | 1/6 |
| Ascomycota `MATsc` | Saccharomycotina | 3/5 | no reference | no reference |
| Ascomycota `PM` | Taphrinomycotina | 2/3 | 2/3 | 2/3 |
| Basidiomycota `redPR` | | 5/5 | no reference | no reference |
| Basidiomycota `bLocus` (pattern family; locus found) | | 22/22 | 22/22 | 21/22 |
| Basidiomycota `aLocus` (locus found) | | 8/8 | 0/8 | 0/8 |
| Basidiomycota `rustHD` / `HD` / `PR` (locus found) | | 4/4, 1/1, 2/4 | 4/4, 1/1, 2/4 | no reference, 0/1, 2/4 |
| Mucoromycota `MAT` | | 9/9 | 9/9 | 4/9 |

- **Mucoromycota independent set, Zygo 23** (current code, never in training): locus on the truth scaffold 23/23 and
  idiomorph correct 23/23 on scaffolds and on contigs. The Mucoromycota held-out rows above are optimistic: the
  idiomorph classifier is prebuilt from every curated record (`classifier_leak` in the rows).
- **Alpha-box and HMG-box profiles, built without the withheld records:** on the 6 heterothallic Pezizomycotina
  records the profile call is correct at every radius (18/18). The three homothallic records are consistent with
  either idiomorph at genus radius and come back undetermined at family and order radius for the two Xanthoria
  genomes (their near references are withheld).
- **Order radius is where it fails, as expected.** With the whole order withheld, Mucorales loci (6 records) are not
  found with the 1 to 2 references of another order that remain (Umbelopsidales, found when Mucorales references
  stay); Serinales `MTL` loci are not found with one reference left (5 of 6). Pezizomycotina `MAT` has references
  across many orders and is 9/9 even then.
- **Misses at every radius:** *Taphrina deformans* `PM` (2 references left, both Pneumocystis). Basidiomycota `aLocus`
  is found at genus radius but never when the taxonomic family is withheld (0/8).
- **`MATsc` wrong calls (2 of 5, genus radius):** the *S. cerevisiae* HMRa and *K. lactis* MATa records. In both genomes
  the active locus is MATalpha and that is what is called; the record is the silent or alternate cassette. The runner
  scores by the record's own idiomorph, so this reads as wrong; I have not checked whether it is a detector fault.

## Performance
210 tasks (detect on one genome plus building and scoring the two profiles where applicable): wall time median
71 s, 90th percentile 239 s, maximum 504 s; peak memory 0.29 GB.

## Caveats
- A first summary counted references across phyla (the locus name `MAT` is shared); the manifest and `summarize
  --manifest` now compare phylum and family. The tables above use the corrected counts.
- Pattern-vocabulary families (`bLocus`, `aLocus`, `HD`) are scored only on whether a locus is found; allele calls are
  not tested here.
- Pezizomycotina heterothallic held-out genomes number six, so the 18/18 profile and 27/27 pipeline figures carry
  little statistical weight. Independent truth (for example BioSample mating-type attributes) would enlarge the set.
