# BioSample truth set under withheld-record held-out runs (2026-10-17)

Code: `scripts/build_truth_holdout_manifest.py`, `scripts/run_asco_holdout.py` (runner patched to take a row's drop list),
manifest `testset/holdout/biosample_truth_manifest.tsv`. Rows and summary: `results/2026-10-17_biosample_truth_holdout/`.
Truth: the 203 Pezizomycotina genomes whose BioSample states MAT1-1 or MAT1-2 (analysis/2026-10-17_biosample-truth-set.md);
none is a curated record. Each genome is run at genus radius (curated records of its genus withheld; 118 of the 203 have no
curated record in their genus, so for them this is no withholding) and order radius (curated records of its order withheld).
`detect --taxid` on the genome, with the withheld records excluded; the alpha-box and HMG-box profiles are rebuilt without them.
The stated mating type is the submitter's, so a disagreement is read, not assumed to be an error. Descriptive; no significance claimed.

## Result (pipeline call against the stated mating type)

| radius | correct | wrong | undetermined | no locus |
|---|---|---|---|---|
| genus | 178 | 21 | 0 | 4 |
| order | 169 | 4 | 23 | 7 |

- **Genus radius, wrong 21:** 17 *Botrytis* (the known reference gap: no *B. cinerea* MAT1-2 record; the truncated MAT1-1-1
  remnant wins the bitscore vote), 2 *Penicillium roqueforti*, 1 *Podospora tetraspora*, 1 Magnaporthales sp. No locus: 2
  *Cryphonectria parasitica*, 1 *P. roqueforti*, 1 *Morchella conica*.
- **Order radius:** the *Botrytis* calls become correct or uncalled when Helotiales references are withheld (the sas56 MAT1-1
  record that caused the wrong vote is gone; 1 *Botrytis* becomes no locus); the 23 undetermined and 4 of the 7 no-locus are
  all *Cryphonectria* (Diaporthales has one curated record, so withholding the order leaves only distant references).
- **Profiles** (alpha box and HMG box, rebuilt without the withheld records), 199 genera-radius loci scored: 194 correct, 4 wrong,
  1 undetermined; order radius: 192 correct, 3 wrong, 1 undetermined. The profile call is correct on all 17 *Botrytis* genomes,
  where the pipeline's vote is wrong. The profile errors are *Penicillium* (2), Magnaporthales (1) and a *Juglanconis*
  (genus radius only), and they overlap the pipeline's own errors there.

## Reading
The pipeline agrees with the submitter's statement on 88% of genomes at genus radius and 83% at order radius, and 96% (178 of
186) without *Botrytis*. The profile agrees on 97% (194 of 199 scored). The largest error group is a database gap that the
profiles already resolve. Order-radius failures come from clades with a single curated record (Diaporthales).

## Caveats
- Truth is submitter-stated, not verified; the *P. roqueforti*, *Podospora* and Magnaporthales disagreements are not examined.
- The genus radius is no withholding for 118 genomes, so genus and order differ little there.
- Scored profile loci exclude genomes with no locus; profile training sets shrink at order radius.
