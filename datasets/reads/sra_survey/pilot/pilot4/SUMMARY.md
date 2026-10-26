# Summary

Runs listed 80; done 79; failed 1; not run 0.

| group | runs | done | MAT1-1 | MAT1-2 | both | none | low_depth | with truth | typed | agree | disagree | untyped | GB | MB/s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| aspergillus | 20 | 20 | 9 | 11 | 0 | 0 | 0 | 16 | 16 | 16 | 0 | 0 | 37.6 | 14.7 |
| botrytis | 20 | 20 | 9 | 9 | 2 | 0 | 0 | 11 | 11 | 11 | 0 | 0 | 35.6 | 13.7 |
| pyricularia | 20 | 19 | 7 | 10 | 0 | 0 | 2 | 16 | 14 | 14 | 0 | 2 | 30.0 | 15.2 |
| verticillium | 20 | 20 | 1 | 17 | 0 | 2 | 0 | 10 | 10 | 10 | 0 | 0 | 72.0 | 14.8 |

Truth = the `detect` idiomorph of a BFD assembly with the same BioSample (or a mating type in the strain name). `typed` = a one-idiomorph call; `untyped` = none or low_depth on a run that has a truth.
Discordant runs: `discordant.tsv`; runs without a result: `failures.tsv`.
