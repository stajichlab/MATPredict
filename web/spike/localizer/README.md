# Localizer evaluation (spike)

Compares candidate replacements for the `tblastn` localization step against native BLAST+ 2.17.0
`tblastn` on the same inputs: the 468 curated reference proteins against `testset/MAT_loci3.fasta`
(366 kb, 189 contigs). Only strong `tblastn` HSPs (bit score >= 33, the pipeline's supported
threshold) are scored.

## Files

| Path | Purpose |
|---|---|
| `bin/tblastn`, `bin/makeblastdb` | Stand-ins that run MMseqs2 and emit the pipeline's 11-column outfmt. Use with `PATH=<localizer>/bin:$PATH`. |
| `compare.py` | miniprot vs tblastn (recall, bit/AS fit, identity and coverage). Expects `tblastn.tsv` and `miniprot.paf` in the working directory. |
| `evaluate.py` | Scores a candidate TSV (`mmseqs` or `blast` format) against `tblastn.tsv`. |

## Results

| Localizer | HSP recall (bits >= 33) | Query recall (209 strong queries) |
|---|---|---|
| miniprot (default) | 3.2% | 40 (19%) |
| MMseqs2 default sensitivity | 76.0% | 187 (89.5%) |
| MMseqs2 `-s 7.5 -e 10` | 86.8% | 208 (99.5%) |
| MMseqs2 `--exhaustive-search` | 75.3% | 202 (96.7%) |

MMseqs2 bit scores run about 5% below BLAST+ on strong hits (median ratio 0.95). It misses 57%
of 20-50-bit HSPs. Exhaustive mode did not improve on the default.

## Pipeline comparison (MAT_loci3, native miniprot and exonerate, native reference)

With MMseqs2 as the localizer, `detect` returns the same single locus (`OQ436909.1`, Ascomycota:MAT,
MAT1-2, medium). Withheld candidates differ (153 vs 262), and some not-detected reasons differ.
Full comparison in the chat record; the reports were written under `~/wasm-build`, not the repo.
