# Uncalled R. arrhizus genomes and the split-locus rule
Status: decided

## Results
- Diagnosis (`results/2026-09-27_rarrhizus_uncalled/NOTE.md`): 12 uncalled
  GL-series genomes. Correction: 11 of the 12 carry sexP at 98.4% alone near a
  contig end; GL5 GCA_011764265.1 has no sexP hit (the note's "all 12" claim
  was wrong — see `results/2026-09-27_rarrhizus_uncalled/best_hits.tsv`).
- Split-locus rule (4174440): one modelled core >=95% within 500 bp of a contig
  end plus >=2 flanks on other contigs -> partial_locus, low, flagged split.
  Mucoromycota 235->247 genomes; Cryptococcus, Serinales sample and
  Dothideomycetes unchanged; Zygo 23/23 (`results/2026-09-27_split_locus/NOTE.md`).

## Limits
Tuned on one Rhizopus sequencing batch.
