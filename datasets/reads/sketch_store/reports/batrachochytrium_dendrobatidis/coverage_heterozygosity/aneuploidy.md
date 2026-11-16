# Bd chromosome depth changes and allele balance (2026-11-14)

Status: exploratory; nothing is curated. Source: `bd_depth_by_seq_JEL423.tsv` (69 isolates against the 20-sequence JEL423 assembly, job set of
2026-11-13), scripts `scripts/bd_aneuploidy.py` and `scripts/bd_aneuploidy_af.py` (tests in `tests/test_bd_coverage.py`). Tables in `tables/`:
`aneuploidy_flags.tsv`, `aneuploidy_af.tsv`, `genomewide_allele_balance.tsv`, `aneuploidy_by_sequence.tsv`.

## Method
Sequences of 500 kb or more (16 chromosomes). The depth ratio of a sequence (median depth over the genome median) is divided by the median ratio of the
same sequence over the isolate's lineage class (CLFT-type 15 isolates, JEL-type 52; the two HET-rich isolates are too few to form a class). This removes
lineage-specific offsets. An isolate-sequence is flagged when the normalised ratio is outside 1 +/- 0.25 and the isolate has a median depth of 20x or more.
Then, for each flag, the heterozygous SNPs (minor allele fraction 0.15 to 0.5) of that sequence are compared with those of the isolate's other chromosomes:
the share near 1/3 (0.27 to 0.40) and near 1/2 (0.43 to 0.50). An extra copy of a chromosome in a diploid should move sites from 1/2 to 1/3.

## Result
- Population-wide offsets (not isolate-specific): median ratio 1.28 on CP161924.1 and 1.23 on CP161925.1 across all isolates, 0.77 on CP161938.1. They are properties of the reference or of the species relative to JEL423, not aneuploidy. The small sequences CP161939 to CP161942 (57 to 141 kb) have ratios of 0, 0, 3.9 and 0.73: repeats or organelles, excluded from the analysis.
- 48 isolate-sequence flags in 28 of 69 isolates: 35 gains (norm > 1) and 13 losses; 10 gains at 1.4 or more, 2 losses at about 0.52.
- Gains (24 with norm 1.3 or more): 15 show the allele-balance shift (1/3 share up and 1/2 share down by more than 0.1 against the isolate's other chromosomes); 9 do not. The clearest cases are in high-depth isolates (50x or more): L2203 CP161926.1 (ratio 1.46; 66% of het sites near 1/3 against 22.5% on other chromosomes, 7% near 1/2 against 50%), KRBOOR_331 CP161926.1 (1.45; 76% against 22%) and RTP5 CP161931.1 (1.38; 68% against 24%). That is the signature of one extra copy of one chromosome in a diploid isolate. The same isolates have 50 to 57% of het sites near 1/2 on their other chromosomes, which is the diploid pattern.
- CLFT065 CP161936.1 (ratio 1.59, 55x) shows the opposite shift (53% of het sites near 1/2 against 30% on the rest), so it is not the trisomy pattern; I cannot say what it is (a duplication, or a different copy state).
- Losses (7 with norm 0.7 or less): in 4 the het SNP density on the flagged chromosome is below 0.75 of the rest, as expected when one copy is lost (CLFT078 and CLFT166 CP161938.1 at ratio about 0.52: 831 and 698 het per Mb against 1,866 and 1,078).
- CP161926.1 is gained in 6 isolates from different lineage classes and countries (Montserrat L2203, South Korea KRBOOR_331, Sweden, and three Brazilian isolates); CP161936.1, CP161937.1 and CP161931.1 are gained in 4 isolates each. I did not test whether these gains are shared by descent.
- Isolates with several flagged chromosomes: CLFT065 (5), CLFT080, CLFT156, CLFT165, CLFT167 (3 each) and the "Generic sample" from Switzerland (4). I did not examine them further.

## What this does not show
- A baseline ploidy for each isolate. The genome-wide allele balance puts 21 isolates in a "third-dominated" group, 15 "half-dominated" and 33 "mixed", but this label is confounded: the third-dominated isolates have a median depth of 21x against 105x for half-dominated ones. At low depth the allele fractions are noisy and reference divergence skews them, so I do not read the label as triploidy. Only the within-isolate comparison (flagged chromosome against the same isolate's other chromosomes) is controlled for depth.
- A whole-genome ploidy change is invisible in a ratio to the genome median.
- The 20x cutoff leaves noise at the lowest depths; 10 isolates are below 20x and were not flagged.
- The two HET-rich isolates are covered in the next section; CLFT039 has a median depth of 16.2x, below the 20x cutoff, so it cannot be flagged.

## The two HET-rich isolates (CLFT160 = SRR10698391, CLFT039 = SRR10698425): hybrid-like or aneuploid?
`scripts/bd_aneuploidy.py --fallback-class JEL-type` compares them with the JEL-type median (two isolates cannot form a class; the JEL-type baseline is a rough
reference). `scripts/bd_hetrich_check.py` gives per-sequence depth, allele balance (JEL423 reference) and parent-of-origin window classes (CLFT044 reference).
Tables: `aneuploidy_flags_with_hetrich.tsv`, `hetrich_by_sequence.tsv`.
- Depth: CLFT160 has two gain flags (CP161931.1 at 1.32 and CP161936.1 at 1.38 against the JEL-type median); CLFT039 (16.2x) is below the flag cutoff, and its ratios lie between 0.77 and 1.24. Neither isolate has a loss flag.
- Parent-of-origin windows on the 14 scaffolds of 500 kb or more: heterozygous windows (JEL423 and CLFT044 alleles both present) are the majority on almost all of them: median share 0.70 in CLFT160 and 0.72 in CLFT039. In CLFT039 the share is 0.47 or more on every scaffold; in CLFT160 it is 0.40 or more on 12 of 14. The two exceptions are JAQIHF010000010.1 (share 0.105; 32% CLFT windows) and JAQIHF010000014.1 (share 0.000; 50% CLFT windows); their depth is normal (0.99, 0.93), so they are not copy-number changes. Windows homozygous for the JEL423 allele are almost absent (median share 0.01 and 0.05); CLFT-class windows are 29% and 20%.
- So the heterozygosity of both isolates is genome-wide and is not tied to a gained chromosome: CLFT160's two gains are on scaffolds with ordinary HET shares. The pattern fits an isolate that carries both lineages in one genome with part of the genome homozygous for the CLFT044 allele (possible loss of heterozygosity on two scaffolds of CLFT160), and does not fit aneuploidy as the cause of the heterozygosity.
- Limits: the allele-balance columns of the same table (median 0.37 and 0.33 of het sites near 1/2) are at depths of 27.5x and 16.2x and are too noisy to read; the "hybrid" reading rests on the window classes, which depend on the 0.5 recovery cutoff set from two control isolates; a third lineage that shares alleles with both references would give the same pattern; the two isolates are the only ones of this type, so no comparison group exists.
