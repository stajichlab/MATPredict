# C. albicans MTL: read depth vs the assembly call

2026-09-26. Data: `depth2.tsv` (map2.slurm), and `depth.tsv` (map.slurm) for two
isolates. Script: `analyze.py`. Table: `zygosity.tsv`. Assembly calls:
`results/2026-09-26_serinales_all_882aa01/runs/*/detection_report.yaml`.

## Question

Do C. albicans assemblies that show one MTL idiomorph hide an a/alpha
heterozygote?

## What was mapped

- References: SC5314 MTLa (AF167162.1) and MTLalpha (AF167163.1) locus
  sequences, and a 20 kb chr1 control (chr1:1,000,001-1,020,000).
- Reads: ENA R1 FASTQ, one WGS run per isolate.
  - Round 2 (14 isolates) sampled about 3M reads at random from the whole file.
  - Round 1 took the first 2M reads. In 4 of its 6 isolates the control was
    49-121x and the MTL was 0x, because the file was sorted by position. Those
    4 are discarded. The 2 round-1 isolates with a normal control (18x) are
    used.
- Mapping: minimap2 `-ax sr`, MAPQ >= 20, primary alignments only.
  Depth = mean `samtools depth -a` over each region.

## Scoring

- The idiomorph ratio is the mean depth of the idiomorph blocks divided by the
  control depth. For a, the blocks are MTLa2-MTLa1 and PAPa-OBPa-PIKa. For
  alpha, they are MTLalpha2, MTLalpha1 and OBPalpha-PIKalpha-PAPalpha.
- An idiomorph is present at ratio >= 0.20 and absent at ratio < 0.05.
- Why these thresholds: in a diploid, a heterozygote should give about 0.5 per
  idiomorph, and a homozygote about 1.0 and 0. The observed values fall in
  two clear groups: 0.43-0.59 for each idiomorph, or 0.90-1.02 against 0.00.
  No isolate falls between 0.05 and 0.20. Any cut-off in 0.05-0.40 gives the
  same result.

## Cross-mapping check

- a against alpha references (minimap2 asm20): the only shared sequence is
  about 1.5 kb at each outer end, at 99% identity. The idiomorph-specific
  cores do not align. Reads in the shared ends get MAPQ 0 and are dropped.
  258 bp of the MTLa2-MTLa1 block lies in the right-hand shared end. That
  makes the a ratio lower, by up to about 15%.
- Measured: the 5 a/a isolates give an alpha ratio of 0.00 (the highest block
  is 0.23x against a 16x control). Cross-mapping is too low to affect any call.

## Result (n = 16)

| assembly call | reads a/a | reads a/alpha |
|---|---:|---:|
| A (one idiomorph) | 5 | 1 |
| alpha (one idiomorph) | 0 | 7 |
| A+alpha | 0 | 3 |

- Of the 11 read-based heterozygotes, 8 have an assembly that shows only one
  idiomorph: 7 alpha-only and 1 A-only.
- No assembly showed both idiomorphs when the reads showed one.
- **Every one of the 7 alpha-only assemblies is an a/alpha heterozygote by
  reads.**
- GCA_000773825.1: the assembly has two A loci on two contigs (AJJA01000022.1
  and AJJA01000023.1). By reads it is a/alpha.

## What can be concluded

- The assemblies in this set collapse the heterozygous MTL. An assembly that
  shows one idiomorph is not evidence of a homozygote.
- The collapse is not random with respect to idiomorph. Of 8 collapsed
  assemblies, 7 kept alpha and 1 kept a.
- MATPredict found the idiomorph that is in the assembly. It did not miss a
  locus that the assembly contains. Its detection is not at fault. The
  assembly is.

## Limits

- n = 16. The isolates were chosen by hand to cover the three assembly call
  classes and to have WGS reads. They are not a random sample of the 134
  C. albicans genomes. So the collapse rate across all 134 (8/11 here) must not
  be used as an estimate.
- Six of the seven alpha-only heterozygotes are one submission series
  (GCA_000447455.1 to GCA_000447555.1, SRR845xxx). They were probably
  assembled the same way. The idiomorph bias may be a property of that
  pipeline.
- One read run per isolate, and 13-18x control depth from 3M sampled reads.
- One 20 kb control region, on chr1. Aneuploidy of chr5 (where MTL lies) is
  common in C. albicans, and it would shift the ratios. The two clean groups
  (about 0.5/0.5 and about 1/0) suggest it did not affect these 16.
- Three isolates in `isolates.json` had no read run and were not tested.

## Curator decision needed

C. albicans genotypes that come from assemblies are not reliable for zygosity.
Two options: (a) report a single-idiomorph call in this species as "one
idiomorph found in the assembly, zygosity unknown"; (b) add a read-based
zygosity check wherever reads exist. With either option, the 134-genome a+alpha
share in the Serinales scan should not be read as a population frequency.

## Curator ruling (2026-09-26)

Report C. albicans single-idiomorph assembly calls as zygosity unknown. A
better treatment (e.g. read-based zygosity) will be needed later; not scoped.
