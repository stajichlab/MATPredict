# 011. Candida albicans assemblies collapse heterozygous MTL

- **Category:** assembly-or-annotation artefact
- **Status:** verified (read depth)
- **Lineage:** Ascomycota, Serinales

## Summary
Read depth shows that 8 of 11 a/alpha heterozygous isolates have an assembly
with only one idiomorph, mostly alpha. So a single-idiomorph C. albicans
assembly is not evidence of a homozygote.

## Evidence
- `results/2026-09-26_calbicans_mtl_reads/NOTE.md`, `zygosity.tsv` (committed
  on polish-scope-cuts `cd08a0e`).
- References: SC5314 MTLa AF167162.1, MTLalpha AF167163.1; chr1 control
  chr1:1,000,001-1,020,000. Idiomorph/control depth ratios: 0.43-0.59 each
  (heterozygote) or 0.90-1.02 vs 0.00 (homozygote); none between 0.05 and 0.20.
- Assembly call vs reads (n=16): A-only 5 a/a + 1 a/alpha; alpha-only 0 + 7;
  A+alpha 0 + 3. Of 8 collapsed assemblies, 7 kept alpha.
- Six of the seven alpha-only heterozygotes are one submission series
  (GCA_000447455.1-GCA_000447555.1). Annotation report B4.1, B4.2.

## Method that found it
minimap2 `-ax sr`, MAPQ >= 20, ~3M sampled reads per isolate (`map2.slurm`,
`analyze.py`).

## Verification done / still open
Done: read depth. Open: a random sample of the 134 genomes; the curator ruled
single-idiomorph C. albicans calls report zygosity unknown (commit `1ac21d9`).

## Limits
n=16, hand-chosen; one control region; chr5 aneuploidy could shift ratios.

## Related
Entry 007.
