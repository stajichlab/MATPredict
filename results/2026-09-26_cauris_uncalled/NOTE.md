# Why 12 C. auris genomes stay uncalled (Serinales scan, code 882aa01)

2026-09-26. Scripts and raw output are in this directory.

## Count

- Scan `2026-09-26_serinales_all_882aa01`: 952 C. auris genomes, 12 uncalled.
  (The earlier `8dec9e1` scan, before the flank roster, had 271 uncalled.)
- The handoff said "six". The real number on 882aa01 is 12.

## Method

- Query: the 13 proteins of the three C. auris records at 882aa01
  (B11205 MTLa, B11243 MTLa, B11221 MTLalpha), including PAP1/OBP1/PIK1.
- tblastn with genetic code 12 (`run_direct.sh`), best hit per protein.
- Controls: B11205, B11220, B11221, B11243, B8441.
- Locus position: blastn of a 55 kb window around the MTL of B11220 (alpha)
  and B8441 (a), then N-run scan of the matching window (`region.py`).
- Assembly method from NCBI Datasets for all 952 genomes (`asm_methods.tsv`).

## Controls behave

Every control hits its own idiomorph at 97-100% over full length. In C. auris
the flanks PAP1/OBP1/PIK1 are idiomorph-specific alleles: the other idiomorph's
allele hits at only 51-68%.

## The 12 uncalled genomes: the idiomorph is not in the assembly

In all 12, no MTL gene and no idiomorph-specific flank allele is present. The
only hits are OBP1/PIK1 paralogs at 28-41%. Detection withholds these, which is
correct. **No detection defect.**

| Genome(s) | Strain | Assembly method | What is at the MTL position | Cause |
|---|---|---|---|---|
| GCA_041026{465,485,505,525,545,565,585,605} (8) | SCO_240..279, Portugal, PRJNA1140147 | CLC 12.0.3, "Chromosome" level; chromosome lengths match B11220 (clade II, alpha) | 8.1-8.2 kb of N runs at 2,562,2xx-2,570,5xx: the exact span of B11220's MTLalpha idiomorph | Reference-guided consensus. Reads did not fill the reference idiomorph |
| GCA_027563775.1, GCA_027563725.1 | 15, 204, England, PRJNA865124 | BWA-MEM 0.7.17 (reference mapping) | 8.8-8.9 kb of N runs at contig 5, 845.9-854.9 kb: the span of the B8441 MTLa idiomorph | Same: reference consensus |
| GCA_902703435.1 | CA1 (CA_KFCM1), Saudi Arabia | SPAdes, 250 contigs | Contig break at the idiomorph: flank on the end of contig 106 and the start of contig 10. The raw NCBI assembly (all 250 contigs) also lacks it | Idiomorph not assembled |
| GCA_046252445.1 | 27-CA, Pakistan, PRJNA1164417 | SPAdes, 12 kb total | Not a genome: 3 contigs of repeat sequence (multi-copy hits across B8441 chromosomes) | Not a genome assembly |

## Inferences, with limits

- In a reference consensus, an N gap exactly over the reference idiomorph
  most likely means the isolate carries the OTHER idiomorph, because its
  reads cannot map across non-homologous idiomorph sequence. So the SCO
  isolates are probably MTLa and the two England isolates probably
  MTLalpha. Low coverage could produce the same gap. Only the reads can
  confirm it. Not tested.
- CA1 is probably MTLa: the only idiomorph-specific sequence left is a
  131-aa PAP1 fragment at the start of contig 10. It is 100% identical to
  the a allele (B11205) and 52% to the alpha allele. One fragment is thin
  evidence.

## Wider risk: reference-guided assemblies

Of 952 C. auris assemblies, 8 are BWA-MEM (PRJNA865124), 22 are CLC (two
projects), and 1 is samtools. When the isolate carries the same idiomorph
as the reference, a reference consensus reports it correctly. But the call
then depends on how the pipeline treats uncovered bases. These pipelines
wrote N, which is safe. A pipeline that fills uncovered bases with
reference sequence would produce a false call of the reference idiomorph.
Not checked beyond these projects. PRJNA865124: 6 called A, 2 uncalled.
PRJNA1140147: 8 uncalled. PRJNA640677 (CLC 21.03): 14 called A.

## Recommended actions

1. Add GCA_046252445.1_ASM4625244v1 to the BFD suppress list: it is 12 kb of
   repeats, not a genome.
2. Do not count the 10 reference-consensus genomes as detection misses.
   Report them as "idiomorph region is an assembly gap". A cheap
   detector-side option: when a flank-anchored locus position holds a long N
   run, report `assembly_gap_at_locus` instead of plain not-detected. Not
   implemented; the curator decides.
3. Optional: resolve the 10 from SRA reads by mapping them to both idiomorphs.
   The C. albicans read work (results/2026-09-26_calbicans_mtl_reads/) shows
   the pattern.
4. Consider flagging reference-guided assemblies (assembly_method BWA,
   samtools, some CLC) in any scoring set.

## Annotation or reference errors found

None new. B11243 MTLa1 is partial (99 aa), as its record already says.

## Curator rulings (2026-09-26)

- Detection will report `assembly_gap_at_locus` when a long N run sits at the
  flank-anchored locus position.
- Reference-built assemblies are marked and excluded from scoring sets.
- Read-based confirmation of the 10 reference-built isolates' mating type is
  deferred. It is not the current focus. Revisit it when writing the paper:
  map each isolate's SRA reads to both C. auris MTL idiomorphs and test the
  prediction that an N gap over the reference's idiomorph means the isolate
  carries the other one (SCO isolates -> MTLa; England isolates -> MTLalpha).
