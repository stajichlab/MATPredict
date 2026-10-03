# 012. Uncalled C. auris genomes: the MTL is an N gap in reference-built assemblies

- **Category:** assembly-or-annotation artefact
- **Status:** verified (direct tblastn; N-run scan)
- **Lineage:** Ascomycota, Serinales (Candidozyma auris)

## Summary
Twelve of 952 C. auris genomes stay uncalled. In ten, the assembly was built
by mapping reads to a reference, and the reference's idiomorph is an N gap.
The isolate most likely carries the other idiomorph.

## Evidence
- `results/2026-09-26_cauris_uncalled/NOTE.md` (committed `cd08a0e`).
- 8 x GCA_041026{465,485,505,525,545,565,585,605} (SCO_240..279, Portugal,
  PRJNA1140147; CLC): 8.1-8.2 kb of N at the B11220 MTLalpha span.
- GCA_027563775.1, GCA_027563725.1 (England, PRJNA865124; BWA-MEM): 8.8-8.9 kb
  of N at the B8441 MTLa span.
- GCA_902703435.1 (CA1, SPAdes): contig break at the idiomorph.
- GCA_046252445.1: 12 kb of repeats, not a genome (suppressed in BFD).
- Detector report `assembly_gap_at_locus` flags 2 of the 10 (commits
  `7fc397c`, `0e36543`; `results/2026-09-26_gap_zygosity_validation/`), because
  the flanks PAP1/OBP1/PIK1 lie inside the idiomorph in C. auris.

## Method that found it
tblastn (genetic code 12) with the three C. auris records; blastn of a 55 kb
window; N-run scan (`region.py`); assembly methods from NCBI Datasets.

## Verification done / still open
Open: map each isolate's SRA reads to both idiomorphs (deferred until paper
writing, curator ruling 2026-09-27).

## Limits
Low coverage could make the same gap; reads not examined.

## Related
Reference-built assemblies are excluded from scoring sets (curator ruling).
