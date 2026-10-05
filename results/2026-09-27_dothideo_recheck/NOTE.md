# Dothideomycetes pilot re-check on current code

2026-09-27. Branch curation-mucor-dothideo rebased onto polish-scope-cuts
4457c3a -> 9932c7c (790 tests pass); frozen run-9932c7c; 99 genomes (the
Tubeufia bin is suppressed). NOT pushed.

Genomes called: d18ec5a 89/100, 52cd292 (round 2) 87/100, 9932c7c 90/99.
Calls vs 52cd292 (compare_vs_52cd292.tsv): 74 same, 3 new (Zymoseptoria
brevis MAT1-1 medium, GCA_019670975.1 MAT1-2 high, Bauco1 MAT1-2 medium),
13 confidence changes (12 medium -> high, 1 high -> medium), 1 class change.

Double calls:
- Nothopassalora personata GCA_059623555.1: GONE. The round-2 MAT1-2 call
  (CM184477.1:308,207-309,146, 940 bp, MAT1-2-1 78.1% + COX13 in the same
  span) is not made; that region now holds only an unadmitted COX13 cluster.
  One call remains: MAT1-1, medium, partial_locus (MAT1-1-1 71.3% modelled,
  APN2 66.7%).
- Cercospora kikuchii GCA_009193115.1 (UBA assembly): REMAINS. Two contigs:
  VTAY01000105.1 MAT1-1 medium partial_locus (MAT1-1-1 64.0% modelled, COX13
  72.6%; weak MAT1-2-1 30.6% ignored by the tier) and VTAY01000089.1 MAT1-2
  medium mat_locus (MAT1-2-1 70.8% modelled, MAT1-1-1 51.2% modelled,
  MAT1-1-3 34.4% unmodelled, APN2 29.0% unmodelled). Both MAT1-1-1 and
  MAT1-2-1 are modelled at distinct loci, and the second locus carries both.
  This is real two-idiomorph evidence in the assembly: a homothallic
  arrangement, a mixed/contaminated assembly, or a paralog. Not resolved.
