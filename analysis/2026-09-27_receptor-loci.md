# Pheromone receptor (B/PR) loci
Status: decided (CAAX calls labelled unverified; review later)

## Results
- Receptor tree (`results/2026-09-27_receptor_explore/NOTE.md`): Agaricomycete
  mating receptors do not separate by tree (3–7 STE3 copies per genome);
  Sporidiobolales A1/A2 do.
- Positional rule (`results/2026-09-27_pheromone_positional/NOTE.md`): strict
  CAAX within 10 kb flags 6/9 mating receptors, 0/25 other copies; random
  windows 2.5%.
- Receptor records Heterobasidion, Trametes, Grifola, Russula nobilis
  (`results/2026-09-27_receptor_curation/NOTE.md`,
  `results/2026-09-27_russulaceae_receptor/NOTE.md`).
- CAAX finder c8412b7 (`results/2026-09-27_caax_precursor/NOTE.md`): Agaricales
  PR calls 14->132, 0 lost, runtime x1.08; 62% in curated mating clades (~1.9x).
- Locus merge a2fe24e (`results/2026-09-27_locus_merge/NOTE.md`): Agaricales
  303->268 calls; Schizophyllum B reported once.
- HD prescreen (`results/2026-09-27_hd_prescreen/NOTE.md`): ~10% saving; not built.
- Negative control (`results/2026-09-28_validation_f3_f4/NOTE.md`): ~13 of 118
  gained calls expected by chance; 0/25 non-mating set too small.

## Curator decisions
Made: all CAAX-dependent calls unverified (2993ee1); review once >=100 labelled
non-mating STE3 loci exist. Subloci: one A and one B call with subloci as
evidence.
