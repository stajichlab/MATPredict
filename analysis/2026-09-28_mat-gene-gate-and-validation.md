# MAT-gene gate, F6 and validation of F3/F4
Status: decided

## Results
- Validation (`results/2026-09-28_validation_f3_f4/NOTE.md`): classifier typing
  0/540 wrong; margin does not separate MAT genes from paralogs; full-protein
  score >=100 bits separates 96/108 vs 9/189.
- Gate 801c55f: model-typed >=100 bits or >=2 flanks modelled at >=40%;
  fragments need flanks. Mucoromycota 258->253 genomes; paralog pairs removed;
  S. racemosum record locus called (`results/2026-09-28_next_fixes/NOTE.md`).
- F6 e41b08b: sub-floor unpolished genes no longer cap confidence (12 calls up).
- All CAAX-dependent calls unverified: 29->102 in the Agaricales panel.
- Assembly accession backfilled for 22/115 records.
