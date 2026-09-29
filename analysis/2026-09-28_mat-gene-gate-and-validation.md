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

## Details (added 2026-09-29)
- Validation (`results/2026-09-28_validation_f3_f4/NOTE.md`): CAAX finder
  6/9 mating receptors vs 0/25 other STE3 copies (95% upper limit 13.7%);
  about 13 false admissions expected by chance among 118 gained calls.
- Gate (076afe4): 17 calls withheld (16 undetermined); 3 gained, including the
  S. racemosum NRRL 2496 record locus. Lichtheimiaceae are called only on the
  score route (no Mucorales flanks); 5 genomes lost their only call, all
  undetermined.
- Later: the gate threshold is now set by each classifier build (99.9 bits on
  PR #9); see the classifier-builds report.

## Curator decisions
Made 2026-09-28: gate option (a); all CAAX-dependent calls unverified, to be
reviewed once >= 100 labelled non-mating STE3 loci exist; Lichtheimiaceae
option (a) for now.
