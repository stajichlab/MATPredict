# Bitscore floor for the flank-carried rule (evaluation, 2026-09-27)

Question: replace the genome-size-dependent E <= 1e-5 floor on a flank-carried
call's strongest core hit (flank_carried.py, b0898d9) with a bitscore floor?

Data: 265 flank-carried loci (no modelled core gene). 176 collected from reports
(collect.py -> loci.tsv): Mucoromycota 4457c3a (13), early-diverging 634dda4 (19),
Serinales 882aa01 (131), Basidiomycota full ad1f865 (11), Dothideomycetes
re-check (2); plus the 89 cap6 Ascomycota calls of the 2026-09-26 audit
(classified3.tsv). Bitscores: fresh genome-wide tblastn with detect settings
(-seg no, evalue 10, report genetic code; hit_one.sh, hits/) for the 159
report genomes; the audit's own tblastn for cap6. Genome lengths recorded for
option (b).

Labels (independent of score): real = 8 unique loci (5 Umbelopsis with the
Mucorales gene order, 2 Mucor irregularis, Trigonopsis variabilis; curated
diagnoses); noise = 52 Serinales second calls outside the PAP1-OBP1-PIK1 span
in genomes with a core-modelled call (Debaryomyces artefact). The audit's
cap6 strong/weak/noise classes are E-based, so only indicative.

See tradeoffs.txt for all tables. Key rows:

| rule | real (unique) | indep. noise | cap6 pass (audit-noise) | Serinales pass |
|---|---|---|---|---|
| E <= 1e-5 (current) | 3/8 | 0/52 | 9/89 (0) | 24/131 |
| E fixed 1e8 <= 1e-5 | 3/8 | 0/52 | 8/89 (0) | 22/131 |
| bits >= 35 | 7/8 | 0/52 | 38/89 (9) | 56/131 |
| bits >= 38 | 7/8 | 0/52 | 27/89 (1) | 36/131 |
| bits >= 39 | 7/8 | 0/52 | 21/89 (0) | 32/131 |
| bits >= 40 | 6/8 | 0/52 | 19/89 (0) | 30/131 |

Independent noise tops out at 30.4 bits; real loci span 33.5-68.2 (the one
below 39 is U. nana, which the classifier leaves undetermined anyway).
Option (b) behaves like the current rule and recovers none of the Umbelopsis.

Recommendation: bits >= 39, one floor for all families. Limits: tuned on 8
real loci; U. vinacea gzUmbVina2 (39.7) passes by 0.7 bits; the extra cap6
passers (12 vs current) are audit-"weak" Dipodascales/Orbiliales calls whose
reality is not established; Basidiomycota HD flank-carried loci all score
>= 42.7 and are unaffected in practice.
