# CAAX precursor scan: measured effect (2026-09-27)

Code: polish-scope-cuts c8412b7; basidio-anchors rebased to f017ecb.
Arms: on = run-f017ecb; off = run-f017ecb-nocaax (same code/db, scan removed
from the Basidiomycota:PR roster, uncommitted); full = 2026-09-26 full
Basidiomycota run (run-ad1f865). Jobs 29156032-35 (short). Script: compare.py.

Panels: 128 Agaricales genomes stratified over 45 families (lists/
agaricales_panel.tsv; seed 20260927; includes Coprinopsis CC3/AmutBmut and
Schizophyllum H4-8/14-01) + 35 Polyporales/Russulales + Cryptococcus JEC21 and
Ustilago 521 (lists/polyrussu_controls.tsv). 162 genomes have both arms.
Timeouts (1 h/genome): GCA_021399455.1, GCA_045781145.1 in both arms;
GCA_049639265.1 in the on arm only.

## Result
- Basidiomycota:PR calls off -> on: 14 -> 132 (118 gained, 0 lost, 0 changed).
  Balpha/Bbeta 10 -> 10. HD-type calls: 0 genomes changed. Versus the full run
  (Agaricales): 95 gained, 0 lost, HD identical.
- All 118 gained calls carry a scan precursor (median 2 ORFs, max 8), all are
  medium, all have a polished receptor (97 polished_disagree, 21 agree),
  receptor identity to curated receptors 31.5-100% (median 62.1).
- 0 of 118 gained calls lie on an HD contig.
- Gained calls per genome: 1 in 62 genomes, 2 in 22, 3 in 4.
- Runtime: 21.34 h on vs 19.73 h off (x1.08).

## Known-answer genomes
- Coprinopsis CC3 and AmutBmut: already called (high) without the scan; the
  same call now also carries 4 scan ORFs. Not a gain.
- Schizophyllum H4-8: gained PR call at NW_026089548.1:190,191, the same
  position as its Balpha and Bbeta calls; receptor nearest neighbour = curated
  bar3 (100%). 14-01: same pattern (JAGVSI010000976.1). This is the real B
  locus, but it is now reported twice (PR and Balpha/Bbeta families).
- Russula nobilis (record genome): now called, medium, 1 ORF, OZ475200.1 at
  1,510,837 (record span 1,514,818-1,529,433).
- Cryptococcus, Ustilago: no PR calls (out of PR scope), unchanged.

## Are gained receptors mating receptors?
Test: nearest neighbour (blastp) of each gained receptor in the step-1 STE3
set, and whether it falls in the small subclades (12-60 tips) around the six
curated Agaricales mating receptors (87 tips).
- Agaricales gained calls: 59/95 (62%) have a nearest neighbour inside those
  subclades; 25 hit a curated mating receptor directly.
- Baseline: 66/197 (33.5%) of all Agaricales STE3 copies in the step-1 tree
  sit in those subclades. Enrichment about 1.9x.
- Only 1 gained Agaricales call has a receptor < 50% identity with a nearest
  neighbour outside the subclades (Inocybaceae GCA_043165475.1, 47.7%, 1 ORF).
- The subclade test is not valid for Polyporales/Russulales (subclades are
  Agaricales-defined): 0/15 and 3/8.

## Families to check (possible non-mating STE3 paralogs)
Lower subclade share among gained calls: Agrocybaceae 8/16 (12 genomes),
Mycenaceae 2/6, Physalacriaceae 2/5, Galerinaceae 0/2. Not proof of paralogs:
the subclades cover only six curated receptors from two families.

## Limits
- No independent truth for most gained calls; the subclade test is indirect.
- The off arm is an uncommitted roster edit in a frozen worktree.
- Double reporting of the Schizophyllum B locus (PR + Balpha/Bbeta) is new.
