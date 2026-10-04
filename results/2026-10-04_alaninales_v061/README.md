# Alaninales re-run on v0.6.1 (2026-10-04)

The 20 BFD Ascomycota genomes (genetic code 26) that failed in the v0.6.0 run
(../2026-10-03_ascomycota_v060) because exonerate has no table 26. Code: frozen worktree
run-ff0e581 (branch fix-gencode-26, v0.6.1): exonerate skipped, miniprot-only models.

## Result (job 29393163, 48 s on 16 CPUs)
- 20/20 reports, 0 failures; 30-35 s per genome. Every genome logged
  "exonerate has no genetic code 26; gene models use miniprot only".
- 4 genomes called, 16 uncalled. There is no curated Alaninales record, so these
  are searched against the whole phylum (reference gap, not measured absence).
  - GCA_001661245.1 (Pachysolen tannophilus): MATsc MATa medium, MTL alpha medium,
    Ascomycota:MAT MAT1-1 high.
  - GCA_030556755.1: MATsc MATa medium. GCA_030567755.1: MATsc MATa high.
    GCA_043388405.1: MATsc MATalpha medium.
- Gene models are `polished_single` (miniprot) or `unpolished`; none use exonerate.
