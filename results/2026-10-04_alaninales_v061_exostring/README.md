# Alaninales re-run, exonerate given the table-26 string (run-b04a77b, 2026-10-04)
Same 20 genomes as ../2026-10-04_alaninales_v061 (miniprot-only, run-ff0e581). Job 29394962.

## Result (2026-10-04)
- 20/20 reports, 0 failures. Every genome logged "exonerate has no built-in
  genetic code 26; passing the NCBI table as a string"; no exonerate error.
- Called 15 of 20 genomes (miniprot-only run: 4 of 20). The 4 earlier calls are
  kept with the same families and idiomorphs; P. tannophilus GCA_001661245.1 adds
  a MATyl A partial_locus. 11 genomes gained a call (mostly one MATsc locus).
- Cause: with exonerate, genes that miniprot left unmodelled now have a model
  (`polished_single`/`polished_disagree`), so the modelled-gene bar passes.
- `polished_disagree` (exonerate and miniprot boundaries differ) is common: 21
  genes in the 15 called genomes (7 in P. tannophilus).
- Cost: median 457 s per genome (miniprot-only: 34 s).
- Caveat: all route `phylum_fallback` (no Alaninales record). Family labels are
  not reliable there: GCA_003706035.2 and GCA_003706045.2 get MATsc MATalpha and
  MTL A in the same genome.
