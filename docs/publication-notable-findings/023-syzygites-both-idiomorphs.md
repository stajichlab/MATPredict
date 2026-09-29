# 023 — Both Syzygites genomes carry both mating-type genes (sexP and sexM)

- Category: biology; homothallism candidate
- Status: candidate
- Lineage: Mucorales, Syzygites (reported homothallic)

## Summary
Both Syzygites genomes in the ZyGoLife LCG set carry a sexP locus and a sexM
locus on different contigs. The sexM protein is 100% identical between the two
genomes. This fits the two-locus arrangement Idnurm 2011 described for S.
megalocarpus. It is a candidate: fusion or heterokaryosis cannot be excluded.

## Evidence
- Syzygites sp. MES 3091 (LCG): Plus and Minus both called, unlinked (different
  contigs), both medium, both classifier-typed on gene models; margins 226.4 and
  160.4. Flanks: glrA at the Plus locus; rnhA and btbA at the Minus locus. Contig
  GC difference 1.87 points. Source: results/2026-09-29_two_idiomorphs/per_genome.tsv.
- Syzygites megalocarpus SC16 (LCG): Plus called. Minus locus found on
  scaffold_38:1-9448 (sexM 100% identity, classifier score 228.3, margin 160.4,
  with rnhA and btbA) but withheld at the 0.5 fraction floor because tptA, algA
  and glrA are not at that locus. Source: results/2026-09-29_sexM_like_paralog/NOTE.md.
- The SC16 sexM is 100% identical to the MES 3091 Minus call (same source).
- Literature: Idnurm 2011 (PMID 21908600, doi:10.1128/EC.05149-11) reports two
  separate sex loci in S. megalocarpus, each with rnhA/glrA copies, one flank
  copy pseudogenised; Schulz et al. 2016 (Endocytobiosis Cell Res 27:39-57)
  infers different chromosomes. Summary: results/2026-09-29_mucoro_homothallism_literature/NOTE.md.

## Method that found it
LCG held-out scan with PR #9 code 076afe4 (--phylum Mucoromycota, no taxid);
HMM idiomorph classifier on modelled proteins; the neutral two_idiomorphs
statement (commit bcd1e0d); direct tblastn and protein comparison in the
sexM-like paralog analysis.

## Verification and open items
- Not verified: single-nucleus state. Syzygites spores are multinucleate and no
  homokaryon was made (Idnurm 2011). Plus/Minus protoplast fusion yields
  homothallic strains in Absidia glauca (Schulz 2016), so fusion cannot be
  excluded from sequence alone.
- Open: the SC16 Minus locus is withheld by the fraction floor; a strong-core
  floor rescue was measured but not adopted (it adds 28 LCG loci, 21 as a second
  idiomorph) pending the two_idiomorphs checks.
- Open: pseudogenised flank copies were not tested (two_idiomorphs lists them
  under not_assessed).

## Limits
Two genomes only; LCG assemblies from an older pipeline; the two Syzygites
genomes may not be independent isolates.

## Related
Entries 004-007 (homothallism candidates); analysis/2026-09-29 reports on
two_idiomorphs and the sexM-like paralog.
