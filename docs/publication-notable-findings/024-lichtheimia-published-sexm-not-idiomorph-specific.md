# 024 — The published Lichtheimia "SexM" is probably a conserved, non-idiomorph HMG gene

- Category: biology; annotation artefact (literature)
- Status: candidate (concerning; curator flagged 2026-10-01)
- Lineage: Lichtheimiaceae sensu lato (Lichtheimia, Circinella group, Rhizomucor, Dichotomocladium, others)

## Summary
Schulz et al. 2016 list a SexM, and no SexP, in three Lichtheimia species; the
L. ramosa protein is CDS03202.1. In our genomes, an ortholog of CDS03202.1 is
present in all seven annotated genera of the family and sits beside the same
conserved neighbours (an Hsp90 gene and two hypothetical genes). In Circinella
it is present in strains whose sex locus types Plus and in strains that type
Minus. A gene found in both mating types is not idiomorph-specific, so
CDS03202.1 is probably not the Lichtheimia sexM.

## Evidence
- Literature: Schulz et al. 2016, Endocytobiosis Cell Res 27:39-57, Table 1
  (local PDF resource/Schulzetal.2016.pdf, pp. 42-43): SexM only in L.
  corymbifera, L. hyalospora, L. ramosa (CDS03202.1); sexM on a different
  scaffold from the tptA/algL/rnhA block (draft assemblies). Summary:
  results/2026-10-01_lichtheimiaceae/LITERATURE.md.
- Classifier: the shipped Mucoromycota classifier scores CDS03202.1 sexM 87.0,
  sexP 72.1, P1 91.5 — below the gate, untyped.
- Orthologs: a full-length CDS03202.1 ortholog in 33/33 annotated Lichtheimia
  genomes (blastp 62-96% to CDS03202.1; 35-40% to Mucorales sexM); annotated
  proteins there score sexM 50-98, none typed.
- Synteny: one ortholog cluster across all 7 annotated genera (73 genomes);
  Hsp90 + two hypothetical genes within 50 kb in all genera (Lichtheimia
  29-31/33, Circinella 20-24/28, Rhizomucor 7/9); no Mucorales flank on its
  contig in any genome.
- Not idiomorph-specific: in Circinella the ortholog is present in genomes with
  Minus-typed and Plus-typed rnhA-adjacent sex loci.
- Source: results/2026-10-01_lichtheimiaceae/NOTE.md (anchors.tsv,
  synteny_clusters.tsv, flank_positions.tsv).

## Method that found it
tblastn + exonerate models scored with the HMM idiomorph classifier
(commit 52b3ff9); DIAMOND all-vs-all clustering of genes within ±100 kb of
anchors; comparison with the published protein.

## Verification and open items
- Not shown for Lichtheimia itself: no Lichtheimia strains of known mating
  type are available, so the "present in both types" test rests on Circinella.
- The true Lichtheimia sex locus has not been found: no sexM/sexP within
  50 kb of rnhA/tptA in any of 51 genomes. Lichtheimia stays discovery-only.
- Family delimitation is unsettled (Hoffmann 2013 vs Walther 2019), so
  "Lichtheimiaceae" in BFD may mix genera with different locus layouts.
- Not added as a paralog class until shown in Lichtheimia of known type.

## Limits
Draft assemblies; short core-gene models (~60-80 aa); single-linkage
clustering may merge paralogs.

## Related
Entries 014 (sexP-type genes in Lichtheimiaceae/Syncephalastraceae),
016 (classifier); analysis/ Lichtheimiaceae report.
