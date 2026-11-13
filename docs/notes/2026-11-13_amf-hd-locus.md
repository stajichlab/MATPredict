# AMF HD1-like / HD2 locus: alleles, controls and a contamination check
Status: in progress (homeodomain tree still running when this note was written)

## Question
Is the Rhizophagus HD1-like/HD2 locus of Ropars et al. 2016 found in other Glomeromycota, how variable are its alleles, and what data exist to genotype it and to measure its expression?

## What the paper says (full text read; resource/nmicrobiol201633.pdf, doi 10.1038/nmicrobiol.2016.33)
- Locus found as a region at half coverage in two heterokaryon-like isolates (A4, A5); PCR on sorted single nuclei supports one allele per nucleus. A third heterokaryon-like isolate, DAOM-240409, was found by Sanger sequencing of 24 isolates.
- Two genes in opposite orientation. HD2 has a Pfam homeodomain hit (PF00046, E = 9.4e-13); the authors could not find a homeodomain motif in HD1-like by manual inspection. Flanks (Fig. 2): choline transporter, P-glycerate mutase, hypothetical protein.
- Amino-acid identity between alleles 44% to 85%; 98.5% for surrounding genes. At least six putative mating types. Allele phylogeny is not congruent with the species phylogeny (isolate MUCL-43196 carries an allele like MAT6 of some R. irregularis isolates).
- The authors state the involvement of the locus in sex has not been shown functionally.

## Results (this repository)
1. Two HD1-like + HD2 loci on separate scaffolds in A4, A5 and DAOM-240409, as in the paper, and in CHC3 (not in the paper). One allele per scaffold, in the same order and orientation. Source: stage-1 candidates.
2. Allele identity (mafft whole-protein alignment, identity over aligned residues, results in hd_alleles/):

| Pairs | HD2 median (n, min) | HD1-like median (n, min) |
|---|---|---|
| AMF same species | 0.830 (303, 0.808) | 0.743 (328, 0.458) |
| AMF same genus | 0.806 (174) | 0.757 (177) |
| AMF different genus | 0.314 (954) | 0.347 (576) |
| Basidiomycota curated, same species | 0.898 (1 pair) | 0.707 (1 pair) |
| Basidiomycota curated, same genus | 0.424 (9) | 0.377 (9) |
| Basidiomycota curated, different genus | 0.207 (35) | 0.193 (35) |

   The AMF same-species range agrees with the paper (44% to 85%). The curated database holds one same-species Basidiomycota pair per gene, so it cannot support a within-species comparison with mushrooms. Literature values for Coprinopsis or Schizophyllum alleles are needed and I have not fetched them.
3. Allele groups at >= 96% identity (single linkage): HD2 30 groups across Glomeromycota, 7 containing R. irregularis; HD1-like 17 groups (one group of 23 Rhizophagus genomes plus Oehlia). This is not the paper method and I make no claim about the number of mating types.
4. Contamination check. Oehlia diaphana SD1: HD1-like, HD2 and the choline-like gene are 100% identical to Rhizophagus sp. MUCL 43196; on the locus contig (QZLH01001367.1, 5 genes) 4 of 5 genes are >= 95% identical to MUCL 43196. I treat this locus as not independent and drop Oehlia from the candidate list. The other 8 non-Rhizophagus genomes with an adjacent pair (Ambispora x2, Diversispora x2, Entrophospora sp. SA101, Geosiphon, Glomus x2) have 0 or 1 genes >= 95% identical to R. irregularis DAOM-197198 on the locus contig (medians 38% to 73%). Script: scripts/contig_origin_check.py.
5. HD1-like HMM (hd1_hmm/AMF_HD1like.hmm; 43 unique training proteins, 314 alignment columns). 52 hits with E <= 1e-10 in 45 of 59 Glomeromycota genomes; 44 sit within 5 genes of a Pfam Homeodomain gene. Hits in other lineages are all weak (best E 5e-8, 16 of 83 Chytridiomycota genomes, 5 of 14 Blastocladiomycota, 1 of 18 Basidiobolomycota, 1 of 9 Entomophthoromycota) and are not taken as HD1-like. The HMM finds none of the 11 curated Basidiomycota HD1 (or HD2) proteins at E <= 10. No sequence link between AMF HD1-like and Basidiomycota HD1 is shown by this profile.

## Data for genotyping and expression (ENA metadata, queried 2026-11-13; not downloaded)
- WGS: 148 Rhizophagus runs from 22 isolates (largest DAOM197198 193 Gb Illumina; HR1, C3, 4401 also have Illumina). Most isolates are the assembled isolates themselves, so this is not a population sample. PopSeq (catalog of 10 species) has no AMF species.
- RNA-Seq: 644 runs, 34 taxa; Rhizophagus irregularis 382 + 122 (DAOM 181602/197198) + 29 (197198w). BFD samples.rnaseq_sra.csv lists no AMF run.
- Plan, not run: map reads to the locus, with each allele as a reference, to count alleles per isolate (half coverage and allele ratio, as in the paper); quantify the four locus genes from RNA-Seq.

## Limits
- Pairwise identity uses funannotate gene models; allele protein sequences were not checked against the Ropars GenBank proteins.
- Single-linkage groups chain through intermediates.
- The tree analysis (scripts/homeodomain_tree.py) is not reported here.

## Files
scripts/hd_allele_identity.py, hd_flank_identity_check.py, contig_origin_check.py, hd1_hmm.py; results/2026-11-13_flank_first_scan/hd_alleles/, results/2026-11-13_flank_first_scan/hd1_hmm/, results/2026-11-13_flank_first_scan/amf_reads/ (ENA tables). Jobs 32068537, 32068576, 32068580, 32068585, 32068593.
