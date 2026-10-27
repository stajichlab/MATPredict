# Neurospora crassa reads-type panel

`neurospora_MAT1-1.fasta` (mat A) and `neurospora_MAT1-2.fasta` (mat a), each the whole MAT locus with its flanks, cut from two BFD assemblies on the coordinates `matpredict detect` reported (v0.7.0 run, both `mat_locus`, high confidence):

| Idiomorph | Assembly | Contig | Span (1-based) | Length | Genes found |
|---|---|---|---|---|---|
| MAT1-1 (mat A) | GCF_000182925.2_NC12 (the reference genome) | NC_026501.1 | 1,851,649-1,872,342 | 20,694 bp | MAT1-1-1, MAT1-1-2, MAT1-1-3, APN2, SLA2, COX13 |
| MAT1-2 (mat a) | GCA_000786625.1_Neucr_trp3_1 | KN389756.1 | 58,545-77,356 | 18,812 bp | MAT1-2-1, APN2, SLA2, COX13 |

Why not the curated records: the two database records for this species (M33876.1 mat A, 5,306 bp; M54787.1 mat a, 3,235 bp) are idiomorph-only deposits with no shared flank, so the single-copy depth control would be zero (the case that left every *C. posadasii* strain `low_depth`). These panels carry about 11 kb of shared flank.

K-mer content (k = 31, `--min-unique-run 100`, `matpredict.reads.panel.Panel`): 7,017 idiomorph-specific k-mers for MAT1-1, 5,344 for MAT1-2, 11,219 shared k-mers (the flank control).

Not yet checked: how the panel types real reads. The first check is the 20-run throughput test of the pilot (`datasets/reads/sra_survey/PILOT_PLAN.md`), which includes strains with an assembly or a mating type in the strain name (FGSC 2489 mat A, 74-ORS-6a mat a). The loci are cut from one strain per idiomorph, so a divergent wild isolate may lose some unique k-mers; if so, the panel is trimmed to the idiomorph core, as for `cocci2`.
