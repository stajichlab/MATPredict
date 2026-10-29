# Taphrinomycotina MAT candidates and why Neolecta and Saitoella are not called (2026-10-28)

Two questions: (1) collect every candidate MAT protein of the Taphrinomycotina assemblies (low confidence included) for a gene tree; (2) why the detection pipeline withholds the *Neolecta* and *Saitoella* loci, and would closer representatives, a consensus or an HMM help.

## Candidates for a tree
Names in the tables and FASTA headers are labelled by reference system (`SchPomMc`, `SchPomPi` ...; columns `ref_gene_label`, `idiomorph_label`; the raw `ref_gene` is kept): see `docs/notes/2026-10-28_idiomorph-naming-proposal.md` and `scripts/idiomorph_labels.py`.
`scripts/taph_mat_candidates.py` aligns 66 reference MAT proteins (13 Taphrinomycotina P/M proteins, 41 alpha-box/HMG-box proteins of other Ascomycota records, one per order and gene name, and 12 Mucoromycota sexP/sexM) to the 104 Taphrinomycotina assemblies of the BFD table with miniprot and translates every alignment.
- `candidates.faa` / `candidates.tsv` (identity >= 25%, query coverage >= 40%): 466 proteins in 93 genomes, 27 species; 79 distinct sequences, 50 species x gene groups. 462 are P/M-type (matPi 116, mat1-Mc 112, mat1-Mi 99, matPc 96, matMc 24 ...), 4 alpha/a-type (MAT1-2-1 in *Protomyces*). 394 are *S. pombe*.
- `loose/` (>= 20% and >= 20%): 468 proteins; adds a *Saitoella coloradoensis* matMc-like protein (215 aa, 40%) and one MAT1-2-1.
- `domain/domain_hits.*`: tblastn of the same references against every locus window of the v0.7.0 detection reports (detected and withheld), one aligned peptide per hit of at least 33 bits (the pipeline's floor): 1,651 hits. Domain-level, mostly HMG-box fragments; many are not MAT genes.
- Each candidate row says whether it lies in a locus the reports list as detected or withheld. Nothing is filtered by confidence. 11 genomes have no candidate (all three *S. japonicus*, *S. complicata* ...).

## Why the loci are withheld
The pipeline needs two modelled genes in one cluster (`MIN_POLISHED_GENES = 2`).
- *Neolecta irregularis* (1,039 contigs, N50 21 kb): v0.7.0 finds matMc on contig `LXFE01003167.1` (5-851, contig edge) and matPi on `LXFE01003595.1` (11,513-12,163), different contigs; the 19 kb window `LXFE01000411.1:5949-25065` with both gene names gets no model from any tool or query set and an HMM score of E = 3e-13 / 1e-6, so it is probably not a PM locus.
- *Saitoella* (both assemblies): one window each with PM-gene hits (S. coloradoensis `AP039747.1:177119-198616`, S. complicata `NW_017567017.1:185699-210006`), 1 modelled gene (matMc), no matPi found in it.

## Probe: does a closer query help? (`probe/probe_table.tsv`, `scripts/taph_polish_probe.py`)
11 windows (3 *Neolecta*, 3 *Saitoella*, 5 *Protomyces*/*Taphrina* controls); per gene (Pi, Mc, Mi), exonerate and miniprot with (A) the pipeline's own references, (B) A plus close representatives from other genera (candidates, identity >= 35%), (C) a consensus of the alignment; a profile HMM from the alignment scanned over the six-frame translation.
- matMc: B beats A where A is weak. *S. coloradoensis* and *S. complicata* windows: A gives a 0.43-identity model covering 40% of the query, B a 0.66-0.68 model covering 99-100%. *Protomyces* control: 0.33 identity covering 62% -> 0.45 covering 100%. *Neolecta* matMc: 0.41 -> 0.46-0.53.
- A consensus is not better than B: partial coverage, and no miniprot model in several windows.
- HMM: E = 1e-30 to 1e-73 where a model exists (matMc: *Neolecta* 4e-73, *S. coloradoensis* 3e-65 / 8e-32); no hit in windows without a gene. It separates real P/M genes from the HMG-box background better than tblastn at 33 bits.
- *Neolecta* matPi (`LXFE01003595.1`): exonerate models it at 0.36-0.39 (miniprot does not); HMM E = 2e-44.
- Saitoella windows: no matPi by any query set or the HMM.

## End to end: detect with 31 close proteins added (`e2e_summary.txt`, `scripts/taph_e2e_test.sh`)
`detect --phylum Ascomycota` on the three genomes with the release database and with a copy whose Taphrina PM record had 31 close matMc/matPi proteins (other genera, >= 40% identity) appended to `proteins.faa`:
- *Neolecta*: modelled genes in the best cluster 0 -> 1 (matPi on `LXFE01003595.1`). Still not called.
- *Saitoella* (both): 1 -> 1. Not called.
- Suppressed PM loci rose from 23-27 to 66-78 per genome: the added references pull in many more HMG-box hits. Naive extra references cost precision and gain no call.

## Reading
1. Polishing is helped by closer representatives (matMc identity and coverage), and an HMG/profile HMM is sensitive and specific for these genes. A consensus query is not.
2. The 2-gene bar is what keeps *Neolecta* and *Saitoella* uncalled: in *Neolecta* the two genes are on different contigs of a fragmented assembly; in *Saitoella* the second gene is absent from the window.
3. Unexplained: the pipeline records 0 polished genes on the *Neolecta* matMc contig although the probe models it (0.41-0.53) with the same reference set. The pipeline's window or acceptance rule differs from the probe's; not traced.
4. Not done: curated Taphrinomycotina records for *Saitoella*/*Protomyces*/*Neolecta* (a proper candidate record, not appended proteins), a per-gene HMM as a MAT-gene-gate profile, and a PM-specific single-gene admission.
