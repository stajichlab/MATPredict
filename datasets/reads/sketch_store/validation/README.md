# Validation of the sketch-store settings (2026-11-07)

The test that chose the sketch store's settings, on the UCR HPCC. Inputs and outputs are in
`/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-11-07_sketch_pcoa_mat/` (the scripts ran with copies of themselves there):
`strains.tsv` (136 Afumigatus_Global strains: 18 lab-strain copies, 78 spread over the source collections, 40 with Starship truth),
`targets/` (bafA and 20 Starship haplotypes from the NII Afumigatus_pangenome study, with `truth_overlap.tsv`), `analysis_solid/`,
`analysis_plain/`, `analysis_followup/` (REPORT.txt, pcoa.png).

| script | what it does |
|---|---|
| `sketch_calib.sh` | size calibration: Af293 genome and one 2M-spot subsample (SRR13579397) at scaled 100 and 1000, k21/k31, with and without the abundance >=2 filter (results/2026-11-07_sketch_calibration) |
| `sketch_test_task.sh` | array task per strain: two independent 30x subsamples of the local FASTQs (pairs 1..N and N+1..2N), each sketched (k21,k31, scaled=100, abund; kept raw and filtered >=2), subsample A typed with reads-type |
| `snp_gt.sh` | genotype matrix of the test strains from the Afumigatus_Global v5 SNP VCF (biallelic, variable among them) |
| `sketch_analyze.py` | same-strain ceiling (A vs B), distances vs SNP counts at scaled 100/500/1000 and k21/k31, lab-strain groups, PCoA vs SNP PCoA, MAT from idiomorph-specific hashes vs reads-type, target fractions vs truth; MODE=solid or plain |
| `sketch_followup.py` | outliers, PCoA without them (Procrustes, nearest neighbours), targets on element-specific hashes, bafA disagreements |

Findings (summarised in ../README.md): plain containment cannot separate *A. fumigatus* strains (noise floor ~4e-4, worth ~12,000
SNPs); solid containment can (same-strain median 0.99996, between-strain p99 0.999); PCoA agreed with the SNP PCoA (Procrustes 0.63 on
2 axes, 0.70 on 3) but fine-scale distances did not (rank correlation 0.34); MAT 129/129 at k21; element-specific Starships near
perfect, the Janus family and Osiris-h3 not; bafA 17/19, with 4 of the 5 disagreements traced by read depth to the truth proxy or the sample.
