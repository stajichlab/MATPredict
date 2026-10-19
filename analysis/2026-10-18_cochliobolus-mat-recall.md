# Why the Cochliobolus clade is almost never called (2026-10-18)

Analysis only: `src/` and `db/` are untouched. Branch `cochliobolus-mat-recall`, data and scripts in `results/2026-10-18_cochliobolus_recall/`.
Base: origin/main 38fad7f (includes #109, #113, #114). The audit's reports (#114) were produced on c9c0ec6; `git diff c9c0ec6..38fad7f` touches no `src/` file
and, for the Ascomycota:MAT family, no db record or roster line (the changes add a new MATkp family and Pichiales records). Re-runs for this note used the 0a76d9b
tree for the clade runs and 38fad7f for a confirmation run of the four C4/C5 genomes (identical result, `main_c45_sgr.tsv`). Pixi env `test`, HPCC, Slurm.
"Measured" = a number from a run or table; "judgement" is labelled.

## 1. Short answer

**The MAT gene is found in every one of the 101 uncalled clade genomes and modelled correctly; the call is lost at the evidence floor because Pleosporaceae MAT
is a single-gene idiomorph with no roster flank within 25 kb.** Measured:

* In the four C4/C5 genomes the curated gene is present at 97.9-100% identity over the full length, the exonerate and miniprot models are complete
  (343 aa MAT1-2-1, 383 aa MAT1-1-1, intron 52-55 bp, no internal stop), the Pfam domain is present. The cluster holds one gene, so `_families_meeting_evidence_floor`
  (min 2 distinct genes) never admits it. With `--min-hits 1` it is admitted, polished (1 gene) and withheld `modelled_gene_bar` (gate 2). Same chain as the
  Parastagonospora note, with the strong-single-gene rule as the intended answer.
* SLA2/APN2/COX13 are 180-290 kb from MAT in Bipolaris/Curvularia/Exserohilum (0 of 111 gene rows within the 25 kb clustering window), so no cluster ever holds a
  flank. The genes that do sit beside MAT are **GAP1 and BGL1** (median 2.6 and 3.6 kb; both on the MAT contig within 25 kb in 92 of 98 Bipolaris/Curvularia/Exserohilum gene rows). Both
  are named in the C4/C5 record notes as "not roster genes".
* **The audit's lead was a paralog.** The `MAT1-1-1 + SLA2` candidate in the C5 report (KB445578.1:1,263,708-1,294,483; JARAAY010000009.1:1,432,303-1,463,093 in 4m) is a 37.9% hit.
  All 35 admitted clusters (the 35 suppressed loci) have best identity <= 50%; the true locus is not among them (`trace4.tsv`).

## 2. The four C4/C5 genomes

BioSample/assembly strain (E-utilities): GCF_028858645.1 CocheC4_6 SAMN33211508 "Cochliobolus heterostrophus C4 Improved Draft" strain C4; GCA_028891585.1 CocheC5_4m
SAMN33211509 strain C5; GCF_000354255.1 CocheC4_1 SAMN00113370 strain C4 (ATCC 48331); GCA_000338975.1 CocheC5_3 SAMN00006316 strain C5. No BioSample states the mating type;
the strain name and the record identity are the evidence. Records: 5016_c5_MAT_MAT1-1 (MAT1-1-1 only, AF029913.1) and 5016_c4_MAT_MAT1-2 (MAT1-2-1 only, AF027687.1).

| genome | N50 | true cluster | identity | genes in cluster | admitted? | flank clusters on the contig | admitted clusters (max id) | gate |
|---|---|---|---|---|---|---|---|---|
| C5_4m | 2.12 Mb | JARAAY010000009.1:1,517,960-1,519,217 | 100.0 | 1 | no | SLA2 89.3% 270 kb away, 75.6% 180 kb | 35 (50%) | evidence floor |
| C5_3 | 1.84 Mb | KB445578.1:1,349,345-1,350,602 | 100.0 | 1 | no | 270 kb, 180 kb | 35 (50%) | evidence floor |
| C4_6 | 2.06 Mb | NW_027093865.1:406,610-407,693 | 97.9 | 2 (2nd hit superseded) | no | 270 kb, 180 kb | 35 (50%) | evidence floor |
| C4_1 | 0.96 Mb | NW_014024937.1:306,993-308,076 | 97.9 | 2 (superseded) | no | 270 kb, 180 kb | 42 (50%) | evidence floor |

Independent confirmation (`confirm4/`, `confirm_locus.py`): blastn of the record's locus.gbk 99.9-100% over 1,087-1,202 bp (one hit), tblastn, miniprot (`-G 5000` and default, identical),
exonerate protein2genome (`--maxintron 5000`): every model complete, starts with M, 0 internal stops; model vs record 100% (C5) and 99.7% (C4, one residue). Pfam (GA): MATalpha_HMGbox
PF04769 225.5 bits on the C5 model; HMG_box PF00505 84.7 and HMG_box_2 37.9 on the C4 model. C5 genomes have no C4 gene and vice versa, as expected for heterothallic
strains. Gene modelling and the 5 kb intron cap are **not** involved. Not run: reads, a mapped-contig check of the locus orientation.

Runs on the four genomes: baseline no call; `--min-hits 1` no call (170-182 suppressed, true locus `modelled_gene_bar`); `--strong-single-gene-policy on`: MAT1-1 (C5) and
MAT1-2 (C4) `idiomorph_gene_only`, low, model identity 100/99.71, coverage 100 (`main_c45_sgr.tsv`).

## 3. The whole clade (103 genomes; 101 uncalled)

Located independently of the pipeline with the audit's six-frame ORF scan (alpha-box profile >= 100 bits = MAT1-1-1; HMG profile >= 80 bits = MAT1-2-1; the HMG gap between
paralogs, <= 70, and real genes, >= 86, was read off the distribution, so the 80 bit threshold is chosen from these data) and a tblastn map of the roster genes.

| class (101 uncalled) | Bipolaris 57 | Curvularia 23 | Exserohilum 10 | Ophiosphaerella 11 |
|---|---|---|---|---|
| MAT-box ORF present in the genome | 57 | 23 | 10 | 11 |
| cause: evidence floor (lone gene, cluster not admitted) | 57 | 16 | 5 | 11 |
| cause: `modelled_gene_bar` (gene in an admitted cluster with APN2 or a paralog, 1 polished) | 0 | 7 | 5 | 0 |
| no hit / no recoverable locus | 0 | 0 | 0 | 0 |
| GAP1/BGL1 beside the MAT gene (<= 25 kb) | 53 (4 mixed: fragmented gigantea) | 23 | 10 | 0 (MAT contig lacks them) |
| N50 < 100 kb / 100 kb-1 Mb / >= 1 Mb | 9 / 14 / 34 | 1 / 17 / 5 | 0 / 2 / 8 | 3 / 8 / 0 |

So 101 of 101 carry the locus and are missed by a gate; none lacks a locus by this test (limit: the scan shows a box ORF, not that the locus is intact or functional).
Fragmentation matters little: 4 B. gigantea genomes (N50 14-45 kb) have the MAT gene on a short contig apart from GAP1/BGL1. Six genomes (five B. gigantea, B. maydis KET7) carry both an
alpha-box ORF and an HMG ORF on different contigs (`clade_gates.tsv`); with fragmented assemblies this cannot separate a homothallic arrangement from contamination or a mixed
isolate (open question). The audit's two clade calls (B. maydis KET7, C. kusanoi 30M1) are unchanged by anything here.
Ophiosphaerella (Phaeosphaeriaceae, N50 <= 0.47 Mb) differs: GAP1/BGL1 are present at 66-87% but not beside MAT, as in Parastagonospora; its flank is SLA2/APN2 on other contigs.

## 4. Remedies (all analysis only; none implemented)

| option | genomes gained of 101 | confidence | FP risk | het. / hom. |
|---|---|---|---|---|
| (a) gene-modelling fix | 0 | n/a | n/a | the model is already complete; there is nothing to fix |
| (b) relax modelled-gene bar when a flank anchors | at most 5 (Exserohilum clusters holding APN2; counted from reports, not run). Curvularia's 7 held clusters pair the MAT gene with a paralog, not a flank | low | the 55% paralog calls of the Parastagonospora note | does nothing for the 90+ lone genes: SLA2 is 270 kb away, outside the 25 kb cluster window |
| (c) strong single gene on (85% id, 90% cov) | **46** (Bipolaris 46 of 57; Curvularia 0, Exserohilum 0, Ophiosphaerella 0). 48 genomes called in all, 0 idiomorph disagreeing with the scan | low by construction | none seen vs the scan; the control in the single-gene note (600 called, 0 changed) was not redone | recovers Bipolaris only because the records are B. maydis (Curvularia/Exserohilum hits are 70-78% to them, under the gate); 5 B. gigantea genomes with both ORFs get one call only |
| (d') GAP1+BGL1 flanks in the roster | **90 of 101** (92 of 103 genomes called, 2 already called: Bipolaris 58/58, Curvularia 24/24, Exserohilum 10/10, Ophiosphaerella 0/11), 12 `high`, 81 `medium`, all idiomorphs agree with the scan | medium-high | **global roster addition is not acceptable**: on 94 called control genomes (Dothideomycetes, seed 7) it lowered 24 from `high` to `medium`, lost a second locus in 1 and added one in 1 (`control_compare.tsv`) because the roster denominator grew. Seven of the new `high` clade loci list a spurious `MAT1-2-1` or `MAT1-1-1` cross-hit as the second core gene | needs both genes in one cluster, so het. fine; a homothallic fused locus would show both genes beside GAP1/BGL1 (not tested here, no such genome identified) |
| (d) more curated records | gains none alone: a lone-gene record does not change the floor | | | |

Judgement. Prototype (d') was a throwaway copy of `db/` (order.yml roster plus GAP1/BGL1 in the two Bipolaris records' proteins.faa), not tuned on the clade except that the 98 rows
were also the test set. Recommended: **(c) turned on now for what it does well, plus a scoped version of (d')**: add GAP1/BGL1 (and ORF1) as flanking genes of
Pleosporaceae records only (a per-order or per-record roster so other families' fractions are untouched; the global form costs 24 of 94 controls their `high`),
and add Curvularia and Exserohilum single-gene records so that (c) also reaches them. The curator ruled (c) off, so the order is: (d') scoped
first (gains 90, medium/high, no ruling needed), (c) as the fallback for lone genes it cannot anchor (Ophiosphaerella).

Candidate reference genomes (N50, from this table; none has a BioSample mating-type statement, so MAT identity rests on the scan and, for B. maydis, the strain/literature):
Curvularia clavata yc1106 GCA_023920165.1 (N50 4.75 Mb, alpha-box), Exserohilum turcicum Et52B GCA_050059805.1 (2.45 Mb, HMG), E. rostratum LWI GCA_026367775.1 (2.2 Mb, HMG),
B. sorokiniana WAI2406 GCA_008452705.1 (chromosome-level, HMG) and ND90Pr GCF_000338995.1, B. oryzae Bo-1 GCA_037177745.1 (2.4 Mb). Ophiosphaerella has nothing above 0.47 Mb. Homothallic Cochliobolus
species are known from Yun et al. 1999 (PubMed 10318929, [doi](https://doi.org/10.1073/pnas.96.10.5592)): homothallic species carry both MAT genes closely linked or fused, heterothallic MAT organisation is
conserved. The abstract names no species and I could not tie a clade genome to a homothallic species; Ulocladium botrytis carries both genes (PubMed 28801599, [doi](https://doi.org/10.1038/s41598-017-08471-3)).
Cited literature only what PubMed returned: 10318929, 28801599.

## 5. Limits and open questions

* Truth is the pipeline-independent scan, not a stated mating type; no clade BioSample states MAT. The 101-of-101 figure depends on the 100/80-bit thresholds.
* (d') was run on the 103 clade genomes and 94 Dothideomycetes controls only; other phyla were not run, and the effect on other families is not known.
* (b) is counted, not run. (c) was not re-run on controls.
* Do the five B. gigantea genomes carrying two box ORFs on separate contigs represent a homothallic arrangement? Reads or a long-read assembly are needed.
* The seven spurious cross-hit second genes in (d') `high` calls should be checked before any confidence ruling.
* Should GAP1/BGL1 live in a Pleosporaceae-scoped roster, or should Pleosporales get its own locus family?

Files: `clade_gates.tsv`, `cause_by_genome.tsv`, `clade_flank.tsv`, `trace4.tsv`, `sgr_clade.tsv` + `sgr_clade_eval.tsv`, `proto_clade*.tsv`, `proto_control.tsv`, `control_compare.tsv`,
`confirm4/`, `audit_run_c45/`, scripts and Slurm files alongside.
