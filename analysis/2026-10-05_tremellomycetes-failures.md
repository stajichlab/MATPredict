# Tremellomycetes failures: rearrangement or lost sensitivity
Status: open (curator review pending). Measurement only; no code or database change.

## Question
Are the mating-type failures in Tremellales and the sister orders Trichosporonales, Filobasidiales and Cystofilobasidiales caused by genome rearrangements (fused or split HD and PR loci, inversions, loci split across contigs), or by the tool lacking sensitivity (distant references, genes not modelled, gates withholding real loci)?

## Data
- v0.6.0 campaign: `results/2026-10-03_basidiomycota_v060/genomes.tsv`, `loci.tsv`, `reports_all.tar.zst` (detection_report.yaml for every genome). 552 genomes in the four orders.
- Assembly quality from BFD `busco_genome.parquet` (complete_pct) and `asm_stats.parquet` (N50_bp, contig_count): `quality.tsv`. Good assembly means BUSCO >= 70, N50 >= 20 kb and contigs <= 5000. 541 of 552 have a BUSCO value; the 11 without one count as not good.
- Failure means uncalled, or called with only a PR call (all PR-only calls are `unverified`, admitted through a strict-CAAX precursor).
- Pipeline-independent scan of all 552 genomes (`scan_ste3_hd.py`; 551 finished, one genome file missing from the library, see Limits).

| Order | Genomes | Good assembly | Uncalled | Uncalled, good | PR-only calls | Other calls |
|---|---|---|---|---|---|---|
| Tremellales | 342 | 327 | 15 | 10 | 0 | 327 |
| Trichosporonales | 127 | 123 | 67 | 63 | 50 | 10 |
| Filobasidiales | 41 | 40 | 22 | 21 | 11 | 8 |
| Cystofilobasidiales | 42 | 40 | 21 | 20 | 14 | 7 |

(Good-assembly counts for the PR-only and other-call groups: Trichosporonales 50 and 10; Filobasidiales 11 and 8; Cystofilobasidiales 14 and 6, with 1 poor assembly among the 7 others. All 50 and 14 PR-only genomes are good.)

## Method
1. `01_report_summary.py`: for each genome, status, calls, and from detection_report.yaml the suppressed_loci (withheld_reason, genes_found, polished_genes) and not_detected entries. Output `report_summary.tsv`.
2. `scan_ste3_hd.py` (uses `scan_genome.py`, `ste3_all.faa`, `hd_queries.faa`): no pipeline code involved. STE3-like loci by miniprot of 1,095 STE3 proteins (the receptor-study set) and by pyhmmer of Pfam PF02076 on a six-frame stop-to-stop translation. Homeodomain loci by Pfam PF05920 (Homeobox_KN, the HD1 class) and PF00046 (Homeobox), and by miniprot of 86 curated HD-class proteins (db/Basidiomycota HD1, HD2, SXI1, SXI2, bE, bW, Y, Z plus Phaffia rhodozyma HD1 and HD2 from KU315762, KU315770, KU315771, KU315777, KU315779). Strict-CAAX short ORFs (motif C[VI][IV][AVMG]) as in the receptor test. `03_scan_summary.py` merges loci (3 kb) and measures the HD to STE3 distance on one contig.
3. `04_join.py`, `06_classify.py`: join and classify failures (rules in the script header). HD1-class gene means a Pfam PF05920 (Homeobox_KN) locus or a miniprot hit with at least 45% identity to a curated HD. A homeodomain (PF00046) locus near STE3 is counted separately, because it is where an HD2-type allele would fall; PF00046 elsewhere in the genome is not counted (it also matches ordinary homeobox genes).
4. `09_withheld_pairs.py`: contents of withheld clusters.
4b. `10_followup.py`: withheld-cluster medians and an HD allele-type cross-tab (`followup_summary.txt`).
5. `07_inorder_ref_test.py`, `08_inorder_eval.py`: direct test of in-order references (protein fragments from called genomes of the same order aligned with miniprot).
6. A re-run on current main (60a3cd0, worktree run-tremello) with `--evidence-diagnostics` was submitted (jobs 29428627 to 29428629). It had not finished when this note was written, so the note uses the v0.6.0 reports only (see Limits).

## Results
### 1. Where the arrangement is, called against uncalled (good assemblies)
Distance from the nearest HD1-class (KN) locus to the nearest STE3 locus on the same contig:

| Order | Group | n | HD and STE3 on one contig | Median distance (bp) | Within 150 kb |
|---|---|---|---|---|---|
| Tremellales | called | 317 | 134 | 21,031 | 113 |
| Trichosporonales | other calls | 10 | 7 | 54,370 | 7 |
| Trichosporonales | PR-only | 50 | 11 | 48,015 | 11 |
| Trichosporonales | uncalled | 63 | 34 | 49,013 | 34 |
| Filobasidiales | PR-only | 11 | 2 | 128,199 | 1 |
| Filobasidiales | uncalled | 21 | 2 | 325,108 | 0 |
| Cystofilobasidiales | PR-only | 14 | 5 | 504,613 | 1 |
| Cystofilobasidiales | uncalled | 20 | 1 | 196,117 | 0 |

- Trichosporonales: called and uncalled genomes have the same arrangement (HD1 and STE3 on one contig about 48 to 54 kb apart where both are found). The arrangement does not explain who fails. It is linked and compact, as expected for the fused locus of Sun et al. 2019, and it is also similar to Cryptococcus, where SXI1 and STE3 are about 60 kb apart in JEC21 (curated record).
- Filobasidiales and Cystofilobasidiales: HD and STE3 are mostly not on one contig or are hundreds of kb apart in both called and uncalled genomes (a tetrapolar-like or unlinked pattern). Again no difference between called and uncalled.
- The strict-CAAX comparison is not evidence. A PR-only call is admitted by the same strict-CAAX rule, so the receptors of PR-only genomes have a CAAX ORF nearby by construction (Trichosporonales 50/50, Filobasidiales 10/11, Cystofilobasidiales 14/14). Among uncalled genomes the share is 13/63 (21%) in Trichosporonales, 4/21 in Filobasidiales and 3/20 in Cystofilobasidiales, which is close to the chance level of the PR #33 CAAX study (about 13% of calls at chance across orders; excess share 0.83 in Trichosporonales). The uncalled receptors therefore carry no more CAAX signal than chance, and the rule says nothing about whether a locus is there.
- Receptor presence: STE3-like locus found in every good-assembly failure (50/50, 63/63, 11/11, 21/21, 14/14, 20/20, 10/10). HD1-class gene found (KN or strong miniprot) in Trichosporonales 15/50 PR-only and 37/63 uncalled, Filobasidiales 7/11 and 19/21, Cystofilobasidiales 13/14 and 11/20, Tremellales uncalled 7/10. Pfam KN misses some HD1 genes, so these are lower bounds.

### 2. What the pipeline did (v0.6.0 reports, good assemblies)
- Every uncalled good genome (10 Tremellales, 63 Trichosporonales, 21 Filobasidiales, 20 Cystofilobasidiales; 114 in all) has a withheld cluster with at least 2 genes found. The best withheld cluster is withheld by `modelled_gene_bar` in every case (10/10, 63/63, 21/21, 20/20). Polishing never models more than one gene in any withheld cluster (the maximum polished_genes over a genome's withheld clusters is 1 in 102 of the 104 uncalled good genomes outside Tremellales and 0 in 2; it is 1 in 7 and 0 in 3 of the 10 Tremellales).
- All 127 + 41 + 42 genomes in the three sister orders were routed through the phylum fallback, called genomes included. The fallback searches every Basidiomycota family, so a genome has about 300 withheld clusters. Medians by group (`followup_summary.txt`): Trichosporonales 332 (other calls, n=10), 314 (PR-only, n=50), 329 (uncalled good, n=63), 316 (all 127); Filobasidiales 289 (all 41); Cystofilobasidiales 258.5 (all 42). A median of 15 to 17 of these per genome are HD-pair clusters (HD1 plus HD2 or bE plus bW from homeodomain paralogs). Tremellales, routed by lineage, has a median of 2 withheld clusters over all 342 genomes, and 3.5 over the 10 uncalled good-assembly genomes only. Fallback routing explains the high noise level and why the modelled-gene bar is the gate that decides. It does not explain which genomes fail, because called and uncalled genomes of an order are routed identically and carry the same number of withheld clusters. In the fallback a real locus cannot be told from the noise, and relaxing the modelled-gene bar would release hundreds of clusters per genome.
- Tremellales uncalled (15): 14 withheld by the modelled-gene bar, 1 below the fraction floor (best cluster 0.33). The best cluster in the 10 good genomes already carries 2 to 4 genes of the family (examples: Naematelia encephala Treen1 FAO1, PAN6 and SXI2; Tremella fuciformis Tfnb_1.0 FAO1, MFa, SXI1 and SXI2; Cryptococcus neoformans ASM1105756v1 PAN6, MFa and STE3 with 0 polished). These are Tremella, Naematelia, Papiliotrema and divergent Cryptococcus lineages, far from the Cryptococcus deneoformans JEC21 references, so polishing models 0 or 1 gene.
- Called Trichosporonales: 6 calls use the MAT family (SXI1 plus STE3, MFalpha) and 3 use bLocus; 50 PR-only calls rest on a CAAX precursor.

### 3. Classification of the failures
Classes from the brief: (a) sensitivity, (b) pipeline gate, (c) architecture, (d) assembly, (e) absence. The gate class (b) applies to every good-assembly failure, because every one has a withheld cluster; the scan evidence splits the same genomes by what is behind the gate. A homeodomain (PF00046) near STE3 is now counted, so HD2-type alleles are no longer reported as "no HD". Counts at 150 kb (`classes_link150k.txt`, `failure_classes.tsv` and `withheld_pairs.tsv`, both regenerated at 150 kb):

| Order | Group | HD1 gene linked to STE3 (<=150 kb) | Homeodomain only, linked | HD1 found, not linked | No HD1 (KN) domain found | Poor assembly (d) | Absent (e) | Total |
|---|---|---|---|---|---|---|---|---|
| Tremellales | uncalled | 0 | 0 | 7 | 3 | 5 | 0 | 15 |
| Trichosporonales | uncalled | 34 | 7 | 2 | 20 | 4 | 0 | 67 |
| Trichosporonales | PR-only | 11 | 7 | 3 | 29 | 0 | 0 | 50 |
| Filobasidiales | uncalled | 0 | 2 | 17 | 2 | 1 | 0 | 22 |
| Filobasidiales | PR-only | 1 | 0 | 6 | 4 | 0 | 0 | 11 |
| Cystofilobasidiales | uncalled | 0 | 1 | 10 | 9 | 1 | 0 | 21 |
| Cystofilobasidiales | PR-only | 1 | 0 | 12 | 1 | 0 | 0 | 14 |
| Total | | 47 | 17 | 57 | 68 | 11 | 0 | 200 |

"No HD1 (KN) domain found" means that no PF05920 or strong miniprot HD gene was found anywhere and no PF00046 locus lies within 150 kb of STE3. It is not a statement that HD is absent. Per-genome distances are in `failure_classes.tsv` (`d_kn_ste3`, `d_hdstrong_ste3`, `d_hb_ste3`, `d_hdany_ste3`).

Link distance. The linked HD1 class is stable from 100 to 300 kb (47 at 100 and 150 kb, 49 at 200 and 300 kb). At 50 kb it falls to 28 because the 40 to 57 kb band, where most Trichosporonales HD1 to STE3 distances lie (median 48 to 54 kb), straddles the cut. 50 kb is not a meaningful boundary.

Reading:
- (e) true absence: 0 genomes. Every failure on a good assembly has a STE3-like locus.
- (d) assembly: 11 genomes (5 Tremellales, 4 Trichosporonales, 1 each in two orders), at N50 below 20 kb, BUSCO below 70 or both.
- (c) architecture: 64 genomes (47 with an HD1 gene and 17 with only a homeodomain) have HD and STE3 on one contig within 150 kb. The family definitions do not express this as one locus unless the MAT family is chosen (see Recommendations). The arrangement is found in both called and uncalled genomes, so in distance and contig there is no evidence that arrangement separates success from failure. Rearrangement was not tested beyond contig and distance: no orientation, no inversions, no FAO1 or PAN6 flank synteny. The conclusion is "no evidence for rearrangement in distance and contig", not "this is not a rearrangement problem". The 57 "HD1 found, not linked" genomes are mostly Filobasidiales and Cystofilobasidiales, where HD and PR are far apart or on different contigs in called and uncalled genomes alike.
- (a) and (b): for the remaining failures detection decides. Where no HD1 gene is found the curated HD references are too distant (strong miniprot hits at 45% or more to a curated HD in 0 of 63 uncalled Trichosporonales and 0 of 21 uncalled Filobasidiales; Cystofilobasidiales 1 of 20 uncalled, against 10 of 14 PR-only, which are close to the Phaffia HD references). For all good-assembly failures the gate (modelled_gene_bar after fallback routing) withholds a cluster that has the genes.

Allele confound (Trichosporonales, good assemblies; `followup_summary.txt`). Call status could track the mating-type allele: one allele carries an HD1-type gene, the other an HD2-type gene, and the receptor allele with a CAAX precursor might sit with one of them. Cross-tab of HD type by status (HD1-type = KN or strong miniprot; HD2-type = no HD1 gene but a homeodomain within 150 kb of STE3; neither otherwise):

| Group | n | HD1-type | HD2-type | Neither |
|---|---|---|---|---|
| Other calls | 10 | 7 | 1 | 2 |
| PR-only | 50 | 15 | 6 | 29 |
| Uncalled | 63 | 37 | 6 | 20 |

- HD1-type with no CAAX near STE3: 33 of 63 uncalled, 2 of 10 other calls, 0 of 50 PR-only. No HD1 gene with CAAX near STE3: 35 of 50 PR-only. The CAAX axis is circular for PR-only (see above); the KN axis is not: KN is found in 15 of 50 PR-only and 37 of 63 uncalled genomes.
- By genus (KN-positive of n): Trichosporon uncalled 17 of 19, PR-only 3 of 21; Pascua uncalled 6 of 6; Cutaneotrichosporon uncalled 13 of 23, PR-only 8 of 13; Apiotrichum uncalled 1 of 15, PR-only 1 of 10.
- Within Trichosporon the split is by species: T. austroamericanum (8) and T. ovoides (8) and T. faecale (1) are all PR-only and all KN-negative; T. asahii has 9 KN-positive and 2 KN-negative uncalled genomes and 2 KN-positive and 1 KN-negative PR-only genomes; the other Trichosporon species (8 genomes) are uncalled and KN-positive. Only 3 species in the order have both an uncalled genome and a called or PR-only genome (T. asahii, Cutaneotrichosporon cutaneum and C. dermatis), and in T. asahii both groups carry HD1-type genomes.
- Result: the call status tracks the clade and the HD type together (KN-positive genomes are the uncalled ones in Trichosporon and Pascua), which is what an allele effect would produce, but it is just as consistent with clade-level divergence of the receptor and precursor, and within T. asahii the HD type does not separate the groups. The existing tables cannot decide, because they have no allele per genome independent of detection, and Pfam KN false negatives and the CAAX rule are mixed in. Data that would decide: genomes of both mating types of one species (or a cross), an HD2 profile or lineage reference to type HD2 genes, an STE3 allele phylogeny (a1 versus a2 receptor clades) to test whether receptor allele predicts HD type, and precursor annotation by in-order tblastn instead of the CAAX rule to test whether the HD1-type receptors have precursors the rule misses.

Exemplars:
- Tremellales, gate with genes present: `GCA_011057565.1_ASM1105756v1` (Cryptococcus neoformans, 14 contigs, BUSCO 74.3, N50 1.19 Mb): best cluster carries PAN6, MFa and STE3 (0 polished genes), plus STE3 and SXI2 elsewhere. A C. neoformans genome failing is a gate failure, not a missing locus.
- Tremellales, divergent lineage: `GCA_002105065.1_Treen1` (Naematelia encephala, BUSCO 93.5, N50 209 kb): cluster of FAO1, PAN6 and SXI2 found, STE3 not found by the pipeline, but the scan finds 4 STE3 loci (2 with a CAAX ORF).
- Trichosporonales, linked and compact: `GCA_977066725.1_gfCryBron1.1` (Cryptotrichosporon brontae, 14 contigs): HD1-class gene 1.9 kb from STE3 on contig OZ370661.1. The pipeline called it as MAT (MFalpha, STE3, SXI1, RPL39) because the Tremellales references are close enough here. This is a real fused-style locus found by the existing MAT family.
- Trichosporonales, fused locus missed: `GCA_001712445.1_ASM171244v1` (Cutaneotrichosporon oleaginosum, uncalled; it is not the ATCC 20508 assembly): STE3 and HD1-class loci on KV757200.1 within 50 kb of each other, 28 CAAX ORFs genome wide; the pipeline reported 334 withheld clusters from the fallback, the best a receptor plus precursor cluster withheld by the modelled-gene bar.
- Cystofilobasidiales, unlinked: `GCA_001007165.2_Xden1` (Phaffia rhodozyma): HD1 (LN483167.1:58,713-60,570, 93% to the Phaffia KU3157 HD1) and HD2 (61,348-62,865, 100%) are adjacent, 0.8 kb apart, while the STE3 loci are on other contigs (LN483157.1, LN483332.1) and at 568 kb on the HD contig. The pipeline reported one PR call. The curated Phaffia HD references (HD1 and HD2 adjacent) are not in the database, so the HD pair is not found as a locus.

### 4. Direct test of in-order references (`inorder_summary.txt`)
Exon-level protein fragments (stop-to-stop six-frame segments with a Pfam KN or STE3 hit) from 4 called Trichosporonales genomes, 1 Filobasidium and 1 Cystofilobasidium genome, aligned with miniprot to all good genomes of that order (identity 0.5 or more, query coverage 0.5 or more):

| Order | Group | n | KN hit | STE3 hit | Linked within 150 kb |
|---|---|---|---|---|---|
| Trichosporonales | uncalled | 63 | 0 | 44 | 0 |
| Trichosporonales | PR-only | 50 | 2 | 17 | 2 |
| Filobasidiales | uncalled | 21 | 4 | 20 | 0 |
| Cystofilobasidiales | uncalled | 20 | 2 | 4 | 0 |

STE3 fragments from the same order align to 44 of 63 uncalled Trichosporonales genomes and 20 of 21 Filobasidiales. The KN fragments hardly align even in the called genomes of their own order (2 of 50), so this test does not measure HD sensitivity: a stop-to-stop segment is a single exon and does not align as a spliced protein. The test of HD sensitivity is therefore left open; full gene models are needed (see Recommendations).

## Limits
- The scan is the independent view, but Pfam KN and the 45% miniprot cut are themselves limited. "No HD1 (KN) domain found" means not detected, not absent. Genes were not modelled with introns, so locus counts can contain fragments. PF00046 matches many homeobox genes; it is used only within 150 kb of STE3.
- Rearrangement was tested only as contig and distance. No orientation or inversion analysis and no FAO1 or PAN6 flank synteny was done for any order. The finding is no evidence for rearrangement in distance and contig, not absence of rearrangement.
- Linked classes are stable from 100 to 300 kb; 50 kb splits the modal 40 to 57 kb band and is not a useful cut.
- 551 of 552 genomes were scanned; one scan failed. The run on current main with evidence diagnostics (jobs 29428627 to 29428629) was not finished at writing, so the new-code behaviour (polish tier, candidates per cluster) was not examined. The v0.6.0 reports predate the identity tier, which matters only for genomes with a strong cluster ranked outside the cap; no genome here has fewer than 2 withheld candidates, so the tier is unlikely to change these calls.
- The in-order reference test used exon-level fragments and is informative for STE3 only.
- PR-only calls are unverified by design (chance level about 13% across orders in the PR #33 CAAX test); any count of "recovered" genomes is an upper bound until verified.

## Recommendations
Ranked by genomes plausibly recovered (good-assembly failures; upper bounds, not run):
1. Order-level scope for the three sister orders (Trichosporonales, Filobasidiales, Cystofilobasidiales) so the genomes are routed by lineage instead of the phylum fallback: removes about 300 noise clusters per genome and makes the modelled-gene bar usable. Affects 104 uncalled and 75 PR-only good genomes (179). This is a precondition for the next two items. It explains the noise and the gate, not which genomes fail; alone it does not add genes.
2. In-order curated records from called genomes: a fused HD+PR record from Cryptotrichosporon brontae (OZ370661.1, HD1-class gene 1.9 kb from STE3; called by the existing MAT family) and C. argae (TUD_CryArg_v1.0, 5 kb), plus a Cutaneotrichosporon oleaginosum or C. curvatum locus, scoped to Trichosporonales. Recoveries (uncalled genomes, upper bound): 34 Trichosporonales genomes with an HD1 gene linked to STE3, plus up to 7 with a linked homeodomain only. Conversions to verified calls (already called as PR-only, so not recoveries): up to 11 linked PR-only genomes plus up to 7 with a linked homeodomain only. Add the Phaffia rhodozyma HD1 and HD2 deposits (KU315762 to KU315779, single genes of 1,860 and 1,239 bp) as an HD record for Cystofilobasidiales: they give 93 to 100% hits in Xden1, and 10 of 14 PR-only genomes have a strong HD hit (conversions), while only 1 of 20 uncalled genomes does. In the Filobasidiales, a Filobasidium wieringae or F. uniguttulatum STE3 record, with an HD1/HD2 pair taken from a called genome, since the arrangement is unlinked.
3. Tremellales: add one or two Tremella and Naematelia loci (the called Tremella and Naematelia genomes, or Papiliotrema) as additional MAT references, and for Cryptococcus ASM1105756v1 and similar check polishing. The 10 good uncalled genomes each have a cluster with 2 to 4 genes found; the cause is modelling failure (0 to 1 polished), so closer references are the fix. Up to 10 genomes. Also consider letting a cluster with at least 2 distinct found genes of a lineage-routed family pass at `low` when flank genes (FAO1 and PAN6 or STE3 plus SXI) are both present; this is a code rule that needs a regression gate.
4. A new record type for a fused HD+PR locus is not needed to describe the Trichosporonales arrangement: the MAT family already has HD (SXI1, SXI2) and receptor (STE3) genes with a 120 kb gap, and brontae is called by it. A code rule is needed only if the curator wants HD1 and HD2 genes inside the PR family called as one locus; that is a design decision.
5. Not fixable by references: the 11 poor assemblies (5 Tremellales, 4 Trichosporonales, 2 others), and HD genes that Pfam also fails to find (68 genomes with STE3 and no HD1 (KN) domain found) until lineage gene models exist. A reference-independent HD1/HD2 detector (Pfam KN plus HD2 homeobox with adjacent gene pair) could be tested as a rescue for these.

Open decisions: (i) approve order scopes for the three sister orders; (ii) choose which Trichosporonales and Phaffia records to curate; (iii) whether PR-only unverified calls stay in the count for Trichosporonales; (iv) whether to run the orientation and FAO1 or PAN6 flank-synteny analysis, which is the missing test of rearrangement (not run); (v) whether to obtain the allele data listed under the allele confound.
