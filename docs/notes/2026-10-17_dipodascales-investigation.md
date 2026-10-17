# Dipodascales: why v0.6.0 misses them, and what would fix it

Branch `dipodascales-refs`. Nothing here is merged and nothing was added to `db/` proper. Drafts are in `db/candidates/Ascomycota/`.
Raw outputs and scripts: `results/2026-10-17_dipodascales/` (`classified.tsv` has one row per Dipodascales genome; `sample24_hits.txt` has the per-genome hits).

## 0. Headline

1. The "rescue" signal in `results/2026-10-05_reference_gap_triage/` is real, but its queries were the wrong ones. Its Dipodascales hits came from *Ogataea* AJ617305 (35-40% identity) and *Y. lipolytica* AJ617307. Using the proteins of the order's own species (Y. lipolytica MATB, *Geotrichum candidum* MATA/MATB) the picture sharpens a lot.
2. The biggest, cheapest win is `Yarrowia`: 28 of 30 uncalled good-assembly Yarrowia genomes carry the **MATB** idiomorph (100% identical to AJ617307 in the CLIB122 reference genome and in strain 22301-5). Only MATA is curated, so every MATB strain is invisible.
3. Geotrichum (the other big genus, 91 genomes, 70 uncalled): the v0.6.0 "calls" are often in the wrong place. The real MAT genes sit next to SLA2 and match the two GenBank deposits `HF558448` (MATA) and `HF558449` (MATB) at 95-100% in 51 genomes (37 of them uncalled).
4. Across the 223 good-assembly uncalled genomes, 128 (57%) have a candidate MAT gene within 15 kb of SLA2 when queried with these proteins; 93 (42%) have SLA2 but nothing detectable beside it; 2 have no SLA2 hit. See section 2.
5. The Wickerhamiellaceae (*Wickerhamiella* 56 + *Starmerella* 34 uncalled, 51 and 9 of them good assemblies) are the real hole: 44 of 51 good-assembly Wickerhamiella and 8 of 9 Starmerella have no recognisable MAT gene by homology to anything we hold. Literature says the MAT genes there are next to a different flank (TFC1), or are pheromone-gene linked.

## 1. What the literature and deposits say

Verified = I opened the record or abstract in this session (NCBI E-utilities or PubMed). "Secondary" = taken from the repo's earlier notes or a search snippet. "Unverified" = could not confirm.

| Taxon | What is known | Evidence |
|---|---|---|
| *Yarrowia lipolytica* | Heterothallic, two idiomorphs MATA and MATB at the same site (SLA2 ... APN2), each carrying two ORFs. MATA: sla2, mata1 (119 aa), mata2 (291 aa, HMG box), apn2. MATB: matb2 and matb1. | MATA locus `AJ617306.1` (CAE84421-24) and `AJ007708` (CAA07612/13), PMID 1588911, 10503550. MATB locus `AJ617307.1` (3,353 bp): matb2 CDS 259..1362 `CAE84425.1`, matb1 CDS 2362..3270 `CAE84426.1`; PMID 14745027 (Butler 2004). All verified by efetch. MATB functional in a MATA background by cassette swap: PMID 18461384 (Rosas-Quijano 2008). Scarless A<->B switching and mating: PMID 32105315 (Han 2020). |
| *Y. lipolytica* assemblies | E122 = MATA, 22301-5 = MATB, CLIB122/E150 = MATB, W29 = MATA. Haplotype-resolved MATA and MATB reference assemblies with Dryad data, not INSDC. | Zali et al., bioRxiv 10.1101/2025.03.24.645116 (read from the PDF: genotypes and the Sla2/Apn2 flank statement). Journal version (NAR Genom Bioinform 7:lqaf175) is secondary (web search). Sequence confirms: CLIB122 `NC_006069.1:979564-987371` and 22301-5 `CP096873.1:1023801-1031608` are 100% to AJ617307. |
| *Geotrichum candidum* (= *Galactomyces candidus*) | Two single-gene deposits: `HF558448.1` MATA, 846 bp, CLIB 918, protein `CCO73812.1`; `HF558449.1` MATB, 1,020 bp, CBS 615.84, protein `CCO73813.1`. Both unpublished (INRA submission, Nov 2012). In my screen MATA and MATB occur in different genomes, never together (section 2), which fits heterothallism. | Records verified by efetch. Heterothallic status itself is **unverified**: a snippet from a review (Perkins et al. 2020, Front Microbiol, PMC7220993) says the species is "self-fertile", and I could not open the full text. McKay et al. 2012 (PMID 22494250) found 1:1 mating-type segregation in the sister species *G. citri-aurantii* by assay. |
| *Magnusiomyces capitatus* | Described as heterothallic, with SNV patterns taken as evidence of recombination. No MAT locus sequence. | PMID 30456648 (Brejova 2018), abstract read. |
| *Dipodascus (Magnusiomyces) starmeri* | Haploid, heterothallic. Type NRRL Y-17816 (= CBS 780.96; BFD genome `GCA_030578655.1`), complementary mating type UCD-FST 81-513.3 (NRRL Y-17817, not in BFD). | PMID 9103614, abstract read. |
| *D. tetrasporeus* | Homothallic. | PMID 18398216, abstract read. |
| *Blastobotrys* | MATalpha1 between SLA2 and APN2 (B. adeninivorans, B. proliferans); MATalpha2 absent in the genus; B. proliferans PHN; B. pratensis lost MATa genes. | Secondary: repo note `2026-09-24_mat-reference-gap-literature-round2.md`, citing Krassowski 2019 (PMID 31353182, abstract read; I did not open the supplement). I confirmed the SLA2-alpha1-APN2 arrangement in *B. adeninivorans* PETY001 (section 2). |
| *Starmerella / Wickerhamiella / Candida apicola* | S. bombicola: SLA2-ORF-TFC1 with no identifiable a1/a2. C. apicola, W. domercqiae: SLA2-MATalpha1-MFalpha-TFC1 (pheromone gene at MAT). | Secondary (same repo note). |
| *Nadsonia fulvescens* var. *elongata* | PHN (pseudo-homothallic, no switching). | Secondary (same note). |
| *Trichomonascus vanleenenianus* CLIB 3782 | RefSeq annotation lists a "putative alpha-factor mating pheromone" (`XM_076682001`) and SLA2 (`XM_076681997`) with adjacent locus tags (TRVA0_001S07646 and a nearby number). A neighbour check on the contig is a cheap follow-up. | Search result only, **unverified** adjacency. Genome paper PMID 41864865 does not discuss MAT in its abstract. |
| Galactomyces, Saprochaete, Dipodascus other than the above, Nadsonia (beyond Krassowski), Sugiyamaella, Magnusiomyces | No MAT locus or protein deposit found (protein and nucleotide searches by organism and title returned only pheromone receptors, STE2 and STE3, for Sugiyamaella, Blastobotrys, Wickerhamiella, Starmerella, Nadsonia). | NCBI esearch, this session. |

Nomenclature trap: "MATA/MATB" means the two *idiomorphs* in Yarrowia (neither is a homeodomain) and the two *proteins* in G. candidum. The G. candidum MATA protein hits Y. lipolytica mata2 (HMG box) at 29%, and the MATB protein hits Ogataea matalpha1 (alpha box) at 35-40%. So by position and sequence, G. candidum looks like an HMG/alpha-box pair at the SLA2 site. Domain content was not checked with Pfam.

BFD assemblies matching characterised strains (from `samples.parquet` strain field):

| Strain | Assembly | Status v0.6.0 |
|---|---|---|
| Y. lipolytica W29 / CLIB89 (MATA) | `GCA_001307795.1`, `GCF_001761485.1`, `GCA_003054345.1` | called MATyl A |
| Y. lipolytica E122 (MATA) | `GCA_059623445.1_CSHL_YarrowiaLipolytica_E122_MATA` | called |
| Y. lipolytica CLIB122 (MATB) | `GCF_000002525.2_ASM252v1` | **uncalled** (MATB 100%) |
| Y. lipolytica 22301-5 (MATB) | `GCA_023272835.1_ASM2327283v1` | **uncalled** (MATB 100%) |
| G. candidum CLIB 918 (MATA deposit) | `GCA_001402995.1_New2.3_08062011` | called, but at a different contig than the MATA gene |
| G. candidum CBS 615.84 (MATB deposit) | not in BFD (no match for "615.84") | - |
| M. starmeri NRRL Y-17816 | `GCA_030578655.1_ASM3057865v1` | uncalled |

## 2. Evidence in the genomes

Method: tblastn (e <= 1e-5) of 35 proteins (`dipo_queries.faa`: the Y. lipolytica MATA and MATB proteins, G. candidum MATA and MATB, plus the earlier Ogataea/Kphaffii/Cfim rescue queries) against all 398 Dipodascales genomes. Job arrays `dipo-blast` (24 tasks) and `dipo-all` (40 tasks), about 12 s per genome. A "candidate" is a core-query hit with >= 30% query coverage on the same contig within 15 kb of the SLA2 hit, excluding the SLA2 gene itself. This is a homology screen; no gene models were built.

Good-assembly uncalled genomes (BUSCO >= 70 and N50 >= 20 kb; my approximation of the triage filter, 223 vs the earlier 224):

| Candidate type | n | Where |
|---|---|---|
| Yarrowia-clade MATB (matb2 >= 60%) | 28 | all Yarrowia but 2 |
| G. candidum-like MATA (HMG), 40-100% | 14 | Geotrichum 13, Magnusiomyces 1 |
| alpha1-like beside SLA2 (any identity >= 21%) | 86 | Geotrichum 27, Blastobotrys 13, Sugiyamaella 10, Magnusiomyces 9, Zygoascus 5, Diddensiella 3, Wickerhamiella 6, Entelexis 3, Spencermartinsiella 2, others |
| SLA2 found, no core hit within 15 kb | 93 | Wickerhamiella 44, Blastobotrys 9, Starmerella 8, Geotrichum 7, Entelexis 6, Nadsonia 4, others |
| no SLA2 hit | 2 | |

Excluding the G. candidum-like hits (>= 85%), the median identity of the 64 remaining alpha1-like hits is 33%, with coverage often 0.35-0.65 of the query. They cluster tightly against SLA2 (often within 0.3-1 kb), which is the main reason to believe they are real.

Called genomes: of 113 called Dipodascales genomes that have an SLA2 hit, 34 have the call more than 30 kb from SLA2 (Geotrichum 18). 16 of the offending loci are `idiomorph_gene_only` with `MAT1-2-1|MAT1-1-3` on a contig without SLA2/APN2. These are probably HMG-box paralogs, not the MAT locus. That is why the triage's Dipodascales called controls (5 with a core hit) had 0 concordant calls. It also means "called" is not equal to "correct" in this order, which matters for evaluation.

Examples (genome id is the BFD file stem; "S" = the SLA2 hit, "A" = APN2 hit; identities from tblastn of the named query):

| Genome | Taxon, strain | Best hit (query, id, cov, contig:start-end) | Flanks | My read |
|---|---|---|---|---|
| `GCF_000002525.2_ASM252v1` | Y. lipolytica CLIB122 | matb2 100/1.0 NC_006069.1:981060-982160; matb1 100/1.0 983163-984068 | A 979564-980799, S 984213-987371 | **present**, MATB, reference genome, uncalled |
| `GCA_023272835.1_ASM2327283v1` | Y. lipolytica 22301-5 | matb2 100/1.0 CP096873.1:1025297-1026397 | A 1023801, S 1028450-1031608 | **present**, MATB |
| `GCA_900518995.1_YAYA0` | Yarrowia yakushimensis CBS 10253 | matb2 82/0.99 ULGW01000001.1:4455568-4456659; matb1 65/0.98 | A 67%, S 96% | **present**, MATB-like in a sister species |
| `GCA_013365045.1_LMA-244_clib` | Geotrichum candidum LMA-244 | HF558449 MATB 97/1.0 QQZM01000017.1:285536-286546 | S 282113-285301; APN2 not adjacent | **present**, MATB |
| `GCA_001402995.1_New2.3_08062011` | G. candidum CLIB 918 | HF558448 MATA 100/1.0 CCBN010000002.1:584221-585063 | S 579067-582264 | **present**, MATA; v0.6.0 called a different contig |
| `GCA_022023575.1_ASM2202357v1` | Geotrichum ghanense CBS 11010 | MATB-like 52/0.53 JAKKCX010000006.1:1653133-1653660; MATA-like 55/0.47 1669226-1669624 | S 1672203-1675394 | **unclear**: both idiomorph-like genes within 20 kb of SLA2 (possible homothallic arrangement; not modelled) |
| `GCA_030563945.1_ASM3056394v1` | G. australiense NRRL Y-17565 | MATB-like 57/0.50 JAJMKR010000227.1:3833-4327; MATA-like 47/0.92 JAJMKR010000162.1:41275-42051 | S on contig 53 | **unclear**: fragmented (N50 169 kb), hits on contigs without SLA2 |
| `GCF_902498895.1_sapIngB` | Magnusiomyces paraingens CBS 517.90 | MATB-like 37/0.80 NW_022474068.1:4138300-4139052 | S 4134847-4138050 | **present** (alpha1-like adjacent to SLA2) |
| `GCA_900497725.1_magCapA` | M. capitatus NRRL Y-17686 | MATB-like 46/0.49 UIDF01000003.1:1939738-1940214 | S 1931067-1934279 | **present**, probable |
| `GCA_030558955.1_ASM3055895v1` | Diddensiella caesifluorescens NRRL Y-48781 | MATB-like 35/0.57 JAKTQR010000001.1:57640-58179 | A 54903-56075, S 58520-61588 | **present**: APN2, alpha1-like, SLA2 triad |
| `GCA_055672235.1_ASM5567223v1` | Blastobotrys adeninivorans PETY001 | matalpha1 29/0.55 JBNHNI010000001.1:2744895-2745203 | S 2741327-2744437, A 2746283-2747347 | **present**, matches Krassowski's SLA2-alpha1-APN2 |
| `GCA_030579815.1_ASM3057981v1` | Sugiyamaella smithiae NRRL Y-17850 | MATB-like 30/0.65 JAKTRV010000001.1:65529-66122 | A 63652-64680, S 66422-69556 | **present** |
| `GCA_030563625.1_ASM3056362v1` | Zygoascus polysorbophilus NRRL Y-27161T | MATA-like 31/0.38 JAJMJP010000004.1:386769-387083 | S 382719-385793, A 387808-388884 | **probable**: HMG-like ORF between SLA2 and APN2, low identity |
| `GCA_001661315.1_Nadfu1` | Nadsonia fulvescens DSM 6958 | MATB-like 28/0.60 KV454408.1:1412060-1412587 (other contig) | S 1137466-1140660, A 1134110-1135360 on KV454406.1 | **unclear/absent** at the SLA2 site; hit is elsewhere. Krassowski reports PHN for the variety |
| `GCA_946409495.1_Wolfe_Bari_UCD613` | Blastobotrys aristatus UCD613 | none | S 3506297-3509341, A 3513054 (cov 0.40) | **unclear**: 3.7 kb SLA2-APN2 gap with no core hit; could be an unrecognised gene |
| `GCA_030583385.1_ASM3058338v1` | Spencermartinsiella ligniputridi NRRL Y-48818 | none | S 10111-13251, A 15682-16797 (2.4 kb apart) | **unclear**: short interval; ORF scan needed |
| `GCA_056749365.1_ASM5674936v1` | Starmerella batistae SB001 | none (best: matalpha2 32/0.46 on another chromosome) | S CM164359.1:2746452-2749571; APN2 580 kb away | **unclear/absent**: BioSample says MATa (`results/2026-10-17_biosample_truth/truth_set.tsv`), so a locus exists but is undetectable with our queries |
| `GCA_056822755.1_Ther_Cv-P4_1` | Wickerhamiella versatilis C0-1-P4 | none | S and A on different contigs | **unclear**: wrong flank (needs TFC1) |
| `GCA_030762955.1_ASM3076295v1` | Entelexis paramagnoliae | none | S, A 2.5 Mb apart | **unclear/absent** |

Flank genes: DIC1 (an Ogataea flank in the earlier queries) is NOT near SLA2 in any Dipodascales genome I looked at; it is on other contigs. COX13 was not queried. APN2 is conserved enough to find in Yarrowia (67-100%) but diverged elsewhere (28-40%, often under 50% coverage), so SLA2 is the usable anchor.

## 3. Proposal

**Family and scope.**
- Option 1 (minimal): extend `MATyl` from `[4952]` to Yarrowia genus `[4951]` and add the MATB record. Fixes 28 genomes plus 2 alpha1-like Yarrowia.
- Option 2 (new family): `MATdc` for the rest of Dipodascales with `taxonomic_scope` on genera or families, idiomorphs `["A","B"]` (A = HMG protein, B = alpha1-like), `flank_carried_window_bp: 20000`, genes `MATA` (core_MAT, A), `MATB` (core_MAT, B), `sla2` and `apn2` (flanking_conserved; apn2 optional, see above). Taxids (verified by esummary): Dipodascales 3243772, Dipodascomycetes 3239873, Dipodascaceae 34353, Trichomonascaceae 410830, Yarrowia 4951, Geotrichum 43987. If scope is the whole order 3243772, `route()` returns MATyl and MATdc for Yarrowia and reports twice (the same problem the Serinales comment in `order.yml` describes), so either list genera, or fold everything into one family. Wickerhamiellaceae has no taxid under that name in NCBI taxonomy (esearch returned nothing); Starmerella 75735 and Wickerhamiella 45787 are the genus ids.
- Curated records (drafts in `db/candidates/Ascomycota/`): `4952_22301-5_MATyl_B`, `43987_clib918_MATdc_A`, `43987_lma244_MATdc_B`. The Geotrichum records are tier-3 (annotated deposit, no paper). Further records that would widen coverage, all to be built by miniprot from the genomes since no deposit exists: *Blastobotrys adeninivorans* alpha1, *Magnusiomyces paraingens* alpha1, *Diddensiella*/*Sugiyamaella* alpha1, *Zygoascus* HMG.

**Profile-HMM stage.** Yes, probably the best next step for the non-Yarrowia, non-Geotrichum 100+ genomes. The alpha1-like proteins sit at 30-40% identity from the query, near the tblastn noise floor, so a single-query search is fragile. Cheap design: (a) build a profile from the ~86 SLA2-adjacent alpha1-like ORFs (plus Ogataea matalpha1) and a second from the HMG ORFs; (b) search only a +-15 kb window around SLA2, with six-frame ORFs so unannotated short genes count; (c) require the hit to be adjacent to SLA2 (this is what keeps precision up; the repo already measured the risk of unanchored HMG hits, since Dipodascales got `MAT1-2-1` paralog calls). Pfam PF04769 (MATalpha1) is mentioned in `order.yml` as hitting curated MTLalpha1s; testing it on the 86 ORFs is a one-hour check I have not done.

**Wickerhamiellaceae need a different anchor.** SLA2 alone gave nothing for 44 of 51 good-assembly uncalled Wickerhamiella. Per Krassowski the locus is SLA2-MFalpha-TFC1 or SLA2-ORF-TFC1; adding TFC1 as a flank and a pheromone-ORF search (the repo already has `2026-09-20_short-pheromone-orf-detection.md`) is the likely route. No deposits exist, so records would be genome-derived.

**Evaluation plan.** There is almost no independent truth for this order, so use layers:
1. BioSample truth set (`results/2026-10-17_biosample_truth/`): 2 Dipodascales genomes (Y. lipolytica W29 MATa, which agrees on A; Starmerella batistae SB001 MATa, currently undetectable). Too small alone.
2. Yarrowia: 80 genomes, all with a MATA or MATB gene pair at SLA2 in my screen. Strain genealogy from Zali 2025 gives independent labels for E122, 22301-5, CLIB122/E150, W29 (and Po1f MATA per the same preprint, via search snippet). Hold out the MATB record, then check that MATB genomes are called MATB. Note 4 "called" Yarrowia are alpha1-like only and 3 are called elsewhere; inspect those.
3. G. candidum: about 70 genomes, bimodal MATA vs MATB at SLA2. Only CLIB 918 has an independent label (the deposit). A useful negative check is that no genome should carry both (G. ghanense is the exception to look at). Many are cheese-rind metagenome bins (`AOP*_surf_bin_*`); treat those as lower-confidence.
4. Species-level expectation from Krassowski 2019 supplement (heterothallic, PHN, PHC per species), which I have not extracted. It would give a "does the tool find a locus where the paper says one exists" metric for ~20 BFD species.
5. Leave-one-genus-out: build the profile from all but one genus and score recovery of the SLA2-adjacent ORF in the held-out genus.
Report two numbers per layer: recall (locus found within 30 kb of SLA2) and concordance (call within 30 kb of the SLA2-adjacent gene), because "called" alone was not a reliable correctness proxy.

## 4. Questions for the curator

1. **Yarrowia scope.** Extend MATyl to the genus (taxid 4951) or keep species-only? Evidence: *Y. yakushimensis* CBS 10253 carries MATB at 82/65% and SLA2 96%, and Yarrowia sp. E02 is in the set; both are outside `[4952]`; 49 Yarrowia genomes are called, 46 of them as MATyl.
2. **One family or two.** MATyl (HMG/unclassified, Yarrowia) and MATdc (HMG + alpha1-like, rest of the order) share the SLA2 site and, in my read, probably homology, but Yarrowia's mata1/matb1/matb2 have no known domains (your ruling 2026-10-01 leaves them unassigned). Merge into one `MATdip` family with several records, as you ruled for Serinales, or keep two with disjoint genus scopes?
3. **Trust level of the Geotrichum deposits.** `HF558448/9` are unpublished and annotated only as "mating type protein". I found 95-100% matches beside SLA2 in 51 genomes with a clean MATA/MATB split, which is strong circumstantial support. Accept as tier 3, or hold until Pfam confirms HMG/alpha box?
4. **idiomorph_class.** If Pfam confirms HMG and alpha1 boxes, set A -> MAT1-2 and B -> MAT1-1 for MATdc? For Yarrowia I did not check and do not recommend a mapping.
5. **Fix v0.6.0's wrong-contig calls?** 34 called genomes (18 Geotrichum) have the call far from SLA2, with 16 `MAT1-2-1|MAT1-1-3` gene-only calls. Do you want these demoted or suppressed once a MATdc family exists (the Pezizomycotina `MAT` family routing via phylum_fallback is the cause)?
6. **Wickerhamiellaceae.** Spend effort on a TFC1-anchored/pheromone route (90 genomes, no deposits, Krassowski only), or declare this clade out of scope for a first release and report "no_reference"?
7. **Homothallic-like patterns.** G. ghanense shows both alpha1-like and HMG-like candidates within 20 kb of SLA2. Should such genomes be reported as `homothallic_candidate` the way the homothallic rule treats Serinales?
8. **Is the 2025 Yarrowia MATA/MATB assembly worth pulling from Dryad** (doi 10.5061/dryad.v15dv427f per search) for a fully annotated MATB record with UTRs, rather than using the 1.1 kb deposit plus the 22301-5 span?
9. **Poor-assembly cutoff.** 59 of 283 uncalled genomes fail BUSCO 70 or N50 20 kb. Several Wickerhamiellaceae pass BUSCO only at 70-73%, since Dipodascomycete BUSCO scores run low in this order. Re-tune the threshold per class or leave it?

## Caveats

- All hits are tblastn homology, no gene models, no Pfam. "Present" means a plausible MAT gene in the expected place, not a validated locus.
- Strain-to-assembly matching used the `STRAIN` field in `samples.parquet`; CLIB 918 is `GCA_001402995.1`, the CBS 615.84 strain was not found.
- Krassowski 2019 details are second-hand via the repo's round-2 note; I did not read the paper's supplement.
- Zali 2025 journal details and the Dryad DOI come from a web search snippet, not a fetched page.
- The working tree at `/Users/jstajich/projects/MATPredict` is on `uncalled-triage` with an untracked `scripts/triage_uncalled.py` and `results/2026-10-17_parastagonospora_check/`; I did not touch them. This branch was made in a separate worktree from `main`.
