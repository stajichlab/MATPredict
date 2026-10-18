# Dothideomycetes MAT1-1 / MAT1-2 calls: audit of current main (2026-10-18)

Analysis only. `src/` and `db/` are untouched, no code is proposed. Branch `dothideomycetes-mat-audit`, data and scripts in
`results/2026-10-18_dothideomycetes_audit/` (every table below is a file under `tables/`; `run_all.sh` rebuilds them from
`loci_audit.tsv`, `genome_status.tsv` and `scan/`). The question: which idiomorph calls MATPredict makes in Dothideomycetes are
right, which are wrong, and which cannot be judged. "Measured" means a number from the runs; "judgement" is marked as such.

## Short answer

Of the 2,730 Dothideomycetes genomes in BFD that have a genome file, current main (`--partial-core-policy coverage_profile`)
gives a single MAT1-1 or MAT1-2 call to 2,474 (90.6%). Against evidence that does not come from the vote itself (alpha-box and
HMG-box profiles built without the genome's own genus where one exists, Pfam domains on the called protein, and an
annotation-free six-frame ORF scan of every genome):

| verdict (per genome, primary locus) | genomes | what it means |
|---|---|---|
| confirmed | 1,969 | vote, profile and domain agree; no sign of a second idiomorph's gene at the locus |
| conflicting | 298 | the profile disagrees, or the other idiomorph's gene is at the same locus (homothallic-like), or two loci carry different idiomorphs, or a core gene was set aside for an optional gene |
| unverifiable | 207 | weak second-gene signal, domain absent on a gene fragment, or profile undetermined |
| no call | 256 | 251 not detected, 3 undetermined only, 2 no report |

The rule for these four classes is in "Classification" below. The calls are right where the genus is heterothallic and the locus
holds one gene (Parastagonospora, Alternaria, Venturia, Zymoseptoria, Ascochyta, Fulvia, Corynespora, Macrophomina, Hortaea,
Dothistroma, Pyrenophora teres, and most others). The **systematic failure is homothallic or fused arrangements reported as one
idiomorph**, plus one optional-gene supersession, one near-total no-call group, and two thinly supported groups:

1. **Homothallic-like loci called as one idiomorph, 298 genomes in 43 genera.** Aureobasidium (119 of 187 conflicting), Pyrenophora
   tritici-repentis (62 of 64), Didymella (14 of 30, all 8 *D. pinodes*), Friedmanniomyces (24 of 26), Elasticomyces (9 of 15),
   Phyllosticta (12 of 56, all 8 *P. capitalensis*), Stemphylium (5 of 5 called). In each, both a complete alpha-box gene and an
   HMG-box ORF sit within 20 kb of the called locus; the label that appears is whichever scores higher (Aureobasidium pullulans
   37 MAT1-1 : 46 MAT1-2, A. melanogenum 32 : 27). No locus in this class carries the class `homothallic_candidate` (0 of 2,592),
   `two_idiomorphs` is empty in all 2,728 reports, and 183 of the 298 conflicting genomes are `high` confidence.
2. **Optional gene supersedes the real MAT1-2-1: 65 loci**, 59 of them *Pyrenophora tritici-repentis* (all 64 genomes are called
   MAT1-1, none MAT1-2). The `MAT1-2-1`/`MAT1-1-3` pair is resolved against the HMG gene. *P. tritici-repentis* is homothallic
   (PubMed 20373969, 29050833), and the scan finds an HMG ORF at the locus in 62 of 64.
3. **Cochliobolus clade almost never called.** Bipolaris 1 of 58, Curvularia 1 of 24, Exserohilum 0 of 10, Ophiosphaerella 0 of
   11, including all four genomes of the strains the curated records come from (C4, C5). Recall, not idiomorph; cannot be judged.
4. **Botryosphaeria dothidea: 140 genomes called MAT1-1, none MAT1-2, cannot be judged.** The alpha-box gene is real (Pfam 134
   bits). 137 of the 140 genomes also carry one HMG-box ORF (238 aa, 69 bits on the HMG profile) at the locus; 126 of the 140
   genomes come from one BioProject. Heterothallic MAT1-1 loci in the baseline groups never have one (0 of 715). It is either a
   MAT1-1-3-like HMG gene or the MAT1-2-1 of a homothallic arrangement; the data here do not separate them.
5. **Macrophomina MAT1-2 (239 genomes) rests on an 11% coverage, 39% identity hit.** For 231 of them the pipeline's gene model has
   no Pfam HMG hit; the scan finds an HMG ORF at the locus (about 100 bits) in all 239 and no alpha-box ORF anywhere in 238 of
   them. I read the call as probably right, but its support is the scan, not the pipeline's own evidence.

Independent truth is thin: BioSample states a mating type for 2 of 2,731 genomes (both called correctly), and 5 further
genomes are matched to a strain in a curated record (IPO323 called correctly; the 4 Cochliobolus C4/C5 genomes not called).
Everything else is internal consistency between methods plus species-level literature.

## What was run

| step | what | where |
|---|---|---|
| taxid | `taxonkit lineage` under Slurm (job 30049146): 147541 is *Dothideomycetes*, rank class | `taxid_147541.txt`, `dothideomycetes_orders_taxonkit.txt` |
| genome list | BFD `ncbi_accessions_taxonomy.csv` CLASS = Dothideomycetes (2,731) union genomes whose taxid lies under 147541 (672 taxids, `taxids_under_147541.txt`); FASTA from `Fungi_BFD_runs/input_clean_genomes` | `genome_list.tsv` |
| run | origin/main `c9c0ec6` (matpredict 0.6.2, database sha `cd4197dc...`, 171 records), `detect --taxid T --partial-core-policy coverage_profile`, lineage routing, every other default; 110 array tasks (job 30049176), 2 CPUs, 8 GB, partition short, at most 14 at a time; per-genome limit 3,600 s; median 127 s, p90 170 s | `reports_all.tar.zst` (HPCC logs and per-genome diagnostics not committed) |
| extraction | every detected locus: report fields, per-gene identity, coverage, bitscore, status, contig-edge distance; Slurm job 30056459 | `loci_audit.tsv` (2,595 loci), `genome_status.tsv` |
| profile | alpha-box (MAT1-1-1, 22 sequences) and HMG-box (MAT1-2-1, 26) profiles rebuilt without the genome's own genus when the genus has a curated record (Bipolaris, Parastagonospora, Plenodomus, Pseudocercospora, Zymoseptoria: 287 loci), else the shipped profiles (2,305 loci); order-out profile (every curated record of the order withheld) as a second view for Pleosporales and Mycosphaerellales; margin 25 bits as in `score_core_profiles.py` | columns `prof_*` |
| Pfam | Pfam 38.2 (2026-01-27) models MATalpha_HMGbox PF04769, MAT1-1-2 PF17043, MATalpha3_N PF28944 and the HMG-box clan CL0114 (13 models), gathering cutoffs, on the called MAT1-1-1 and MAT1-2-1 proteins | `make_pfam_subset.py`, columns `pfam_*` |
| genome-wide scan | six-frame translation of all 2,730 genomes, stop-to-stop ORFs of at least 50 aa, scored with the two profiles and Pfam; reports every alpha-box or HMG-box ORF of at least 40 bits with its position (28 array tasks, job 30053543) | `scan/`, `scan_genome.py` |
| truth | BioSample text of all 2,717 genomes that have one (`biosample/`); strain-matched literature statements | `truth.py` |
| comparison runs | v0.6.0 campaign (`v060/`) and the 2026-10-05 with-records run (`v060/genomes_records_run_20261005.tsv`) | |

No NCBI e-mail address from the user's configuration was used for the BioSample fetch.

## 1. Genome set and call rates

The CLASS column holds 2,731 genomes; taxonkit adds none (no genome with a blank CLASS lies under 147541) and removes none, apart
from one genome labelled Dothideomycetes whose taxid is not under 147541 (GCA_059628725.1, *Cladosporium* sp. MA8_LABAP, taxid
3751182; it was run). One genome has no FASTA (GCA_982436715.1, *Aureobasidium pullulans* apul_EXF-150). Two produced no report:
GCA_986281625.1 (*Cladosporium halotolerans*, empty FASTA) and GCA_016165855.1 (*Pseudocercospora ulei*, 3,600 s limit); both were
also missing from the 2026-10-05 run. The set is 24 orders (BFD ORDER, blank ones filled by taxonkit), 209 genera, 602 species.

Per order (full table `tables/note_orders.md`; `v0.6.0 called` counts genomes with a MAT1-1 or MAT1-2 locus in the v0.6.0 campaign):

| order | genomes | report | called | MAT1-1 | MAT1-2 | two-idiomorph loci | not detected | v0.6.0 called | confirmed | conflicting | unverifiable |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Pleosporales | 1010 | 1010 | 873 | 500 | 373 | 3 | 136 | 854 | 759 | 99 | 15 |
| Botryosphaeriales | 807 | 807 | 758 | 464 | 294 | 0 | 49 | 753 | 604 | 13 | 141 |
| Mycosphaerellales | 416 | 415 | 376 | 205 | 171 | 35 | 39 | 347 | 308 | 53 | 15 |
| Dothideales | 202 | 202 | 195 | 114 | 81 | 3 | 5 | 197 | 46 | 120 | 29 |
| Venturiales | 112 | 112 | 109 | 45 | 64 | 0 | 3 | 109 | 108 | 0 | 1 |
| Cladosporiales | 88 | 87 | 85 | 41 | 44 | 5 | 2 | 85 | 77 | 8 | 0 |
| Myriangiales | 28 | 28 | 28 | 7 | 21 | 0 | 0 | 28 | 24 | 0 | 4 |
| unclassified Dothideomycetes order | 26 | 26 | 22 | 13 | 9 | 1 | 4 | 20 | 18 | 3 | 1 |
| Trypetheliales | 13 | 13 | 8 | 6 | 2 | 0 | 5 | 8 | 7 | 1 | 0 |
| Capnodiales | 6 | 6 | 4 | 0 | 4 | 0 | 2 | 3 | 4 | 0 | 0 |
| 14 further orders with 1 to 3 genomes each (Gloniales, Mytilinidiales, Hysteriales, Aulographales, Neophaeothecales, Eremomycetales, Lichenotheliales, Microthyriales, Phaeotrichales, Acrospermales, Valsariales, Patellariales, Lineolatales, Tubeufiales) | 22 | 22 | 16 | 7 | 9 | 0 | 6 | 16 | 14 | 1 | 1 |
| **all** | 2730 | 2728 | 2474 | 1402 | 1072 | 47 | 251 | 2420 | 1969 | 298 | 207 |

The 14 small orders are listed one by one in `tables/note_orders.md` (Aulographales 2 genomes, 1 called; Hysteriales 2, 1; Tubeufiales 1,
none; ...); none has enough genomes to say anything about a group. Genera with at least 10 genomes (`tables/note_genera_ge10.md`; all 209 in
`tables/note_genera_all.tsv`):

| order | genus | genomes | called | MAT1-1 | MAT1-2 | not detected | confirmed | conflicting | unverifiable | group label |
|---|---|---|---|---|---|---|---|---|---|---|
| Botryosphaeriales | Macrophomina | 491 | 490 | 251 | 239 | 1 | 490 | 0 | 0 | confirmed |
| Pleosporales | Alternaria | 279 | 267 | 159 | 108 | 12 | 265 | 0 | 2 | confirmed |
| Pleosporales | Parastagonospora | 188 | 184 | 92 | 92 | 4 | 184 | 0 | 0 | confirmed |
| Dothideales | Aureobasidium | 187 | 180 | 105 | 75 | 5 | 32 | 119 | 29 | **conflicting** |
| Botryosphaeriales | Botryosphaeria | 144 | 143 | 143 | 0 | 1 | 4 | 0 | 139 | unverifiable |
| Venturiales | Venturia | 106 | 104 | 41 | 63 | 2 | 104 | 0 | 0 | confirmed |
| Pleosporales | Pyrenophora | 90 | 90 | 76 | 14 | 0 | 28 | 62 | 0 | **conflicting** |
| Cladosporiales | Cladosporium | 80 | 79 | 38 | 41 | 0 | 74 | 5 | 0 | mixed |
| Mycosphaerellales | Hortaea | 78 | 68 | 37 | 31 | 10 | 66 | 1 | 1 | mixed |
| Mycosphaerellales | Zymoseptoria | 76 | 76 | 41 | 35 | 0 | 71 | 0 | 5 | confirmed |
| Pleosporales | Corynespora | 58 | 58 | 27 | 31 | 0 | 50 | 0 | 8 | mixed |
| Pleosporales | Bipolaris | 58 | 1 | 1 | 0 | 57 | 1 | 0 | 0 | confirmed, call rate 2% |
| Botryosphaeriales | Phyllosticta | 56 | 51 | 25 | 26 | 5 | 39 | 12 | 0 | **conflicting** |
| Pleosporales | Ascochyta | 43 | 43 | 16 | 27 | 0 | 43 | 0 | 0 | confirmed |
| Botryosphaeriales | Neofusicoccum | 42 | 5 | 3 | 2 | 37 | 3 | 1 | 1 | mostly not called (call rate 12%) |
| Botryosphaeriales | Lasiodiplodia | 39 | 36 | 17 | 19 | 3 | 35 | 0 | 1 | confirmed |
| Mycosphaerellales | Cercospora | 36 | 35 | 14 | 21 | 1 | 31 | 3 | 1 | mixed |
| Mycosphaerellales | Fulvia | 31 | 31 | 20 | 11 | 0 | 31 | 0 | 0 | confirmed |
| Pleosporales | Didymella | 30 | 29 | 21 | 8 | 1 | 15 | 14 | 0 | **conflicting** |
| Mycosphaerellales | Sphaerulina | 28 | 8 | 2 | 6 | 20 | 8 | 0 | 0 | confirmed, call rate 29% |
| Mycosphaerellales | Friedmanniomyces | 26 | 26 | 21 | 5 | 0 | 1 | 24 | 1 | **conflicting** |
| Myriangiales | Elsinoe | 26 | 26 | 7 | 19 | 0 | 22 | 0 | 4 | mixed |
| Pleosporales | Curvularia | 24 | 1 | 0 | 1 | 23 | 1 | 0 | 0 | confirmed, call rate 4% |
| Botryosphaeriales | Diplodia | 21 | 19 | 13 | 6 | 2 | 19 | 0 | 0 | confirmed |
| Mycosphaerellales | Dothistroma | 20 | 20 | 6 | 14 | 0 | 20 | 0 | 0 | confirmed |
| Mycosphaerellales | Elasticomyces | 15 | 14 | 14 | 0 | 1 | 2 | 9 | 3 | **conflicting** |
| Pleosporales | Stemphylium | 14 | 5 | 5 | 0 | 9 | 0 | 5 | 0 | **conflicting**, call rate 36% |
| Pleosporales | Plenodomus | 14 | 14 | 6 | 8 | 0 | 13 | 1 | 0 | mixed |
| Mycosphaerellales | Teratosphaeria | 13 | 13 | 8 | 5 | 0 | 10 | 2 | 1 | mixed |
| Pleosporales | Epicoccum | 12 | 12 | 5 | 7 | 0 | 12 | 0 | 0 | confirmed |
| Mycosphaerellales | Pseudocercospora | 12 | 10 | 4 | 6 | 1 | 8 | 2 | 0 | **conflicting** |
| Pleosporales | Paraphaeosphaeria | 11 | 11 | 6 | 5 | 0 | 11 | 0 | 0 | confirmed |
| Pleosporales | Ophiosphaerella | 11 | 0 | 0 | 0 | 11 | 0 | 0 | 0 | no call |
| Pleosporales | Exserohilum | 10 | 0 | 0 | 0 | 10 | 0 | 0 | 0 | no call |

(The group label for Neofusicoccum in the generated table is "conflicting (low call rate)"; 1 of 5 called genomes is conflicting and
37 of 42 are not called, so I read it as a no-call group. Label rule below.) By family (taxonkit family, `tables/group_family.tsv`),
Saccotheciaceae, Teratosphaeriaceae, Phyllostictaceae and eight smaller families are conflicting, 48 are confirmed (243 genomes),
14 are mixed (2,056 genomes, the large families that contain both kinds of genus: Botryosphaeriaceae, Pleosporaceae,
Mycosphaerellaceae, Phaeosphaeriaceae, Didymellaceae, Cladosporiaceae).

How the calls changed since v0.6.0 (genomes present in both; the crosstab is derivable from `genome_audit_final.tsv` columns
`v060_idiomorphs` and `idiomorphs_called`): 1,307 MAT1-1 and 952 MAT1-2 unchanged; MAT1-1 to MAT1-2 in 89 (87 are the
*Parastagonospora* correction by the records, 2 others), MAT1-2 to MAT1-1 in 20; 47 genomes gained a call, 4 lost one; 17 gained a
second idiomorph. Against the 2026-10-05 with-records run 41 genomes gained a call (35 MAT1-1, 6 MAT1-2; 16 of them *Zymoseptoria*,
the polish-cap loss reversed), 1 lost one and 2 changed idiomorph, so the records and the later rules left the Dothideomycetes labels
almost unchanged.

## 2. What the calls rest on (measured)

- 2,592 called loci (single idiomorph): `mat_locus` 2,281, `partial_locus` 269, `idiomorph_gene_only` 42; confidence `high` 2,116,
  `medium` 473, `low` 3. No `homothallic_candidate` and no `two_idiomorphs` statement anywhere in the class.
- The called idiomorph's own core gene is present in 2,591 of 2,592, but **only 582 (22%) are full-length** by the pipeline's own
  measure (coverage of the reference 90% or more, modelled). 1,995 loci (77%) cover under 80% of their reference and 628 are under
  40% identical. The reason is the reference set: median own-gene coverage is 93% in Pleosporales (curated Bipolaris, Parastagonospora,
  Plenodomus records) and 41% in Botryosphaeriales, 39% in Dothideales, 31% in Cladosporiales and Myriangiales, where the
  nearest reference is far. 58 calls rest on an own gene under 80% covered and within 500 bp of a contig end (31 *Corynespora*,
  9 *Hortaea*, 6 *Aureobasidium*).
- **The `coverage_profile` partial-core policy never fired in this class** (0 loci carry an `idiomorph_vote_note`). Its trigger
  needs the winning core hit under 80% covered and the other idiomorph's core gene full-length (90%+); with distant references that
  combination does not occur. The policy neither helps nor harms Dothideomycetes; the Botrytis-type error is not what is wrong here.
- The vote was decided by an optional gene (`MAT1-1-4`, `MAT1-1-5`) in 5 loci (4 loci in 3 genomes of Cryomyces and Oleoguttula
  with MAT1-1-5 at 40 bits, and one *Botryosphaeria agaves* locus that has only `MAT1-1-4`, which is a call on no core gene at
  all). So the optional-gene vote problem is rare as a *vote* winner here.
- It is common as a *resolution* loser: in 65 loci a core gene was set aside in favour of an optional gene (the
  `MAT1-1-3 > MAT1-2-1` pair), 59 in *P. tritici-repentis*, 2 *Ramularia coccinea*, 1 each in *Ramulariopsis*, *Cladosporium
  austrohemisphaericum*, *Neofusicoccum cordaticola*, *Plenodomus libanotidis*. In another 1,475 loci the resolution went to the
  core gene, which is correct for the Parastagonospora case once the records exist.

## 3. Independent checks

### 3.1 Vote against the alpha-box / HMG-box profile
Over the 2,592 called loci the profile agrees with the vote in 2,445 (94.3%), is undetermined in 129 (both genes present or both weak)
and disagrees in 18. All 18 disagreements are loci where both an alpha-box and an HMG-box gene are present (Aureobasidium 13,
Recurvomyces, Meristemomyces, Elasticomyces, Catenulostroma, Setomelanomma): the profile prefers the HMG side, the vote the alpha side.
That is two methods disagreeing because the locus holds both genes, not one of them being wrong. Among loci with a determined profile
the agreement is 99.3%. Genus-out profiles (the five genera with a curated record) agree in 282 of 287 loci (5 undetermined,
all *Zymoseptoria*, none disagree); the order-out profile agrees in 862 of 878 Pleosporales loci (2 disagree, 14 undetermined) and 446 of 466
Mycosphaerellales loci (5 disagree, 15 undetermined). *Not independent evidence of truth*: the profile and the vote both derive from
the same curated records. It is evidence that the vote's gene identity matches the protein it picked.

### 3.2 Pfam on the called protein
MAT1-1 calls: the MAT1-1-1 protein has MATalpha_HMGbox (PF04769) in 1,390 of 1,402 primary calls (5 more by the scan ORF; 6 have
neither; 1 not scored). MAT1-2 calls: the MAT1-2-1 protein has an HMG-box domain in 813 of 1,072 (54 of them also have an alpha-box
protein at the locus); in 248 more the pipeline's gene model is a fragment with no Pfam hit and the scan finds an HMG-box ORF at
the locus; 11 have neither. The domain check is independent of the curated database and of the pipeline's search. It cannot separate MAT1-2-1
from MAT1-1-3 (both HMG), which matters for homothallic calls (below). Seventeen primary calls have no domain on the called gene
(8 *Corynespora cassiicola* MAT1-2 at 26% coverage, and single genomes of Teratosphaeria, Neocatenulostroma, Hortaea, Cercospora,
Trematosphaeria, Nothophaeocryptopus, Mycosphaerellales sp. BF13, Meristemomyces, Catenulostroma).

### 3.3 Genome-wide ORF scan: is the other idiomorph's gene at the locus?
For every genome with a call, the best alpha-box ORF (profile >= 40 bits) and HMG-box ORF within 20 kb of the called locus, whatever
the pipeline did with them. **Baseline**: in the 12 groups (11 genera and *Pyrenophora teres*) whose literature or sampling says
heterothallic or single-idiomorph (Alternaria, Parastagonospora, Venturia, Zymoseptoria, Ascochyta, Corynespora, Hortaea,
Cercospora, Fulvia, Dothistroma, Macrophomina, Pyrenophora teres), no MAT1-1 locus has an HMG ORF of 40 bits within 20 kb and no
MAT1-2 locus has an alpha-box ORF there: 0 of 1,400 primary calls (715 MAT1-1, 685 MAT1-2), all scanned. That zero is what makes the
positives interpretable:

| group | calls scanned | HMG or alpha-box ORF of the other kind at the locus | strong (>= 90 bits) |
|---|---|---|---|
| Aureobasidium | 180 | 145 (81%) | 77 |
| Pyrenophora tritici-repentis | 64 | 62 | 62 |
| Botryosphaeria dothidea | 140 | 137 | 0 (all 69 bits) |
| Friedmanniomyces | 26 | 23 | 2 |
| Didymella pinodes | 8 | 8 | 8 |
| Stemphylium (5 called of 14) | 5 | 5 | 5 |
| Phyllosticta capitalensis | 8 | 8 | 8 |
| Elasticomyces | 14 | 12 | 3 |
| the 12 baseline groups (P. teres, Parastagonospora, Alternaria, Venturia, Zymoseptoria, Macrophomina, Corynespora, Ascochyta, Hortaea, Fulvia, Cercospora, Dothistroma) | 1,400 | 0 | 0 |

The pipeline's Pfam-based flag and the scan agree on direction: of 2,474 primary calls the two give "strong" together in 95 genomes,
the scan alone in 82 more and Pfam alone in 5 (`second_gene_pfam` against `second_gene_scan` in `genome_audit_final.tsv`).

### 3.4 Independent truth (what exists)
- **BioSample `mating_type`**: 2 of 2,731 genomes (`biosample/`: the BioSample was fetched for 2,604 genomes, 127 had none; a text
  scan of every attribute, title and description for mating-type wording found no other statement; the `propagation:
  asexual/sexual` hits are life-cycle notes, not mating types). *Ascochyta rabiei* ArDII stated MAT1-2, called MAT1-2 (high; profile
  MAT1-2, Pfam HMG yes). *Venturia effusa* albino stated MAT1-1, called MAT1-1 (medium; profile and Pfam agree). 2 of 2, no inference
  possible from n = 2. The submitter's statement is not verified.
- **Strain-matched literature, in the curated records** (not independent of the vote for these genomes, because the record is a
  reference): *Zymoseptoria tritici* IPO323 is MAT1-1 (PubMed 11929216); the genome `GCF_000219625.1` is called MAT1-1 (`partial_locus`,
  medium) and the Zymoseptoria-out profile also says MAT1-1. The four *Bipolaris maydis* genomes of strains C4 (MAT-2) and C5 (MAT-1)
  from the original Turgeon cloning (PubMed 8479433) are **not detected**.
- **Population ratio** for species with at least 10 called genomes (`tables/species_ratio.tsv`, exact binomial against 1:1; heterothallic
  samples should be near 1:1 unless the sample is clonal). Balanced: *Parastagonospora nodorum* 87:88, *Zymoseptoria tritici* 33:28,
  *Macrophomina phaseolina* 237:228, *Hortaea werneckii* 36:31, *Corynespora cassiicola* 27:30, *Venturia inaequalis* 36:51,
  *Aureobasidium pullulans* 37:46, *A. melanogenum* 32:27. Skewed: *Alternaria alternata* 99:46 (p = 1e-5), *Phyllosticta citricarpa* 5:20,
  *Ascochyta rabiei* 0:16, *Friedmanniomyces* 21:4, *Botryosphaeria dothidea* 140:0, *Pyrenophora tritici-repentis* 64:0, *Elasticomyces* 14:0.
  `species_skew_by_bioproject.tsv` shows that the *Alternaria alternata* skew sits in one BioProject (PRJNA770139, 52:9 of 61; the other
  BioProject with 29 genomes is 15:14), *P. citricarpa* in PRJNA636278 (15:0), *A. rabiei* in PRJNA729126 (12:0), *Friedmanniomyces* in
  PRJNA666634 (20:4 of 24) and *B. dothidea* in PRJNA777748 (126:0). So the skews of the heterothallic species are sampling. The two
  that are not sampling are the homothallic ones (*P. tritici-repentis* 4 BioProjects, 0 MAT1-2 in each; *Didymella pinodes* 8:0 in 2).
  Replacing the Parastagonospora 178:4 of v0.6.0 by 92:92 is the same test, passed.

### 3.5 Literature (from PubMed abstracts; DOIs in `literature.tsv`)
According to PubMed: *Pyrenophora tritici-repentis* is self-fertile, each isolate carries both idiomorphs and deleting either gene
abolishes pseudothecia (10.1016/j.fgb.2017.10.004; PMID 20373969 for the homothallic nature). All four sequenced *Aureobasidium pullulans*
varieties contain a homothallic MAT locus (10.1186/1471-2164-15-549). *Didymella clematidis* and *Peyronellaea (Didymella) pinodes* are
homothallic by a crossover (and an inversion) of MAT1-1 and MAT1-2 and *Pey. pinodella* is heterothallic (10.1111/j.1364-3703.2011.00751.x).
*Eutiarosporella darliae* and *E. pseudodarliae* carry both genes in one genome (10.1016/j.fgb.2017.10.005). Botryosphaeriaceae: nine
species homothallic and seven heterothallic, with MAT fragments flanking the loci (10.1016/j.fgb.2018.03.003; the abstract does not say
which are which). *Parastagonospora (Phaeosphaeria) nodorum* heterothallic with both mating types in field populations
(10.1016/s1087-1845(03)00062-8). *Zymoseptoria (Mycosphaerella graminicola)* heterothallic, equal allele frequencies
(10.1006/fgbi.2001.1322). *Cercospora beticola* is heterothallic but every isolate carries exon fragments of both MAT1-1-1 and MAT1-2-1
at other loci (10.1016/j.fgb.2013.10.011). *Cladosporium (Fulvia) fulvum* has MAT1-1-1 or MAT1-2-1 per strain
(10.1016/j.fgb.2006.11.004). *Venturia effusa* 1:1 MAT1-1:MAT1-2 (10.1094/PHYTO-07-17-0233-R); *Diplodia pinea* heterothallic
(10.1016/j.fgb.2013.10.013); *Leptosphaeria maculans* (10.1007/s00294-003-0391-6); *Cochliobolus heterostrophus* heterothallic with
homothallic Cochliobolus species carrying both genes (10.1007/BF00279556, 10.1073/pnas.96.10.5592); *Alternaria alternata* populations carry one idiomorph
per isolate (PMIDs 26666175, 23441968, abstracts only).
Links (PubMed, abstracts read): [20373969 homothallic P. tritici-repentis](https://pubmed.ncbi.nlm.nih.gov/20373969/),
[10.1016/j.fgb.2017.10.004](https://doi.org/10.1016/j.fgb.2017.10.004), [10.1186/1471-2164-15-549](https://doi.org/10.1186/1471-2164-15-549),
[10.1111/j.1364-3703.2011.00751.x](https://doi.org/10.1111/j.1364-3703.2011.00751.x), [10.1016/j.fgb.2017.10.005](https://doi.org/10.1016/j.fgb.2017.10.005),
[10.1016/j.fgb.2018.03.003](https://doi.org/10.1016/j.fgb.2018.03.003), [10.1016/s1087-1845(03)00062-8](https://doi.org/10.1016/s1087-1845(03)00062-8),
[10.1006/fgbi.2001.1322](https://doi.org/10.1006/fgbi.2001.1322), [10.1016/j.fgb.2013.10.011](https://doi.org/10.1016/j.fgb.2013.10.011),
[10.1016/j.fgb.2006.11.004](https://doi.org/10.1016/j.fgb.2006.11.004), [10.1094/PHYTO-07-17-0233-R](https://doi.org/10.1094/PHYTO-07-17-0233-R),
[10.1016/j.fgb.2013.10.013](https://doi.org/10.1016/j.fgb.2013.10.013), [10.1007/s00294-003-0391-6](https://doi.org/10.1007/s00294-003-0391-6),
[10.1007/BF00279556](https://doi.org/10.1007/BF00279556), [10.1073/pnas.96.10.5592](https://doi.org/10.1073/pnas.96.10.5592),
[26666175 Alternaria alternata cryptic sex](https://pubmed.ncbi.nlm.nih.gov/26666175/), [23441968](https://pubmed.ncbi.nlm.nih.gov/23441968/).

**From memory or from the repo notes, not re-opened:** *Stemphylium* homothallic species (Inderbitzin 2005, in
`docs/notes/2026-09-24_mat-reference-gap-literature.md`), *Pyrenophora teres* heterothallic (Rau 2005, same note), *Phyllosticta
capitalensis* both genes in one deposit KT708824.1 (same note, unpublished GenBank record), *Setosphaeria turcica* heterothallic.
**Searched, nothing usable:** *Botryosphaeria dothidea* and *Macrophomina phaseolina* mating system (two PubMed searches; the 2017
pathogen-profile review has no accessible full text). **Not searched:** Friedmanniomyces, Elasticomyces, Cryomyces, Corynespora, Hortaea,
Elsinoe, Cladosporium s.s. mating systems.

### 3.6 HMG-box paralogs
For the 1,072 MAT1-2 primary calls, an off-locus HMG ORF scores at least 10 bits higher on the MAT1-2-1 profile than the one at the
called locus in 82 (17 by 30 bits); it concentrates in Dothistroma (14 of 14), Fulvia (11 of 11), Lasiodiplodia (18 of 19) where
the locus ORF scores 64 to 73 (medians) and an off-locus ORF 85 to 96 (`tables/hmg_paralog_check.tsv`). These are exons scored by a
full-length profile, so I do not read them as wrong calls; they show that a profile score alone cannot tell the MAT1-2-1 HMG ORF
from a paralog, which is why the verdict never uses it alone. Eight MAT1-2 loci have no HMG ORF of 40 bits at the locus at all.

## 4. Classification

Per genome, primary locus = highest confidence, then locus class, then bitscore (`finalize.py`):
- **conflicting**: the profile disagrees with the vote; or a strong second-gene signal (alpha-box Pfam on a MAT1-2 locus with HMG,
  or an HMG-box ORF of 150+ bits by Pfam/profile or 90+ bits by the scan at a MAT1-1 locus) sits at the locus; or loci of both
  idiomorphs are called; or a core gene was set aside for an optional gene.
- **unverifiable**: a weak second-gene signal (HMG ORF 40 to 90 bits at a MAT1-1 locus, below the MAT1-2-1 range), no domain on the
  called gene and none by scan, or the profile undetermined or not scored.
- **confirmed**: profile agrees and a domain is there (Pfam, or the scan ORF at 60+ bits where the pipeline's gene model is a
  fragment) and nothing above flags.
- Group (genus, family, order) label over genomes with a call: confirmed if 90% confirmed and none conflicting; conflicting if 20% or
  more conflicting; unverifiable if under 50% confirmed and under 20% conflicting; mixed otherwise; "low call rate" when under half
  of the group's genomes have a call. The thresholds are mine (judgement); the counts are in the tables.

Results over 2,730 genomes: 1,969 confirmed, 298 conflicting, 207 unverifiable, 256 no call. By genus: 138 genera (1,658 genomes)
confirmed, 37 conflicting (543 genomes), 10 mixed, 9 unverifiable, 15 no-call (36 genomes). Confirmed calls are mostly `high`
(1,682) but 285 are `medium`; conflicting calls are 183 `high` and 115 `medium`, so the confidence field does not separate them.
Of the 298 conflicting: second gene at locus alone 169, plus a core set aside for an optional gene 60, two idiomorph loci 25, second
gene plus two loci 22, profile disagrees plus second gene 16, core set aside alone 5, profile disagrees alone 1. Of the 207
unverifiable: weak second gene alone 186 (137 *B. dothidea*, 29 *Aureobasidium*, 20 others), profile undetermined alone 8,
domain absent with profile undetermined 8, domain absent alone 3, weak second gene with another reason 2.

## 5. Suspect groups and the evidence (`tables/suspect_groups.tsv`, `tables/note_suspects.md`)

| group | genomes (called) | verdict | evidence |
|---|---|---|---|
| Aureobasidium, 187 (180) | A. pullulans 37:46, A. melanogenum 32:27, A. subglaciale 9:0, others 27:2 | **wrong or arbitrary**: homothallic | complete alpha-box and HMG-box genes in one locus in 145 of 180 (scan), Pfam alpha and HMG on one locus in 91; profile undetermined in 71 and disagreeing in 13 because both genes are there; 15 loci called MAT1-2 have `MAT1-2-1 + APN2 + MAT1-1-4` and no MAT1-1-1 found (partial); literature homothallic (10.1186/1471-2164-15-549) |
| *Pyrenophora tritici-repentis*, 64 (64) | 64 MAT1-1, 0 MAT1-2 | **wrong**: homothallic, real MAT1-2-1 set aside | `MAT1-2-1`/`MAT1-1-3` resolved against the HMG gene in 59 (2 `high`, 62 `medium`); HMG ORF at the locus (90+ bits) in 62 of 64; literature self-fertile. *P. teres* in the same genus, 11:13, no second gene, all 24 confirmed |
| Didymella pinodes 8 (8), other Didymella | 8:0; other Didymella 5 of 13 strong | **wrong** for *D. pinodes* (homothallic, PMID 22014305); other species 5 of 13 with strong second gene | scan HMG ORF at locus 8 of 8 (*pinodes*), 1 of 8 (*pinodella*, heterothallic, correct) |
| Friedmanniomyces 26, Elasticomyces 15, Cryomyces/Recurvomyces/Meristemomyces 13 | 21:4, 14:0, 4:6 | **likely wrong**: two loci or both genes | 26 genomes carry loci of both idiomorphs (Friedmanniomyces 20 of 26, Elasticomyces 5 of 15, Cryomyces 1); 24 and 9 conflicting; no literature found; black rock fungi, mostly one BioProject (PRJNA666634, 24 Friedmanniomyces and 14 Elasticomyces genomes) |
| Phyllosticta capitalensis 9 (8) | 7:1 | **wrong** (both genes; the repo note says homothallic, unpublished) | scan 8 of 8 strong, profile undetermined 8 of 8; *P. citricarpa* (25 called) is clean: 5:20, 25 confirmed |
| Stemphylium 14 (5) | 5:0 | wrong for the 5 called; 9 not called | strong second gene 5 of 5; homothallic from memory/repo note only |
| Botryosphaeria dothidea 141 (140) | 140:0 | **cannot be judged** | alpha-box real; one 238 aa HMG ORF of 69 bits at the locus in 137 of 140; the pipeline's own MAT1-2-1 hit is 33% identical over 27% of the reference; both a MAT1-1-3-like gene and a homothallic MAT1-2-1 fit; one BioProject carries 126 genomes |
| Macrophomina 491 (490) | 251:239 | **probably right, thin evidence** | MAT1-2 gene model 11% coverage, 39% identity, 231 of 239 with no Pfam hit; scan: HMG ORF about 100 bits at the locus in 239 of 239, no alpha-box ORF in the genome; MAT1-1: alpha-box Pfam 109 bits, 251 of 251; balanced ratio |
| Cercospora 36 (35) | 14:21; 3 genomes with two loci | MAT1-2 call in *C. kikuchii* **probably wrong** | *C. kikuchii* (3 genomes): MAT1-2 `high` at 16% coverage, 71% identity, no domain; plus a real MAT1-1-1 locus; matches the off-locus MAT fragments of 10.1016/j.fgb.2013.10.011 |
| 47 genomes with two loci, other genera | Cladosporium 3, Teratosphaeria 2, Baudoinia, Hortaea, Pseudocercospora, ... | **likely fragment calls** | in 43 of 47 one locus has an own-gene coverage under 35%, in 15 under 25%; off-locus fragments (as in Cercospora) are the working explanation |
| Corynespora 58 (58) | 27:31 | right, 8 unverifiable | 8 MAT1-2 at 26% coverage with no domain (a fragment at a contig end: 31 loci), the rest confirmed |
| Bipolaris, Curvularia, Exserohilum, Ophiosphaerella, Neofusicoccum, Sphaerulina, Stemphylium | 1 of 58, 1 of 24, 0 of 10, 0 of 11, 5 of 42, 8 of 28, 5 of 14 | **cannot be judged** (no call) | the C5 report has 35 suppressed unpolished loci; the `MAT1-1-1 + SLA2` candidate on KB445578.1 is `withheld: modelled_gene_bar`; not diagnosed further |

Right and well supported: Parastagonospora (184 of 184 confirmed, 92:92, v0.6.0 was 178:4), Alternaria (265 of 267 confirmed),
Venturia (104 of 104), Zymoseptoria (71 confirmed, 5 profile-undetermined), Ascochyta (43), Fulvia (31), Lasiodiplodia (35 of 36),
Diplodia (19), Dothistroma (20), Epicoccum, Paraphaeosphaeria, Coniothyrium, Laburnicola, Neoscytalidium (9:0, 6 of the 9 from
one BioProject), Pyrenophora teres. These are right in the sense that the gene named at the locus is the idiomorph's own gene by profile, domain and
scan, with nothing of the other idiomorph beside it. For *Alternaria alternata* the 99:46 skew is one BioProject; nothing in the audit
flags the calls.

## 6. What cannot be judged

- Every genome where the call rests on one gene with a distant reference and no profile margin (the 207 unverifiable).
- All no-call groups (256 genomes): recall, not idiomorph.
- Whether *Botryosphaeria dothidea* is homothallic (needs the annotated Botdo1 or the Nagel 2018 tables).
- Whether the 40 to 90 bit HMG ORFs at Aureobasidium (29), Friedmanniomyces (21), Elasticomyces (9) loci are MAT1-2-1 exons or MAT1-1-3-like genes.
- Any species with no literature and no BioSample statement, which is almost all of them (only 2 genomes have a BioSample statement).

## 7. Limits

- **Independence.** The profile and the vote share the curated records, so profile agreement is a consistency check. Pfam and the scan do
  not use the records, but the scan scores with the shipped profiles (including, for the five reference genera, the genus's own
  records), and Pfam GA cutoffs and my scan thresholds (40, 60, 90, 150 bits) were set by looking at these data: the baseline
  of 0 positives in heterothallic genera is what calibrates them, there is no held-out set.
- **The baseline assumes those genera are heterothallic** (literature for Parastagonospora, Zymoseptoria, Venturia, Ascochyta,
  Cercospora, Fulvia, Alternaria alternata in the populations cited; assumed for Hortaea, Corynespora, Macrophomina, Dothistroma).
  A real MAT1-1-3 HMG gene beside MAT1-1-1 would also count as a second gene; none appeared in the baseline groups.
- **Genomes are not independent.** 126 of 140 *B. dothidea*, 39 of 64 *P. tritici-repentis*, 61 of 145 *A. alternata* come from one
  BioProject each, so a group's genome count overstates the evidence; every count in this note is a count of assemblies.
- The verdict is per primary locus; genomes with several loci are summarised by the best one plus the two-idiomorph flag.
- Scan windows: 20 kb either side of the pipeline's locus; ORFs are stop-to-stop, so a gene with introns appears as several
  exon-sized ORFs (a 69 bit HMG ORF can be one exon of a gene that scores 150 whole).
- The reference set is uneven (five Dothideomycetes genera, two orders). The order-out profile (all curated records of the order
  withheld) was computed for Pleosporales and Mycosphaerellales only and agrees with the vote in 97%, but it was not used in
  the verdict, so Botryosphaeriales, Dothideales, Venturiales and the rest are scored by profiles that have no close relative.
- BioSample truth is the submitter's statement, 2 genomes. Strain-matched literature statements were chosen by name (C4, C5, IPO323)
  and not cross-checked against the assemblies' provenance.
- The v0.6.0 and records-run comparisons use the campaign's own labels; v0.6.0 vs now mixes the effect of records, the
  partial-core policy and the other merges since, which were not separated here.
- Two genomes have no report; 256 have no call.

## 8. Open questions (for the curator; no code proposed)

1. **Homothallic loci.** Should a Dothideomycetes locus with a complete alpha-box gene and an HMG ORF at the same locus be reported
   as one idiomorph at all? Today `homothallic_candidate` and `two_idiomorphs` never appear (0 of 2,592 and 0 of 2,728) while 145
   Aureobasidium, 62 *P. tritici-repentis*, 14 Didymella, 23 Friedmanniomyces genomes show the pattern at the locus (scan). The tracker already notes the
   exclusion in `full_length_models` that blocks `homothallic_candidate`.
2. **MAT1-1-3 supersedes MAT1-2-1** in 65 loci, 59 of them one species; the optional-gene note and the dropped PR concerned exactly
   this resolution step. Does the curator want the Pyrenophora result counted as a case for it, given the campaign found
   `core_first`/`resolution_only` changes unverified elsewhere?
3. **Botryosphaeria dothidea** (140 genomes, 126 one BioProject): homothallic or MAT1-1 with an HMG-box companion? Needs the annotated
   genome (JGI Botdo1) or the Nagel 2018 supplement.
4. **Cochliobolus clade.** Why are Bipolaris/Curvularia/Exserohilum not called even for C4 and C5, the isolates of the records
   (35 suppressed loci, `modelled_gene_bar`)? Not diagnosed in this audit.
5. **Truth.** Which Dothideomycetes species could give strain-level mating types (Z. tritici populations, Parastagonospora
   populations, A. rabiei, C. beticola, V. inaequalis) for a real held-out test? BioSample has 2 statements for 2,731 genomes.
6. **Macrophomina MAT1-2**: is a call on an 11%-coverage HMG fragment acceptable at `high` confidence? The scan supports it, the
   pipeline's own evidence does not.
7. **Off-locus fragments** (43 of the 47 two-idiomorph genomes have a locus under 35% own-gene coverage, Cercospora among them): should a
   second locus of that kind be allowed to add an idiomorph?
8. Are the 129 `profile_undetermined` calls (71 Aureobasidium) a feature of the 25-bit margin, which was borrowed from Mucoromycota and
   not calibrated for Ascomycota?

## 9. Files (`results/2026-10-18_dothideomycetes_audit/`)

`genome_list.tsv` (2,731), `genome_status.tsv`, `loci_audit.tsv` (raw per-locus extraction), `loci_audit_flagged.tsv` (plus flags),
`genome_audit_final.tsv` (one row per genome with the verdict), `scan/` (genome-wide ORF scan), `reports_all.tar.zst` (all 2,728
detection reports), `bioproject.tsv`, `biosample/` (text scan of BioSample), `v060/` (v0.6.0 and records-run per genome),
`literature.tsv`, `tables/` (`by_order`, `by_genus`, `by_family`, `group_*`, `suspect_groups`, `species_ratio`,
`species_skew_by_bioproject`, `genomewide_*`, `hmg_paralog_check`, `truth_comparison`, `summary.json`, `note_*.md`), scripts
(`make_list.py`, `audit_worker.slurm`, `audit_extract.py`, `scan_genome.py`, `scan_worker.slurm`, `analyze.py`, `tables.py`,
`truth.py`, `scan_integrate.py`, `finalize.py`, `suspects.py`, `paralog_check.py`, `species_skew_by_bioproject.py`, `note_tables.py`,
`run_all.sh`; `adhoc/` holds the one-off queries quoted above). Slurm jobs: taxonkit 30049146 and 30050120, run 30049176, Pfam subset
30049468, extraction 30056459, scan 30053543 (chunks 7 and 23 resubmitted, jobs 30053090 range, 30056836, after a >100,000-residue ORF
crashed the first pass; resumable). Not committed: the per-genome `evidence_diagnostics.jsonl` and logs on the HPCC under
`.../run-c9c0ec6/results/2026-10-18_dothideomycetes_audit/`.
