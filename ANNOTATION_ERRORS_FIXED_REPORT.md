# Annotation errors found and fixed

A running record of every error found in **gene assignment**, **gene models**
or **protein matches**, whether in this project's curated database, in
deposited sequence records (GenBank / RefSeq), or in the literature's
metadata. The purpose is to know what needs updating upstream and what our
records work around.

Started 2026-09-25 at the curator's request. Add to it whenever a record is
curated or corrected. Newest entries at the end of each section.

## How to read an entry

**Evidence codes**

* **V** -- verified in this session by a measurement we can re-run (the
  evidence path is given).
* **A** -- measured by a literature/database agent and reported in a notes
  file, not re-checked independently. Treat as a lead until re-checked.
* **M** -- from a prior session's curation notes (project memory), not
  re-checked on 2026-09-25.

**Status**

* **FIXED** -- corrected in this project's database or docs (commit given).
* **WORKED AROUND** -- the upstream record is still wrong; our record avoids
  the error (and says so in its definition_note).
* **OPEN (upstream)** -- the deposited record should be updated; we have not
  contacted NCBI or the submitters.
* **OPEN (ruling)** -- needs a curator decision before we act.

---

## A. Errors in this project's curated database

| # | record | error | ev. | status |
|---|---|---|---|---|
| A1 | `4959_cbs767_MTL_A` (*D. hansenii* CBS767) | Said "no MTLalpha1 gene is present in this strain's assembly" and `system: heterothallic`. XP_460134.1 at NC_006047.2:1,587,954-1,588,592, ~1.5 kb from MTLa1, is MTLalpha1: the only PF04769 protein in the proteome, reciprocal best hit of *C. albicans* MTLalpha1 (E=2.8e-25), CTG-clade placement in a gene tree, inside the MTL locus. | V `results/2026-09-24_dhansenii_alpha1/EVIDENCE.md` | FIXED `ca619b6` (record_version 2; homothallic candidate) |
| A2 | `28985_nrrl-y-1140_MATsc_MATa` (*K. lactis*) | Evidence block copied from the *S. cerevisiae* S288C record: cited only Astell 1981 plus Kostriken 1983 (the *S. cerevisiae* HO endonuclease; *K. lactis* switches through alpha3), and a boundaries sentence about S288C coordinates on NC_001135.5. | V | FIXED `25f47ba`: Astell 1981 kept as structural basis; Astrom 2000 (PMID 10978277), Butler 2004 (PMID 14745027) added |
| A3 | `28985_nrrl-y-1140_MATsc_HMLalpha` (*K. lactis*) | Same copied block. | V | FIXED `25f47ba`, plus Barsoum 2010 (PMID 20008928) |
| A4 | `381046_cbs-6340_MATsc_HMLalpha` (*L. thermotolerans*) | Same copied block, and it claimed tier 1 although the idiomorph and boundaries come from genome sequence only. | V | FIXED `25f47ba`: Souciet 2009 (PMID 19525356); evidence now tier 2 |
| A5 | `order.yml`, `MTL` roster | `MTLalpha2` has no roster slot, so every hit to it is silently dropped. It exists in *C. albicans*, *C. dubliniensis* and *C. tropicalis* but is biologically absent in Metschnikowiaceae (Munoz et al. 2018), so adding it naively depresses scores there. | M, V (deposits AF167163.1, AY622606.1) | OPEN (ruling): add as `optional`, or add a per-clade restriction |
| A6 | `order.yml`, `MTL` gene_class | We set `MTLA2: HMG_box` on 2026-09-25 from the literature. Measured afterwards: Pfam HMG_box (PF00505) does not hit the curated *C. albicans* or *C. lusitaniae* MTLa2 proteins. The class rests on the literature, not a domain hit. | V `results/2026-09-25_cauris_mtl/` | FIXED (documented) `8dec9e1` |

## B. Errors in deposited gene models and annotations

### B1. Core MAT genes present in the genome but NOT annotated

| # | assembly / record | gene | where it actually is | ev. | status |
|---|---|---|---|---|---|
| B1.1 | *C. auris* B8441, RefSeq `GCF_002759435.1` | MTLa2 and MTLa1 | NC_140807.1, between PIK1 (ends 1,494,353) and B9J08_03706 (1,495,702): tblastn *C. lusitaniae* a2 at 1,494,442-1,495,046 (two HSPs, E=4.7e-11), *C. albicans* a1 at 1,495,483-1,495,620 (E=2.7e-5). exonerate and miniprot build no model even with relaxed thresholds (30-43% identity). | V `results/2026-09-25_cauris_mtl/` | WORKED AROUND (curated from other clades' annotations, B2.2/B2.3); OPEN (upstream) |
| B1.2 | *C. auris* B11205, `GCA_016772135.1` | MTLa2 | Not annotated; the ~780 bp between MTLa1 (WZC25182.1) and PIK1 is where it would sit. Not modelled. | V | OPEN (upstream) |
| B1.3 | *C. tropicalis* MYA-3404 scaffold `GG692408.1` | MTLa2 | tblastn *C. albicans* MTLa2: 7,360-7,638 and 7,692-8,021, E=4e-41. | A round-2 notes | OPEN (upstream); not curated |
| B1.4 | *C. tropicalis* MYA-3404 `GG692402.1` | MTLalpha2 | tblastn *C. albicans* alpha2 at 185,474-186,082, 45% id, E=6e-47. | A | OPEN (upstream) |
| B1.5 | *L. thermotolerans* CBS 6340 `NC_013082.1` | MATa1 at the MAT position | ORF at 282,329-282,664 (-), 111 aa, unannotated. The only accessioned a1 (XP_002554220.1, silent cassette) lacks the homeodomain because a1's 3' end runs from the Y box into a non-homologous X box -- real switching biology, not an annotation error. | M | OPEN (upstream); curation deferred by curator 2026-09-21 |
| B1.6 | *Mycosarcoma (Ustilago) maydis* RefSeq `NC_026482.1` | mfa1 (41 aa pheromone precursor) | Absent from the modern annotation; present in the 1995 locus deposit U37795. A curation fork had relabelled Rba1 as "mfa". | M (annotation-gap lesson) | WORKED AROUND (curated from U37795) |

### B2. Gene models that are truncated or split

| # | record | error | ev. | status |
|---|---|---|---|---|
| B2.1 | *C. tropicalis* MYA-3404 `EER30103.1` ("Mtla1p", 148 aa) | Lacks MTLa1's first exon: *C. albicans* a1 aligns at 8,179-8,367 (aa 1-66) and 8,421-8,912 (aa 64-206). Do not curate EER30103.1 as the a1 protein. | A | OPEN (upstream) |
| B2.2 | *C. auris* B11243 `PSK75932.1` ("mating_type_MTLa1", 99 aa) | A single-exon model with the homeodomain at residues 46-92, against 166 aa and residues 99-154 in the three-exon clade I models (e.g. WZC25182.1). The 5' exons are missing. | V | WORKED AROUND: curated as `completeness: partial` in `498019_b11243_MTL_A` (`8dec9e1`); the full-length a1 curated separately from B11205 |
| B2.3 | *C. albicans* WO-1 MTLalpha2 model | Missed a 59 bp intron: 167 aa instead of 186. | M | WORKED AROUND (not curated; SC5314 deposit AF167163.1 used) |
| B2.4 | *M. importuna* `KY782629.1` / `KY782630.1` | APN2 is annotated as TWO adjacent CDS in each deposit (AVI60802.1 + AVI60803.1; AVI60822.1 + AVI60823.1), almost certainly one gene split by the annotation. | V | WORKED AROUND: both kept under one gene name APN2 (`81ea22a`) |
| B2.5 | *L. thermotolerans* MATALPHA2 | 108 aa over two exons against 210 aa (*S. cerevisiae*) and 223 aa (*K. lactis*); possibly N-terminally truncated. No continuous alternative ORF longer than 63 aa. | M (record note) | OPEN: treat as the record's weak gene |

### B3. Genes annotated without their name, or with a wrong/uninformative product

| # | record | error | ev. | status |
|---|---|---|---|---|
| B3.1 | *D. hansenii* CBS767 `XP_460134.1` | MTLalpha1 annotated only as "DEHA2E19096p"; a name search misses it (this is what caused A1). | V | OPEN (upstream) |
| B3.2 | *C. albicans* SC5314 `AF167163.1` | Every CDS has product "unknown"; the alpha2 CDS (AAD51408.1) also has no /gene name. | V | WORKED AROUND (named in our record) |
| B3.3 | *C. dubliniensis* `AY622606.1` | OBP product spelled "OPB alpha"; PAPalpha CDS 3'-partial; no /gene qualifiers. | V | WORKED AROUND |
| B3.4 | *C. auris* B11205 `WZC25182.1`, and FDK38_003635 `QRG39207.1` | Full-length MTLa1 (166 aa, homeodomain E=4.8e-19, next to PIK1) annotated only as "hypothetical protein". | V | WORKED AROUND (curated `498019_b11205_MTL_A`); OPEN (upstream) |
| B3.5 | *K. lactis* NRRL Y-1140 RefSeq | The MAT-position MATa CDSs are unnamed; idiomorph was established by blastn of the AF195067.1 cassette. | M (record note) | WORKED AROUND |
| B3.6 | *M. importuna* `KY782629.1` | SLA2 annotated as "Endocytosis protein end4" (the *S. pombe* name). | V | WORKED AROUND (named SLA2) |
| B3.7 | Mortierellomycota NCBI proteins, e.g. GJJ68563.1 "mating-type protein A1", OAQ36521.1 "mating type protein 2, partial" | Names assigned by similarity pipelines, not characterized MAT genes. | A | OPEN: do not use as references |

## C. Errors in deposit metadata or literature statements

| # | record / source | error | ev. | status |
|---|---|---|---|---|
| C1 | *Y. lipolytica* `AJ617307.1` | Declares `strain=W29`, but the W29 reference genome contains no matb sequence (best tblastn E=0.82; matb1 no hit at E<=10). | M | WORKED AROUND (strain attribution not copied); OPEN (upstream) |
| C2 | *Lobaria pulmonaria* `JX520967.1` / `JX520966.1`; *Microbotryum* `JQ423666.1`, `JQ423663.1`, `JQ423661.1`; `KX832965.1`; `ON315861.1` | Flagged UNVERIFIED in GenBank; several carry no CDS. | A | Not usable as references |
| C3 | *P. tritici-repentis* `AM884596`-`AM884619` | Titled "mat1-1 gene and partial mat1-2 gene" -- unusual for a heterothallic MAT1-1 deposit. | A | OPEN: check before use |
| C4 | *Morchella* literature | The Morchella locus deposits are from Chai et al. 2017 (Mycol Prog 16:743, doi:10.1007/s11557-017-1309-x), not Du et al. as first recalled. | A, V (Crossref) | FIXED in our record |

## D. Assignment traps: real biology that looks like an annotation error

Not errors, but each causes a wrong call if the gene is taken at face value.
Recorded so a detector rule or a curator does not "fix" them.

| # | lineage | trap | ev. | consequence |
|---|---|---|---|---|
| D1 | Sordariomycetes, Leotiomyceta | MAT1-1-3 (HMG box) in the MAT1-1 idiomorph cross-matches MAT1-2-1. | V (panels; resolution events "MAT1-2-1 beat MAT1-1-3") | Looks like both idiomorphs. Excluded from the homothallic rule (`338ee18`). |
| D2 | Ophiostomatales | MAT1-2 idiomorphs carry a truncated MAT1-1-1 remnant (e.g. a 266 bp fragment ~1.2 kb from MAT1-2-1 in *Leptographium procerum*). | V | Looks like both idiomorphs. Excluded by the full-length (>= 50%) rule (`338ee18`). |
| D3 | *Colletotrichum* / *Glomerella* | The MAT1-2 HMG box is present in BOTH partners of fertile crosses; MAT1-1-1 is not the determinant (Menat et al. 2012, PMID 22223174). | A | A MAT1-1-1 vs MAT1-2-1 idiomorph call is wrong in this genus. |
| D4 | CTG clade (Metschnikowiaceae) | MTLalpha2 is biologically absent (Munoz et al. 2018). | M | Not a missing annotation. |
| D5 | *C. parapsilosis* | MTLa1 is a pseudogene (Logue et al. 2005, PMID 15947193). | A | A missing a1 is expected. |
| D6 | *L. elongisporus*, *C. sojae* | No MAT genes at all. | A | A no-call is correct. |
| D7 | *Hanseniaspora* | MAT genes are present (Krassowski et al. 2019); the claim of MAT loss is not supported. | A, V (synteny: 15/18 no-call genomes carry a core MAT hit next to SLA2) | A no-call is a miss. |
| D8 | *Heterobasidion*, *Microbotryum* | a2 / HD pseudogenes in some strains or species. | A | A disrupted HD is not a detection failure. |

## E. Errors in this project's own documents

| # | document | error | ev. | status |
|---|---|---|---|---|
| E1 | `docs/decoy-feasibility.md` | Listed PF08800 as "MATA_HMG". PF08800 is "BT4734-like_N". | V (InterPro HMM NAME line) | FIXED (row corrected with a note; count not re-measured) |
| E2 | `docs/holdout-benchmark.md` | Recall table scored on reports written before the routing fix; `no_reference` columns assigned by hand; any overlap counted as a hit. *S. pombe* mat1-M called "structurally unfindable" and the two Saccharomycetaceae `a` cassettes called "genuine misses" -- both are bar losses. | V | FIXED `3de6022`, `40cd34f` (correction appended) |

---

## Log

* **2026-09-25** Report started. Entries A1-A6, B1-B3, C1-C4, D1-D8, E1-E2
  compiled from the 2026-09-24/25 sessions, the two literature rounds
  (`docs/notes/2026-09-24_mat-reference-gap-literature*.md`) and project
  memory.
