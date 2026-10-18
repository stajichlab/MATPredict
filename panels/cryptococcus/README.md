# Cryptococcus neoformans / deneoformans / gattii-complex idiomorph panels for `matpredict reads-type`

Seven panels, each a pair `<panel>_MATa.fasta` / `<panel>_MATalpha.fasta` with one record per idiomorph
(record id = `source accession:coordinates`). Use one panel per lineage:

```
matpredict reads-type --idiomorph MATa=panels/cryptococcus/neoformans_MATa.fasta \
                      --idiomorph MATalpha=panels/cryptococcus/neoformans_MATalpha.fasta \
                      --reads R1.fastq.gz R2.fastq.gz --out call.tsv
```

## Why one panel per lineage (evidence)

The C. neoformans, C. deneoformans and gattii-complex MAT loci are 84-89% identical (blastn, aligned parts, mean identity of
HSPs >= 200 bp), which leaves almost no exact 31-mers in common:

| pair of alleles | mean identity of aligned blocks |
|---|---|
| neoformans MATalpha (H99) vs deneoformans MATalpha (JEC21) | 88.8% |
| neoformans MATa (125.91) vs deneoformans MATa (JEC20) | 88.1% |
| neoformans vs gattii (WM276 / E566) | 83.4-85.0% |
| gattii MATalpha lineages (WM276 VGI, R265 VGII, NIH312 VGIII, IND107 VGIV, 7685027 decagattii) | 93-97% (VGIII and decagattii/VGIV closest) |
| gattii MATa (E566, B4546 VGIII, CBS 10090 VGII) | 94-95% |

Consequences, tested on 39 genome assemblies and on real reads (`validation_assemblies.tsv`, `validation.tsv`):

* A panel built on one species gives `low_depth` or `none` on the other species (shared depth 0.1-4, breadth <= 0.06). It never
  gave the wrong idiomorph: 0 wrong a/alpha calls in 151 typings of assemblies (39 in-lineage + 112 cross-lineage) and 42 typings of real reads (13 in-lineage + 29 cross-lineage).
* Inside the gattii complex a MATa call crosses lineages (breadth ~0.3, one panel gives `a` for VGII/VGI MATa strains), a MATalpha
  call mostly does not (breadth ~0.2-0.46, below the 0.5 x expected rule), so lineage panels are needed for MATalpha.
* A combined gattii-complex panel (3 alpha + 2 a records in one FASTA per idiomorph, tested only on assemblies; `gattiiX`, not shipped)
  fails: of the 7 gattii-complex strains typed 3 were called correctly (WM276, CA1280, CA1873), 2 gave `low_depth` (R265, CBS 10090) and
  2 `none` (IND107, 7685027); shared k-mers vanish when several lineages are pooled and breadth falls to 0.2-0.4. Not recommended; lineage-specific panels are.

So: **neoformans, deneoformans, gattii (VGI), deuterogattii (VGII), bacillisporus (VGIII), tetragattii (VGIV), decagattii**.
Pick the panel by species/lineage of the strain (the classifier does not need the answer to be exact: a wrong panel returns
`none`/`low_depth`, so running all panels and taking the one with a call is safe).

## Why the FASTA is the whole locus and not a trimmed "specific region"

`Panel.from_sequences` defines idiomorph-unique k-mers as those absent from the other idiomorph, and `shared` k-mers as those present in both;
`shared_depth` (the single-copy control and the `low_depth` gate) is computed from the shared k-mers. A FASTA trimmed to the specific region
would have no shared k-mers and every strain would be `low_depth`. The records therefore include the shared flanks; the panel code removes
the flank k-mers from the idiomorph-specific sets (k-mers made unique by SNPs in flanks are dropped by `--min-unique-run 100`; residual
unique k-mers inside the flanks, where a flank could be delimited (deneoformans, gattii, deuterogattii): 0 gattii, 0 deneoformans MATalpha,
287 deneoformans MATa, 547 / 341 deuterogattii MATa / MATalpha; neoformans and the HQ3961xx VGIII records have no delimitable 10 kb flank in the blast alignment).

The idiomorph-specific region was located by blastn of MATa vs MATalpha record of the same species (E-value 1e-10, HSPs >= 100 bp);
the specific part is the complement of the union of aligned blocks, 10 kb inside from the left end:

| species | MATa record: unaligned bp (of length) | MATalpha record: unaligned bp (of length) | left flank | right flank |
|---|---|---|---|---|
| neoformans | 88,522 (153,973) | 60,248 (124,734) | first 10,000 bp | ~10 kb |
| deneoformans | 67,364 (148,067) | 52,998 (133,474) | 1-10,000 (99.2%) | 130,253-148,067 / 115,658-133,474 (99.5%) |
| gattii | 77,804 (144,400) | 66,167 (132,746) | 1-10,000 | 134,401-144,400 / 122,747-132,741 (99.8%) |

Unaligned stretches are interleaved with diverged (66-98% identity) homologous genes (STE20, STE12, NOG2, RPL22 ... alleles), which are
also allele-specific as exact 31-mers. Hence the unique k-mer sets (below) cover far more than the unaligned bp.

## Records (provenance)

| file | record id | source |
|---|---|---|
| neoformans_MATalpha | AF542529.2:1-124734 | C. neoformans var. grubii H99 MAT alpha locus (GenBank, Fraser/Heitman lab) |
| neoformans_MATa | AF542528.2:1-153973 | C. neoformans var. grubii 125.91 MAT a locus (GenBank) |
| deneoformans_MATalpha | AF542531.2:1-133474 | JEC21 MAT alpha locus (GenBank; same record as db 40410_jec21_MAT_alpha) |
| deneoformans_MATa | AF542530.2:1-148067 | JEC20 MAT a locus (GenBank; db 40410_jec20_MAT_a) |
| gattii_MATalpha | AY710430.1:1-132746 | C. gattii WM276 MATalpha locus (GenBank) |
| gattii_MATa | AY710429.1:1-144400 | C. gattii E566 MATa locus (GenBank; E566 is serotype B, molecular type not verified) |
| deuterogattii_MATalpha | DQ096809.1:1-120440 | C. gattii R265 MAT alpha locus (GenBank) |
| deuterogattii_MATa | CM153009.1:131910-260575 | CBS 10090 (C. deuterogattii) MATa, cut from assembly GCA_056632935.1 (flanks located with blastn of the WM276 flanks; forward strand) |
| bacillisporus_MATalpha | HQ396150.1:1-116198 | C. gattii VGIII NIH312 mating type locus (GenBank; MATalpha, 99.8% to CA1280 locus) |
| bacillisporus_MATa | HQ396149.1:1-114300 | C. gattii VGIII B4546 mating type locus (GenBank; MATa, 99.9% to CA1873 locus) |
| tetragattii_MATalpha | NC_089950.1:130631-246237 | C. tetragattii IND107, GCF_000835755.1 |
| tetragattii_MATa | HQ396149.1:1-114300 | PROXY: VGIII B4546 MATa. No VGIV MATa was found in any public assembly. |
| decagattii_MATalpha | NC_088593.1:129784-256537 | C. decagattii 7685027, GCF_036417295.1 |
| decagattii_MATa | HQ396149.1:1-114300 | PROXY: VGIII B4546 MATa. No decagattii MATa is public. |

Allele (a vs alpha) of every assembly used was assigned independently of the k-mer panels by tblastn of the JEC21 SXI1 (alpha-only) and
JEC20 SXI2 (a-only) proteins (`db/.../40410_jec21_MAT_alpha/proteins.faa`, `..._jec20_MAT_a/proteins.faa`): the MATalpha strains hit SXI1
(bitscore 400-620, 65-87% id) and have no SXI2 (<= 40); MATa strains the opposite (SXI2 466-889).

Assemblies checked (accession: strain, allele): neoformans GCF_000149245.1 H99 alpha, GCA_002216725.1 KN99 (alpha), GCA_010233705.1 Ze90-1 alpha,
GCA_056632715.1 AD1-7a alpha, GCA_056632795.1 Bt63 a, GCA_056632815.1 Bt206 a, GCA_056632835.1 Bt1 alpha, GCA_056622635.1 T4 alpha,
GCA_043513845.1 Bt210 alpha, GCA_056632885.1 125.91 a; deneoformans GCF_000091045.1 JEC21 alpha, GCA_056621545.1 JEC20 a, GCA_002216205.1 XL280 alpha,
GCA_023523875.1 XL280alpha, GCA_056632855.1 NIH12 alpha, GCA_056632945.1 NIH433 a, GCF_000149385.1 B-3501A alpha; gattii complex GCF_000185945.1 and
GCA_036428525.1 WM276 alpha, GCF_002954075.1 R265 alpha, GCA_000499585.1 CBS 7750 alpha, GCA_000835815.1 2001/935-1 alpha, GCA_000836375.1 Ram5 alpha,
GCA_000836455.1 99/473 alpha, GCA_000875795.1 CA1014 alpha, GCA_000875855.1 MMRL2647 alpha, GCA_012922885.1 Chr10-7 alpha, GCA_000836315.1 LA55 a,
GCA_056632935.1 CBS 10090 a, GCA_000835745.1 EJB2 alpha, GCA_000836355.1 Ru294 alpha, GCA_000935105.1 NT-10 alpha, GCA_031841975.1 SR061 alpha,
GCA_000875815.1 E566 a, GCA_014964225.1 WM1802 a, GCF_000836335.1 CA1280 alpha, GCA_000855695.2 CA1873 a, GCF_000835755.1 IND107 alpha, GCF_036417295.1 7685027 alpha.
No MATa isolate was found for VGIV (tetragattii) or decagattii; the only public MATa gattii assemblies are E566, WM1802, LA55, CBS 10090, CA1873.
No KN99a assembly exists in the NCBI taxon 5206 genome list (KN99 GCA_002216725.1 is alpha).

## k-mer counts (k = 31, default `--min-unique-run 100`)

| panel | MATa unique | MATalpha unique | shared (single-copy control) |
|---|---|---|---|
| neoformans | 109,791 | 87,697 | 23,862 |
| deneoformans | 88,903 | 79,571 | 41,229 |
| gattii | 97,070 | 89,958 | 31,434 |
| deuterogattii | 89,028 | 84,219 | 25,946 |
| bacillisporus | 91,261 | 95,534 | 12,511 |
| tetragattii | 96,091 | 98,818 | 6,670 |
| decagattii | 92,644 | 106,931 | 10,264 |

bacillisporus, tetragattii and decagattii have a small shared set because the GenBank VGIII records (HQ396149/50) carry almost no flank
and the tetragattii/decagattii pairs mix lineages.

## Validation

`validation.tsv`: 13 real-read runs with independently known allele, typed with their own lineage panel (`in-lineage`) plus the other
panels (`cross-lineage`). Reads: first 2M pairs (8M FASTQ lines, ~4M reads) for runs whose FASTQ is not position-sorted
(contig share of the first 100k reads 0.10-0.13); for old runs (HiSeq 2000 era, SRR576xxx, SRR797/798xxx) the FASTQ is position sorted
(99% of the first 100k reads map to one H99 contig), so `head` is useless and the whole file was streamed and every Nth read kept
(1/10 to 1/12); two of those (CA1280, CA1873) had a CRC error in file 1, leaving 76% / 58% of file 1 (partial genome, MAT locus
was covered) and CBS 10090 file 2 has 96%. Reads were deleted after typing.

Result (own-lineage panel): **13/13 correct** (6 alpha, 7 a), no `low_depth`, no `both`/`none`.

| lineage | strain (run) | expected | call | breadth a / alpha | shared depth |
|---|---|---|---|---|---|
| neoformans | H99 (SRR30210006) | alpha | alpha | 0.02 / 1.00 | 21.5 |
| neoformans | Bt1 (SRR576318) | alpha | alpha | 0.03 / 0.76 | 5.4 |
| neoformans | 125.91 (SRR35626592) | a | a | 1.00 / 0.02 | 22.6 |
| neoformans | Bt63 (SRR576308) | a | a | 0.87 / 0.04 | 6.2 |
| deneoformans | JEC21 (SRR39693470) | alpha | alpha | 0.05 / 1.00 | 18.9 |
| deneoformans | XL280 (SRR1505177) | alpha | alpha | 0.04 / 0.99 | 5.5 |
| deneoformans | NIH12 (SRR15309084) | alpha | alpha | 0.05 / 1.00 | 24.3 |
| deneoformans | JEC20 (SRR11576540) | a | a | 1.00 / 0.07 | 19.9 |
| deneoformans | NIH433 (SRR15309083) | a | a | 1.00 / 0.08 | 28.5 |
| deuterogattii | R265 (SRR10500147) | alpha | alpha | 0.00 / 1.00 | 19.6 |
| deuterogattii | CBS 10090 (SRR798761) | a | a | 0.99 / 0.00 | 5.3 |
| bacillisporus | CA1280 (SRR797769) | alpha | alpha | 0.01 / 0.81 | 2.6 |
| bacillisporus | CA1873 (SRR798818) | a | a | 0.94 / 0.02 | 4.9 |

Caveats on independence: H99, 125.91, JEC20, JEC21, R265 are the strains the panel records came from; CBS 10090 supplied the deuterogattii MATa
record; CA1280/CA1873 loci are ~99.8% identical to the NIH312/B4546 records. Independent strains: Bt1, Bt63, XL280, NIH12, NIH433 (the last
two are the MAT parents of JEC21/JEC20), but they are the same lineage. 29 cross-lineage typings gave `none` (18) or `low_depth` (11), never a wrong allele.

`validation_assemblies.tsv`: all assemblies listed above converted to 100 bp reads tiled every 25 bp (4x, no errors; this tests allele
divergence only, not read errors or depth) and typed with every relevant panel. In-lineage: 37/39 correct. The two non-calls are
data-quality/labelling issues: SR061 (NCBI "C. gattii") is typed alpha by the deuterogattii panel (0.95), i.e. it is VGII and the label is
wrong; WM1802 (a) matches the VGIII a-allele better (breadth 0.62 with bacillisporus panel, 0.22 gattii panel). IND107 (VGIV) is correctly called
alpha by its own tetragattii panel (1.00) but not by the other panels (and by `gattii` at 0.43); 7685027 (decagattii) alpha is called by
decagattii (1.00) and by bacillisporus (0.46).

## Not verified / limitations

* No MATa allele exists for VGIV or decagattii in public data; tetragattii/decagattii panels use the VGIII B4546 MATa as proxy. A true MATa
  strain of these lineages is likely to return `none` or a weak `a` (a-alleles cross lineages at breadth ~0.3, < 0.5 x expected).
  They are validated for MATalpha on assemblies only (no real tetragattii/decagattii reads were typed).
* AD hybrids (neoformans x deneoformans, serotype AD) were not tested; `both` would require reads for both lineages' panels.
* KN99a, JEC20-derived cross progeny, and manifest rows with an uncertain label (e.g. SRR1520841 "R265" labelled VGIII bacillisporus; "H99" runs
  labelled deneoformans) were not used as truth. Cross-progeny runs ("Progeny of KN99a and C##") were not used because their MAT allele is not recorded.
* E566's molecular type (VGI vs VGII) and the lineage of WM1802 were not verified.
* Within-lineage diversity: neoformans MATalpha breadth was 0.76 (Bt1) and 0.67 (T4, assembly) versus 1.00 for H99: other VNI/VNII/VNB alleles
  work, but divergent strains are closer to the 0.5 threshold; neoformans MATa Bt206 breadth 0.75 (assembly).
