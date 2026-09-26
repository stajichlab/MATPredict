# MAT reference gap lineages: literature and accession survey

Date: 2026-09-24. Scope: primary-source MAT locus characterizations for lineages
with no curated MATPredict record. Nothing in `db/` was changed.

## Method and what "verified" means

- Literature: PubMed (E-utilities and the PubMed tool), plus one web search.
- Accessions: every nucleotide accession below marked **Y** was fetched from
  NCBI nuccore (efetch, GenBank flat file) on 2026-09-24. The length, definition
  line, reference, and CDS list quoted come from that fetch.
- **G** = genome-derived. No locus-specific deposit exists. The coordinates were
  found by me with local `tblastn` (BLAST+ 2.17.0) of curated-type proteins
  against the named NCBI sequence. They are a measurement, not a published
  coordinate. Re-check before curation.
- **N** = no verified accession. Paper exists but deposits are genome-only or I
  did not locate them.
- NCBI search backends returned HTTP 500/429 errors late in the session, so a
  few protein-database searches (noted) did not complete.
- A title-field search misses some records. Example: the *Magnaporthe* MAT
  deposits did not match a `MAT1[title]` query and were found only via the
  PubMed-to-nuccore link. "No deposit found" therefore means "none found by
  title search plus PubMed link", not proof of absence.

## Summary table

| lineage | BFD genomes | best candidate record (species, idiomorph, accession) | verified? | notes |
|---|---|---|---|---|
| Morchellaceae | 60 | *Morchella importuna* YPL6-1 MAT1-2 `KY782629.1` (75.7 kb); YPL6-3 MAT1-1 `KY782630.1` (63.6 kb) | Y | Both idiomorphs with flanks: APN2 (2 CDS), SLA2 (end4), SDH2, MBA1, COX6A. Chai et al. 2017 Mycol Prog 16:743 (no PMID) |
| Discinaceae | 72 | none deposited; Dirks et al. 2025 (PMID 39788220) genome-derived | N | Mostly heterothallic; *G. esculenta* CBS101906 may carry both idiomorphs (needs confirmation) |
| Pyronemataceae | 9 | *Pyronema confluens* CBS100304: MAT1-1-1 `HF935381.1`, MAT1-2-1 `HF935431.1` | Y (scaffolds) | Homothallic. The two genes are on different scaffolds |
| Pezizaceae | 8 | Terfezia claveryi 2026 (PMID 42095936) | N | Accessions not checked |
| Rhizinaceae | 5 | *Phymatotrichopsis omnivora* MAT1-1 `MW357818.1` / `MW357819.1`; MAT1-2-1 `MW357816.1` | Y | Two MAT1-1-1 paralogs in tandem. No flanks |
| Tarzettaceae / Ascobolaceae | 3 / 1 | none | N | No published locus found |
| Orbiliaceae | 67 | none deposited | N | *A. oligospora* reported heterothallic (PMID 34576814); no locus deposit |
| Pleosporales | 1,012 | *Bipolaris maydis* C5 MAT-1 `AF029913.1`; C4 MAT-2 `AF027687.1`; *L. maculans* `AY174048.1`/`AY174049.1` | Y | Flanks GAP1, BGL1, ORF1, APN2 fragment. No SLA2. Didymellaceae homothallic fusions (PMID 22014305) |
| Botryosphaeriales | 807 | *Diplodia sapinea* CMW190 MAT1-1 `KF551229.1` (30.7 kb); CMW39103 MAT1-2 `KF551228.1` (30.0 kb) | Y | APN2, COX6A, APC5 flanks. MAT1-1-4 and MAT1-2-5 present. *Phyllosticta* `KT708823.1`/`KT708824.1` |
| Mycosphaerellales | 411 | *Z. tritici* MAT1-2 `AF440398.1`, MAT1-1 `AF440399.1`; *C. beticola* `KC960688.1`/`KC960689.1` | Y | APC5, APN2 flanks. *Cercospora* carries MAT exon fragments elsewhere in the genome |
| Dothideales | 204 | *Aureobasidium zeae* JIMHK-9 `MK416192.2` (both genes) | Y (unpublished) | Homothallic (Gostinčar 2014) |
| Venturiales | 112 | *V. inaequalis* MNH120 MAT1-1 `MG818328.1`; 1639 MAT1-2 `MG818329.1` | Y | APN2 flank |
| Cladosporiales | 91 | *Fulvia fulva* MAT1-1 `DQ659350.2`, MAT1-2 `DQ659351.2` | Y | NCBI puts *Fulvia* in Mycosphaerellaceae. No *Cladosporium* s.s. locus deposit found |
| Myriangiales | 28 | *Elsinoe australis* MAT1-1-1 `MN607981.1`, MAT1-2-1 `MN607982.1` | Y (unpublished) | Short (3.9/4.2 kb), no flanks |
| Glomerellales | 573 | *V. dahliae* MAT1-2 `AB505214.1`, MAT1-1 `AB505215.1`; *Glomerella* MAT1-2 `AY357890.1` | Y | *Colletotrichum*: atypical, MAT1-2 in both partners (PMID 22223174) |
| Magnaporthales | 544 | *M. oryzae* 70-6 MAT1-1 `AB080670.2`; 70-14 MAT1-2 `AB080671.2` | Y | Four MAT1-1 genes, three MAT1-2 genes. No flanks |
| Xylariales | 257 | none | N | Robinson & Natvig 2019: no canonical MAT regions in the order |
| Microascales | 164 | *Ceratocystis fimbriata* CBS114723 `KF033902.1` (unswitched) / `KF033903.1` (switched) | Y | Unidirectional switching by deletion. SLA2, APN2, APC5, COX6 flanks |
| Amphisphaeriales | 98 | none | N | Only genome-predicted mRNAs |
| Chaetothyriales | 223 | *Exophiala dermatitidis* CBS 132752 MAT1-1 `MH341450.1` | Y | MAT1-1-4 present |
| Thelebolales | 40 | *Pseudogymnoascus destructans* MAT1-2 `KJ938434.1`, MAT1-1 `KJ938437.1`; *P. roseus* homothallic `KJ938436.1` | Y | `KJ938436.1` spans APN2 to SLA2 |
| Umbilicariales / Caliciales | 47 / small | none | N | No locus deposit found |
| Peltigerales | small | *Lobaria pulmonaria* `JX520967.1` / `JX520966.1` | Y, but UNVERIFIED flag in GenBank; no CDS | Not usable as protein reference |
| Serinales (CTG) | 2,368 | *C. lusitaniae* CL143 MTLa `FJ524850.1`; *M. guilliermondii* MTLalpha `FJ524851.1`; *C. parapsilosis* MTLa `AY961981.1` | Y | a1 pseudogene in *C. parapsilosis*. *C. auris*: genome loci (G) |
| Pichiales | 579 | *K. phaffii* CBS 7435 `FR839631.1` (a1 4058..4735; alpha1 c139129..139701); *O. polymorpha* `AJ617305.1` | Y | Flip/flop inversion switching (138 kb / 19 kb) |
| Phaffomycetales, Ascoideales, Saccharomycodales, Dipodascales (non-Yarrowia), Lipomycetales, Trigonopsidales | 188 / 80 / 71 / 398 / 59 / 16 | no locus deposits found | N | Krassowski et al. 2019 (PMID 31353182) inferred MAT systems from 332 genomes. *Ascoidea rubescens* a2 `XP_020044210.1` (protein only) |
| Taphrinales (beyond *T. deformans*) | 22 | none | N | *Protomyces* genomes exist (PMID 33741074). No MAT analysis found |
| *Saitoella* / *Neolecta* | 2 / 1 | none | N | No published MAT locus found |

## Tier 1: Pezizomycetes and Orbiliomycetes

### Morchellaceae (60 genomes)

Papers:
- Chai H et al. 2017. "Characterization of mating-type idiomorphs suggests that
  *Morchella importuna*, Mel-20 and *M. sextelata* are heterothallic." Mycol Prog
  16(7):743-752. Not in PubMed. Cited in the GenBank records.
- Chai H et al. 2019. Mycologia 111:551-562. PMID 31251705,
  doi:10.1080/00275514.2019.1628553. Mel-20 idiomorphs: MAT1-1 7.8 kb (MAT1-1-1
  plus a second MAT1-1 gene), MAT1-2 7.5 kb.
- Chai H et al. 2022. J Fungi 8:746. PMID 35887501,
  doi:10.3390/jof8070746. "Unconventional integration" of MAT loci.
- Du XH & Yang ZL 2021 review. MMBR. PMID 34319143, doi:10.1128/MMBR.00220-20.

Verified accessions (all Chai et al. 2017):
- `KY782629.1`, 75,706 bp. "*Morchella importuna* strain YPL6-1 mating type locus
  MAT1-2-1 and flanking genes genomic sequence." 25 CDS. They include MAT1-2-1
  (AVI60809.1), two APN2 CDS (AVI60802/3), SLA2/end4 (AVI60811.1), SDH2, MBA1,
  COX6A, APC-related.
- `KY782630.1`, 63,555 bp. YPL6-3, MAT1-1-1 (AVI60816.1). Same flank set.
- `KY782631.1`, 24,180 bp. *M. semilibera* M115 MAT1-1-1.
- `KY782632.1`, 23,538 bp. *M. semilibera* M115 MAT1-2-1, with end4 (SLA2).
- Single-gene records: MAT1-1-10 `MG681026.1` and MAT1-1-11 `MG681027.1`
  (*M. importuna* YPL6). These show extra MAT1-1 genes in *Morchella*.

Architecture: heterothallic. The flank neighbourhood includes SDH2 and MBA1 next
to MAT, not only SLA2/APN2. **Best record: KY782629.1 + KY782630.1.** Note the
coordinator wrote "Du et al."; the locus deposits are from Chai et al.

### Discinaceae (72 genomes; *Gyromitra*)

- Dirks AC et al. 2025. Mol Phylogenet Evol 205:108286. PMID 39788220,
  doi:10.1016/j.ympev.2025.108286. 75 draft genomes. Mating-type loci identified
  from genomes. "Predominantly heterothallic." One colocalized MAT1-1/MAT1-2 case
  in *G. esculenta* CBS101906, which the authors say needs confirmation.
- No locus-specific nuccore deposit found (title and all-field searches, 0 hits).
- The genome-to-locus coordinates were not extracted. **Unverified.**

### Pyronemataceae (9 genomes)

- Traeger S et al. 2013. PLoS Genet 9:e1003820. PMID 24068976,
  doi:10.1371/journal.pgen.1003820. *Pyronema confluens* CBS100304 genome.
- Proteins in that assembly:
  - CCX07902.1, 340 aa, "Similar to mating type protein MAT1-1-1," coded on
    `HF935381.1` (33,197 bp scaffold, verified) 5918..7124, 3 exons.
  - CCX30263.1, 274 aa, "Similar to Mating-type protein a-1" (HMG, MAT1-2-1
    type), on `HF935431.1` (116,018 bp scaffold, verified) 101354..102335, 4 exons.
- Architecture: homothallic. Both genes are present, on separate scaffolds. I
  did not check whether the scaffolds join.
- *Ascodesmis nigricans*: genome in Lütkenhaus et al. 2019 Genetics (PMID
  31604798). No MAT analysis found.

### Pezizaceae (8 genomes)

- Andreu-Ardil L et al. 2026. "*Terfezia claveryi* MAT locus characterization…"
  Mycorrhiza. PMID 42095936, doi:10.1007/s00572-026-01266-3. I did not read it
  or check its accessions. **Unverified.** Worth reading first for this family.

### Rhizinaceae (5 genomes)

- Mattupalli C et al. 2021. Phytopathology. PMID 33728936 (genome resources);
  Mattupalli et al. 2022 Plant Dis, PMID 35156845 (heterothallic mating system).
- `MW357818.1`, 8,233 bp, *P. omnivora* NFPo20. "MAT1-1-1 gene, MAT1-1-1-1
  allele, complete cds; and MAT1-1-1 gene, MAT1-1-1-2 allele, partial." Two
  MAT1-1-1 CDS in opposite orientation.
- `MW357819.1`, 8,611 bp, NFPo30. Same structure.
- `MW357816.1` (1,057 bp) / `MW357817.1` (1,083 bp): MAT1-2-1 genes only.
- No SLA2/APN2 in these deposits. Two MAT1-1-1 copies are real biology per the
  deposit. Do not collapse them.

### Tarzettaceae (3), Ascobolaceae (1)

- No MAT locus paper or deposit found. The Tarzettaceae hit was a *Geopyxis
  carbonaria* WGS scaffold, which is not a MAT deposit.

### Orbiliaceae (67 genomes)

- Zhou D et al. 2021. Microorganisms 9:1919. PMID 34576814,
  doi:10.3390/microorganisms9091919. *A. oligospora* mating types typed by PCR
  in 239 isolates. Reported as heterothallic (84 MAT1-1 : 113 MAT1-2 in a
  subset). This is PCR typing, not a locus characterization.
- No nuccore locus deposit found. Protein search for *Arthrobotrys* MAT proteins
  by title returned 0.
- **No published characterized locus found.** A reference would have to be built
  from a genome.

## Tier 2: Dothideomycetes

### Curator-supplied sources (priority)

**1. Woudenberg JHC, de Gruyter J, Crous PW, Zwiers LH. 2012.** Mol Plant
Pathol 13(4):350-362. PMID 22014305, doi:10.1111/j.1364-3703.2011.00751.x.
Didymellaceae (**Pleosporales**). Species: *Phoma clematidina*, *Didymella
vitalbina*, *D. clematidis*, *Peyronellaea pinodes*, *P. pinodella*, *Phoma
herbarum*.

Architecture:
- *D. clematidis*: homothallic. It arose from a single crossover between MAT1-1
  and MAT1-2.
- *Pey. pinodes*: homothallic. It arose from a crossover plus an inversion of the
  fused MAT1/2 locus.
- The others are heterothallic.
- Flanks: only a DNA lyase (APN2) fragment. UniProt lists 16-205 aa "DNA
  lyase-like" fragments.

Verified nuccore records (UniProt protein IDs from lit_pubmed:22014305):

| accession | bp | taxon / strain | content | UniProt MAT proteins |
|---|---|---|---|---|
| `JF815534.1` | 5,402 | *D. clematidis* CBS 123705 | fused MAT1-2-1 + MAT1-1-1, DNA lyase | G8EEZ7 (MAT1-2-1, 345 aa), G8EEZ8 (MAT1-1-1, 363 aa) |
| `JF815533.1` | 4,862 | *Pey. pinodes* CBS 235.55 | DNA lyase, MAT1-2-1, MAT1-1-1 | G8EEZ4 (338 aa), G8EEZ5 (357 aa) |
| `JF815526.1` | 4,329 | *P. herbarum* CBS 615.75 | MAT1-2 idiomorph | G8EEX4 (351 aa) |
| `JF815530.1` | 3,661 | *P. clematidina* CBS 102.66 | MAT1-2 | G8EEY6 (352 aa) |
| `JF815528.1` | 3,605 | *P. clematidina* CBS 196.64 | MAT1-1 | G8EEY0 (364 aa) |
| `JF815532.1` | 3,112 | *D. vitalbina* CBS 123706 | MAT1-2 | G8EEZ2 (349 aa) |
| `JF815527.1` | 3,005 | *D. vitalbina* CBS 123707 | MAT1-1 | G8EEX7 (363 aa) |
| `JF815531.1` | 3,090 | *Pey. pinodella* CBS 108.31 | MAT1-2 | G8EEY9 (350 aa) |
| `JF815529.1` | 2,680 | *Pey. pinodella* CBS 110.32 | MAT1-1 | G8EEY3 (357 aa) |

`JF815533.1` CDS are AER26938.1 (lyase), AER26939.1 (MAT1-2-1), and AER26940.1
(MAT1-1-1). `JF815534.1` CDS are AER26941-AER26944.

**2. Bolton MD et al. 2014.** Fungal Genet Biol 62:43-54. PMID 24216224,
doi:10.1016/j.fgb.2013.10.011. *Cercospora beticola* (**Mycosphaerellales**).

Architecture:
- Heterothallic at MAT1.
- Every isolate also carries exon fragments of BOTH MAT1-1-1 and MAT1-2-1 at
  other loci across the genome. These are homogenized by concerted evolution.
  The same pattern occurs in *C. zeae-maydis*, but not in *Z. tritici*.
- **Detection risk**: these fragments can give off-locus MAT hits and a false
  "both idiomorphs present" call in *Cercospora* genomes.

Flanks in the deposits: zinc-type alcohol dehydrogenase, a putative integral
membrane protein, and hypothetical proteins. No SLA2/APN2 is annotated.

Verified records:
- `KC960688.1`, 12,840 bp, MAT1-1-1 = AHX24198.1 (UniProt X2F8Y1-type).
- `KC960689.1`, 13,103 bp, MAT1-2-1 = AHX24203.1 (UniProt X2FEQ0).
- Haplotype single-gene records: `KC960675-79.1` (MAT1-1-1, 2,859 bp) and
  `KC960680-87.1` (MAT1-2-1, about 3,092 bp).
- Fragment loci ("ROI A-F"): `KF225566.1`-`KF225571.1`, 1,102-3,564 bp. Do NOT
  curate these as MAT loci. They are the retroposed fragments.

**3. University of Pretoria repository PDF.** This is the author manuscript of
Bihon W, Wingfield MJ, Slippers B, Duong TA, Wingfield BD. 2014. "MAT gene
idiomorphs suggest a heterothallic sexual cycle in a predominantly asexual and
important pine pathogen." Fungal Genet Biol 62:55-61. PMID 24220137,
doi:10.1016/j.fgb.2013.10.013. *Diplodia sapinea* (= *D. pinea*,
**Botryosphaeriales**). Genome project `AXCF00000000`.

Architecture:
- Heterothallic.
- The MAT1-1 idiomorph carries MAT1-1-1 (444 aa, 1 intron) and MAT1-1-4.
- The MAT1-2 idiomorph carries MAT1-2-1 (389-398 aa, 1 intron) and a novel
  MAT1-2-5.
- Each idiomorph contains a partial copy of the other idiomorph's accessory gene.
- APN2 (DNA lyase) sits on the 5' side, about 13 kb from MAT. Between them are
  COX6A (COX13), APC5, CIA30, a mitochondrial carrier, and SAICAR synthetase.
  The paper compares this to *Grosmannia clavigera*.

Verified records:
- `KF551229.1`, 30,723 bp, strain CMW 190, MAT1-1. CDS: DNA lyase AHA91685.1,
  Cox-VIa AHA91686.1, APC5 AHA91687.1, CIA30, mito carrier, SAICAR, MAT1-1-4
  AHA91691.1, MAT1-1-1 AHA91692.1, integral membrane, hypothetical.
- `KF551228.1`, 29,983 bp, strain CMW 39103, MAT1-2. Same flanks. MAT1-2-1
  AHA91681.1, MAT1-2-5 AHA91682.1.
- UniProt: V5NSF4 (MAT-1, 444 aa), V5NSV2 (MAT1-2-1, 389 aa), V5NSW1 (MAT1-1-4),
  V5NSE7 (MAT1-2-5), V5NSF0 (APC5), V5NU92 (APN2), V5NSV7 (COX VIa).

### Order mapping and citation follow-up

**Pleosporales (1,012).** Classic, well-anchored deposits. All are Y:
- Turgeon BG et al. 1993. Mol Gen Genet 238:270-284. PMID 8479433. *Cochliobolus
  heterostrophus* (= *Bipolaris maydis*).
  - `AF029913.1`, 12,572 bp, strain C5. GAP1 (c), ORF, MAT-1 (AAB82945.1), BGL1.
  - `AF027687.1`, 12,449 bp, strain C4. GAP1, ORF, MAT-2 (AAB84004.1), BGL1.
- Cozijnsen AJ & Howlett BJ 2003. Curr Genet 43:351-357. PMID 12679880.
  *Leptosphaeria maculans* (*Plenodomus lingam*).
  - `AY174048.1`, 9,823 bp. GAP1, ORF, MAT1-1 (AAO37757.1), DNA lyase (partial).
  - `AY174049.1`, 10,124 bp. MAT1-2 (AAO37761.1).
- Bennett RS et al. 2003. Fungal Genet Biol 40:25-37. PMID 12948511.
  *Parastagonospora nodorum*. `AY212018.1` (6,193 bp, MAT1) and `AY212019.1`
  (6,401 bp, MAT2).
- Rau D et al. 2005. Genome 48:855-869. PMID 16391692. *Pyrenophora teres*.
  `AY950585.1` (2,673 bp, MAT-1) and `AY950586.1` (2,830 bp, MAT-2). *P.
  tritici-repentis* `AM884596-619` (about 4.4 kb) are titled "mat1-1 gene and
  partial mat1-2 gene", which is unusual. Check before use.
- *Alternaria alternata* `AB009451.1` (2,667 bp, MAT1). Arie et al. database-only
  deposit. Short, no flanks.
- *Stemphylium* `AY335164`-series and `AY340940.1` (Inderbitzin et al. 2005).
  Several homothallic fused MAT1-1/MAT1-2 deposits.

Pleosporales flanks are GAP1, BGL1, and ORF1. SLA2 is not in these deposits. If
detect requires SLA2 for Pleosporales, the requirement will fail.

**Botryosphaeriales (807).**
- *D. sapinea* above.
- Nagel JH et al. 2018. Fungal Genet Biol 114:24-33. PMID 29530630,
  doi:10.1016/j.fgb.2018.03.003. *Botryosphaeria dothidea*, *Macrophomina*,
  *Diplodia*, *Lasiodiplodia*, *Neofusicoccum*. Loci were taken from genome
  assemblies. No nuccore links. 9 species homothallic, 7 heterothallic. MAT gene
  fragments flank the MAT regions, which is the same detection risk as
  *Cercospora*.
- Lopes A et al. 2017 (*Neofusicoccum*, PMID 28317541) and 2018 (*Diplodia*,
  PMID 29880198): single-gene deposits, e.g. `KX505932.1` MAT1-1-1.
- *Phyllosticta* (Phyllostictaceae): `KT708823.1` (13,405 bp, *P. citricarpa*
  Gc12, APN2 + MAT1-2-1) and `KT708824.1` (16,940 bp, *P. capitalensis* Gm33,
  APN2 + MAT1-1-1 + MAT1-2-1, homothallic). Both unpublished in GenBank. Also
  `KX280782-83` (*P. citricarpa*, APN2 + MAT1-1-4 + MAT1-1-1).

**Mycosphaerellales (411).**
- Waalwijk C et al. 2002. Fungal Genet Biol 35:277-286. PMID 11929216.
  *Zymoseptoria tritici*. `AF440398.1` (14,878 bp: pal1-like, APC, DNA lyase,
  mat1-2 AAL30836.1) and `AF440399.1` (7,531 bp: DNA lyase partial, mat1-1
  AAL30838.1).
- Conde-Ferraez L et al. 2007. Mol Plant Pathol 8:111-120. PMID 20507483.
  *Pseudocercospora fijiensis*. `DQ787016.1` (9,799 bp: APC, DNA lyase,
  MAT1-2-1) and `DQ787015.1` (5,243 bp, MAT1-1-1). *P. musicola* `GU057991.1` /
  `GU057992.1` and *P. eumusae* `GU046393.2` / `GU046394.2` (11-12.5 kb) also
  exist. Their paper was not checked.
- *Teratosphaeria*: `MN119556.1` (16,186 bp, *T. zuluensis* MAT1-1: APC5, APN2,
  MAT1-1-10, MAT1-1-1, COX6a) and `MN119559.1` (16,001 bp, *T. gauchensis*
  MAT1-2). Both unpublished in GenBank.
- *Cercospora beticola* above.

**Dothideales (204).**
- Gostinčar C et al. 2014. BMC Genomics 15:549. PMID 24984952,
  doi:10.1186/1471-2164-15-549. All four *Aureobasidium* genomes have a
  homothallic MAT locus. No locus deposit from this paper.
- `MK416192.2`, 7,831 bp, *A. zeae* JIMHK-9. APN2 (+), MAT1-2-1 (c), MAT1-1-1
  (+), PH-domain protein. Direct submission only, no paper.

**Venturiales (112).**
- Paper in the record: "Evidence for sexual reproduction: identification,
  frequency and spatial distribution of *Venturia effusa*…", Phytopathology 2018.
  PMID 29381450. I did not check the author list.
  *V. inaequalis* `MG818328.1` (12,722 bp, MNH120, MAT1-1-1 + APN2) and
  `MG818329.1` (13,133 bp, isolate 1639, MAT1-2-1 + APN2). *V. effusa*
  `MF167363.1` / `MF167364.1` come from the same work. *V. carpophila*
  `MN562204-05`.

**Cladosporiales (91).**
- Stergiopoulos I et al. 2007. Fungal Genet Biol 44:415-429. PMID 17178244.
  *Cladosporium fulvum* (*Fulvia fulva*). `DQ659350.2` (5,433 bp, MAT1-1) and
  `DQ659351.2` (6,344 bp, MAT1-2). NCBI taxonomy places *Fulvia* in
  Mycosphaerellaceae, not Cladosporiaceae.
- For *Cladosporium* sensu stricto, no locus deposit was found. That part of the
  order has no primary locus reference.

**Myriangiales (28).**
- *Elsinoe australis* `MN607981.1` (3,870 bp, MAT1-1-1 + hypothetical) and
  `MN607982.1` (4,185 bp, MAT1-2-1 + hypothetical). Direct submission, Nanjing
  Forestry University, no paper found. No flanking genes.

## Tier 3: Pezizomycotina orders without a record

**Glomerellales (573).**
- *Verticillium dahliae* (Plectosphaerellaceae): Usami T et al. 2009. J Gen
  Plant Pathol 75:422-427. No PMID.
  - `AB505214.1`, 8,496 bp, TV103, MAT1-2 idiomorph + flank.
  - `AB505215.1`, 5,849 bp, MAT1-1.
  - Also Short DP et al. 2014, PMID 25383550.
- *Colletotrichum*: `AY357890.1`, 11,592 bp, *Glomerella cingulata* MAT1-2.
  Unpublished GenBank title "…from the homothallic ascomycete *Glomerella
  cingulata*."
- **Architecture flag**: Menat J et al. 2012, Mycologia 104:641-649, PMID
  22223174, doi:10.3852/10-265. The MAT1-2 HMG box is present in BOTH partners of
  fertile crosses in *Glomerella*. MAT1-1-1 is not the mating determinant. An
  idiomorph call based on MAT1-1-1 vs MAT1-2-1 will be wrong in *Colletotrichum*.
  Many `LC936xxx` records are only apn2-MAT1-2 IGS markers.

**Magnaporthales (544).**
- Kang S et al. 1994. Genetics 138:289-296. PMID 7828813. No linked deposits.
- Kanamori M et al. 2007. Gene 403:6-17. PMID 17881155,
  doi:10.1016/j.gene.2007.06.015.
  - `AB080670.2`, 4,666 bp, isolate 70-6. MAT1-1-1, MAT1-1-2, MAT1-1-3a,
    MAT1-1-3b.
  - `AB080671.2`, 3,736 bp, isolate 70-14. MAT1-2-1, MAT1-2-2a, MAT1-2-2b.
  - Extra isolates: `AB080668-73`.
- These deposits carry no flanks. SLA2 is annotated in 70-15 (`XM_003720710.1`,
  MGG_02949).

**Xylariales (257).**
- Robinson AJ & Natvig DO 2019. Fungal Genet Biol 122:47-52. PMID 30557613,
  doi:10.1016/j.fgb.2018.12.004. 35 genomes, 15 genera. No MAT1-1-1 or MAT1-1-2
  candidates were found in any member. The MAT1-2-1/MAT1-1-3-like HMG genes found
  are highly divergent or are non-MAT HMG paralogs.
- **Architecture flag**: expect no call in this order. That is biology, not a
  miss.

**Microascales (164).**
- Wilken PM et al. 2014. PLoS ONE 9:e92180. PMID 24651494. *Ceratocystis
  fimbriata* CBS114723.
  - `KF033902.1`, 46,416 bp, unswitched. COX6, APN2, APC5, SLA2, MAT1-1-1,
    MAT1-2-7, MAT1-2-1, MAT1-1-2, importin-beta.
  - `KF033903.1`, 42,835 bp, switched. MAT1-2 genes deleted.
  - Unidirectional switching by DNA loss gives self-sterile progeny. A genome can
    legitimately carry either state.
- Aylward J et al. 2016. Fungal Genet Biol 96:47-57. PMID 27720822.
  *Knoxdaviesia proteae* `KX832966.1` (29,542 bp: SLA2, MAT1-2-7, MAT1-2-1, APN2,
  APC5, COX6a, Coq4). `KX832965.1` is flagged UNVERIFIED in GenBank.
- Wilson AM et al. 2018. Fungal Genet Biol 113:32-41. PMID 29409964.
  *Thielaviopsis punctulata* `KX989056.1` (16,967 bp). MAT1-1-2 sits in the
  MAT1-2 idiomorph. Also TPA records `BK010318-21`.

**Amphisphaeriales (98).**
- No locus paper found. Only genome-predicted *Apiospora* mRNAs, e.g.
  `XM_066825547.1` "mating type 1-2", partial.

**Chaetothyriales (223).**
- Metin B et al. 2019. Fungal Genet Biol. PMID 30611834,
  doi:10.1016/j.fgb.2018.12.011. *Exophiala dermatitidis* CBS 132752
  `MH341450.1` (12,011 bp: MAT1-1-4, MAT1-1-1, hypotheticals). GenBank lists the
  reference as "Unpublished". The PubMed paper has the same title.
- Teixeira MM et al. 2017. Stud Mycol. PMID 28348446. Genome survey of the
  order; loci not deposited separately.

**Thelebolales (40, Pseudogymnoascaceae in NCBI).**
- Palmer JM et al. 2014. G3. PMID 25053709. *Pseudogymnoascus destructans*.
  - `KJ938434.1`, 7,551 bp, MAT1-2: MAT1-2-1, MAT1-2-5.
  - `KJ938437.1`, 7,703 bp, MAT1-1: MAT1-1-3, MAT1-1-6, MAT1-1-1.
  - *P. roseus* `KJ938436.1`, 13,208 bp. Homothallic, APN2 through SLA2.

**Umbilicariales (47), Caliciales (small).**
- No MAT locus paper or deposit found. Only WGS/MAG scaffolds matched.

**Peltigerales (small).**
- Singh G et al. 2012. PLoS ONE 7:e51402. PMID 23236495. *Lobaria pulmonaria*
  `JX520967.1` (2,664 bp) / `JX520966.1` (1,360 bp). Both are flagged UNVERIFIED
  in GenBank and have no CDS. They are not usable as protein references.

## Tier 4: Saccharomycotina outside Saccharomycetaceae

The existing queue (C. albicans MTLa/MTLalpha, K. lactis, L. thermotolerans,
Yarrowia MATB) is in `matpredict_saccharomycotina_curation_queue.md` and is not
repeated here.

**Serinales (2,368).**
- Reedy JL, Floyd AM, Heitman J 2009. Curr Biol 19:891-899. PMID 19446455,
  doi:10.1016/j.cub.2009.04.058.
  - `FJ524850.1`, 14,633 bp, *Clavispora lusitaniae* CL143 MTLa. PAP1, OBP1,
    PIK1, a2 (ACS29266.1), a1 (ACS29267.1, 164 aa, 1 intron), RCY1.
  - `FJ524851.1`, 14,353 bp, *Meyerozyma guilliermondii* NRRL Y-2075 MTLalpha.
    HIP1, OBP1, PIK1, alpha1 (ACS29273.1, 210 aa), PAP1, RCY1.
  - Neither idiomorph carries alpha2. This matches the recorded MTLalpha2 clade
    absence.
- Logue ME et al. 2005. Eukaryot Cell 4:1009-1017. PMID 15947193.
  `AY961981.1`, 15,704 bp, *C. parapsilosis* CLIB214 MTLa. PAPa, OBPa, PIKa,
  MTLa2 (AAY33181.1). **MTLa1 is a pseudogene** (no CDS).
- Sai S et al. 2011. Eukaryot Cell. PMID 21335529. *C. orthopsilosis*
  `HQ696681.1` (16,160 bp) and *C. metapsilosis* `HQ696678.1` (14,971 bp). MTLa
  with GAP1 flank.
- *C. dubliniensis* `AY622606.1` (8,593 bp, MTLalpha2, OBP, PIK, alpha1, PAP
  partial). Pujol et al. 2004, PMID 15302834.
- *Candida auris* (Muñoz JF et al. 2018, Nat Commun 9:5346, PMID 30559369,
  doi:10.1038/s41467-018-07779-6). Clades carry either MTLa or MTLalpha. No
  locus deposit exists. My tblastn measurement (G):
  - MTLa, B8441 (clade I), `NC_140807.1` (chr 3): PAP1 1,487,923-1,489,548 (78%
    id to *C. lusitaniae* PAP1), OBP1 1,489,775-1,491,049, PIK1
    1,491,537-1,494,350, a2 1,494,442-1,495,046 (E=2e-11). a1: only a weak hit
    at 1,495,483-1,495,620 (E=0.005, 46 aa). **a1 presence is unverified.**
  - MTLalpha, B11221 (clade III), `CP126635.2` (chr 3): OBP 1,357,062-1,358,228,
    PIK 1,358,433-1,361,207, alpha1 1,361,372-1,361,791 (E=9e-17), PAP
    1,361,839-1,363,506. No alpha2 hit.
- *C. tropicalis*: Porman 2011 (PMID 22158989) and Xie 2012 (PMID 22544905)
  report mating biology. No locus deposit found. WGS proteins exist, e.g.
  KAK6891836.1 "Mating-type-like protein A1" on `JAIZWC020000234.1` (17,973 bp);
  alpha1 KAK6886132.1. Reference genome MYA-3404 scaffold `GG692408.1`. Not
  curated. **Unverified as a locus.**

**Pichiales (579).**
- Hanson SJ, Byrne KP, Wolfe KH 2014. PNAS 111:E4851-8. PMID 25349420,
  doi:10.1073/pnas.1416014111.
- Maekawa H & Kaneko Y 2014. PLoS Genet. PMID 25412462.
- *Komagataella phaffii* CBS 7435, `FR839631.1` (chromosome 4, 1,820,458 bp,
  verified). a1 CCA40182.1 at 4058..4735 (near the telomere). alpha1 CCA40251.1
  at complement(139129..139701). This is the 138-kb invertible region: the
  orientation puts one of MATa or MATalpha beside the telomere.
- *Ogataea polymorpha*:
  - `AJ617305.1`, 6,239 bp, *O. angusta* CBS 4732. yol077C, dic1, matalpha2,
    matalpha1, mata1, sla2 partial. Butler et al. 2004, PMID 14745027.
  - In BY4329 (Maekawa 2014 genome) `DF933571.1`: DIC1 913,308-914,168 (100% id),
    alpha2 915,290-915,808 (94%), alpha1 915,817-916,401 (100%). This is a G
    measurement.
  - The a1 and SLA2 proteins from `AJ617305.1` gave no hit (E<1) in
    880-960 kb. **Unresolved.** The a1 side lies about 19 kb away on the other
    end of the invertible region per Hanson 2014, but I did not locate it.
- *O. minuta* `LC373261.1` (10,766 bp, MAT1 and MAT2 regions, both flanked by
  SLA2 copies). PMID 30064813.
- **Architecture flag**: flip/flop inversion. Both idiomorphs are present in
  every haploid genome, and one is silenced. Detect will always find both.

**Phaffomycetales (188), Ascoideales (80), Saccharomycodales (71), Dipodascales
(398), Lipomycetales (59), Trigonopsidales (16).**
- Krassowski T et al. 2019. Curr Biol 29:2555-2562. PMID 31353182,
  doi:10.1016/j.cub.2019.06.056. MAT systems inferred from 332 genomes.
  Switching arose at least 11 times, flip/flop in at least 10 groups. This is the
  source to consult for per-order locus structure. I did not extract per-species
  coordinates from it.
- Wolfe KH & Butler G 2022. MMBR 86:e0000721. PMID 35195440. Review.
- Specific findings:
  - *Ascoidea rubescens* DSM 1968: a2 protein `XP_020044210.1` (305 aa, genome
    annotation, Riley et al. 2016 PNAS, PMID 27535936).
  - *Geotrichum candidum* (Dipodascaceae): `HF558448.1` (846 bp, MATA, CLIB 918)
    and `HF558449.1` (1,020 bp, MATB, CBS 615.84). Single genes, unpublished.
  - *Hanseniaspora*: Steenwyk JL et al. 2019. PLoS Biol. PMID 31112549. Loss of
    cell-cycle and DNA-repair genes. I did not confirm from the text whether MAT
    genes are lost. **Unverified.**
  - *Wickerhamomyces*, *Lipomyces*, *Trigonopsis*, *Blastobotrys*: no MAT paper
    or deposit found. Protein search returned 0, or failed with NCBI 500 errors
    for some genera.

## Tier 5: Taphrinomycotina

- Existing record: *Taphrina deformans* 5011 (Almeida et al. 2015, PMID
  25587012; Cissé et al. 2013, PMID 23631913). Fused M+P locus, homothallic.
- *Protomyces*: Wang K et al. 2021. IMA Fungus 12:8. PMID 33741074. Genomes of
  all available species. The abstract does not mention MAT. No MAT deposit found.
- *Saitoella complicata*: genome papers only (PMID 21914972, 26021914). No MAT
  analysis found.
- *Neolecta irregularis*: Nguyen TA et al. 2017. Nat Commun. PMID 28176784.
  The WGS hit `LXFE01003687.1` contains only a swi5 homolog. It is not a MAT
  locus. No MAT analysis found.

## Detection-relevant architecture flags (collected)

1. **Off-locus MAT fragments.** *Cercospora* (Bolton 2014) and Botryosphaeriaceae
   (Nagel 2018) have them. They can create false "both idiomorphs" calls.
2. **Homothallic fusions.** Didymellaceae (inversion in *Pey. pinodes*),
   *Stemphylium*, *Aureobasidium*, *P. capitalensis*, *Pyronema*,
   *Pseudogymnoascus roseus*.
3. **Accessory idiomorph genes.** MAT1-1-4 and MAT1-2-5 (*Diplodia*), MAT1-1-10
   (*Teratosphaeria*, *Morchella*), MAT1-2-7 (Microascales), MAT1-1-6
   (*Pseudogymnoascus*). None of these has a roster slot today, to my knowledge.
   I did not check the roster.
4. **Unidirectional switching.** *Ceratocystis*: the same strain can be MAT1-1
   only after switching.
5. **Idiomorph not defined by MAT1-1-1/MAT1-2-1.** *Colletotrichum*. MAT loci
   are absent or divergent in Xylariales.
6. **Flip/flop inversion.** *Komagataella* and *Ogataea*: both idiomorphs are in
   every genome.
7. **Gene loss.** MTLalpha2 is absent in *Clavispora* and *C. auris*. MTLa1 is a
   pseudogene in *C. parapsilosis*.
8. **Flank sets differ by class.** Pleosporales use GAP1/BGL1 (no SLA2).
   Botryosphaeriales and Mycosphaerellales use APN2/APC5/COX13. Morchellaceae
   have SDH2/MBA1 beside MAT. A detection rule that requires SLA2 fails in these
   lineages.
