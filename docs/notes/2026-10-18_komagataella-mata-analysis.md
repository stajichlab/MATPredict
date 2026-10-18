# Why the Komagataella MATa genes have no good profile match, 2026-10-18

Branch `komagataella-curation`. Question from the curator: "analyze the MATa genes more: why can't we get a good match for these?" This note answers it for MATa1 and MATa2 of *K. phaffii* CBS 7435 and draws a recommendation (not a change) for `gene_class` and `idiomorph_class` of `MATkp` and `MATdc`. Scripts and raw outputs are in `2026-10-18_komagataella-mata-analysis_files/` (`an1.py` profiles and jackhmmer, `an2.py` iteration/alignment/ORF/tblastn, `an3.py` profile building and leave-one-lineage-out; `an*.out` outputs; `fs_summary.txt` Foldseek; `af_summary.txt` AlphaFold geometry; `idmat.txt` identities). They ran on the HPCC (jobs `kom-an1..3`, pyhmmer 0.12.3, MAFFT).

## Short answer

1. **MATa1 is fine**: it carries a homeodomain (Pfam PF00046 51.7 bits, GA 22.9), like MATa1 of every other yeast tested. The profile problem is **MATa2 only**.
2. **MATa2 carries a divergent HMG box.** Three independent lines, none of them conclusive alone: (a) iterative search (jackhmmer, 6 iterations) from Komagataella MATa2 pulls in the HMG-box proteins of every yeast tested (Ixr1, Abf2, Rox1, Hmo1, Rfg1, "HMG box domain-containing protein" of Lachancea, Kluyveromyces, Candida, Vanderwaltozyma, Ogataea, Pichia; the 5-6 iteration runs differ by a few bits) and the Pezizomycotina MAT1-2-1 proteins of the database at 100-128 bits; (b) the AlphaFold models of the Komagataella, *K. lactis* and *C. albicans* MATa2 each show one structured 75-90 residue region of three helices inside long disordered termini, and Foldseek against PDB100 returns HMG-box structures (HMGB1/HMGB2, SOX, LEF1) as the top hits for all three (E 0.02-0.6, so borderline individually); (c) the MAFFT alignment of the MATa2 core to Nhp6a, SRY, Ixr1 and Aspergillus MAT1-2-1 lines up the HMG-box start (RPxNxF/RPxNAF) and the helix-2 `S/K..S..W` motif.
3. **Why the existing profiles miss it**: the box is short and degenerate and the orthologues have diverged from each other. Core identity between MATa2 orthologues is only 18-43% (table below); the first phmmer iteration from Komagataella MATa2 finds only *C. albicans* MTLA2 (38 bits) among the other MATa2 proteins (not *K. lactis*, Geotrichum or Yarrowia). Pfam HMG_box scores Komagataella MATa2 at 13.6 bits (GA 25.4) and finds nothing (E<=1000) in *K. lactis*, *C. albicans*, *D. hansenii*, *C. lusitaniae*, Geotrichum or Yarrowia MATa2; HMG_box_2 and HMG_box_5 find nothing. The curator MAT1-2-1 profile (trained on Pezizomycotina) finds nothing except Geotrichum MATA (39.9 bits). It is not truncation and not a different fold (below).
4. **A Saccharomycotina MATa2 profile works**: built from 7 MATa2 proteins, leave-one-lineage-out recovers every held-out MATa2 at 22-52 bits (E 1e-3 to 1e-13) where Pfam and the curator profile give nothing, and the profile built on all 7 separates them from generic HMG proteins (139-262 bits versus 33 bits for the best HMG protein). The margin in leave-one-out is thin for the most divergent lineages (see the caveat on the table).

## 1. Profile scores on yeast MAT proteins

Best-domain bit score / E-value (hits with E<=10 only; "-" = no hit). Profiles: Pfam from the InterPro API (PF00046 Homeodomain GA 22.9; PF00505 HMG_box GA 25.4; PF04769 MATalpha_HMGbox GA 28.1; PF05920 Homeobox_KN GA 24.0; PF09011 HMG_box_2 GA 22.4; PF14887 HMG_box_5 GA 25.0) and the curator profiles `MAT_core/MAT1-1-1.hmm`, `MAT1-2-1.hmm` (no GA). Sequences: Kp = CBS 7435 (db records); Sc, Kl, Ca, Yl from UniProt (P0CY11 HMRA1, P0CY13 HMRA2, P0CY06, P0CY08, Q9HG12, Q9HG13, Q08398, Q9UW23, Q71U11, Q9UW19, Q9UW22, Q707Y0) and the curated db records (Dh, Cl, Yl mata1/mata2/matb1/matb2, Geotrichum MATA/MATB); controls Nhp6a (P11632), SRY (Q05066), Ixr1 (P33417) and Aspergillus fumigatus Af293 MAT1-2-1 (db).

| protein | Homeodomain PF00046 | HMG_box PF00505 | MATalpha_HMGbox PF04769 | Homeobox_KN PF05920 | HMG_box_2 PF09011 | HMG_box_5 PF14887 | curator MAT1-1-1 | curator MAT1-2-1 |
|---|---|---|---|---|---|---|---|---|
| Kp_MATa1 | 51.7 2e-16 | - | - | 24.2 9e-08 | - | - | - | - |
| **Kp_MATa2** | - | 13.6 3e-04 | - | - | - | - | - | - |
| Kp_MATalpha1 | - | - | 37.5 5e-12 | - | - | - | 30.5 1e-09 | - |
| Kp_MATalpha2 | 27.0 1e-08 | - | - | 43.7 7e-14 | - | - | - | - |
| Sc_HMRA1 | 50.9 4e-16 | - | - | 32.7 2e-10 | - | - | - | - |
| Sc_HMRA2 (alpha2-like remnant) | 38.7 3e-12 | - | - | 35.7 2e-11 | - | - | - | - |
| Sc_MATALPHA1 | - | - | 51.5 3e-16 | - | - | - | 33.5 2e-10 | - |
| Sc_MATALPHA2 | 36.4 1e-11 | - | - | 33.3 1e-10 | - | - | - | - |
| Kl_MATA1 | 56.3 8e-18 | - | - | 23.2 2e-07 | - | - | - | - |
| **Kl_MATA2** | - | - | - | - | - | - | - | - |
| Kl_ALPHA1 | - | - | 125.5 6e-39 | - | - | - | 42.4 3e-13 | - |
| Ca_MTLA1 | 51.8 2e-16 | - | - | 25.5 3e-08 | - | - | - | - |
| **Ca_MTLA2** | - | - | - | - | - | - | - | - |
| Ca_ALPHA1 | - | - | 42.3 2e-13 | - | - | - | 33.7 1e-10 | - |
| Ca_ALPHA2 | 32.9 2e-10 | - | - | 39.2 2e-12 | - | - | - | - |
| Dh_MTLA1 | 47.0 7e-15 | - | - | 36.0 2e-11 | - | - | - | - |
| **Dh_MTLA2** | - | - | - | - | - | - | - | - |
| Cl_MTLA1 | 48.9 2e-15 | - | - | 21.6 6e-07 | - | - | - | - |
| **Cl_MTLA2** | - | - | - | - | - | - | - | - |
| Yl_mata1 | 34.0 8e-11 | - | - | 24.8 5e-08 | - | - | - | - |
| **Yl_mata2** | - | - | - | - | - | - | - | - |
| Yl_matb1 (db) | - | - | 15.6 3e-05 | - | - | - | - | - |
| Yl_matb2 (db) | 12.6 4e-04 | - | - | - | - | - | - | - |
| **Geo_MATA** | - | - | 13.9 1e-04 | - | - | - | - | 39.9 1e-12 |
| Geo_MATB | - | - | 35.5 2e-11 | - | - | - | 22.6 3e-07 | - |
| Nhp6a (HMG control) | - | 91.4 1e-28 | - | - | 71.1 4e-22 | - | - | 33.5 1e-10 |
| SRY (HMG control) | - | 93.1 1e-29 | - | - | 35.7 4e-11 | - | - | 78.8 2e-24 |
| Ixr1 (HMG control) | - | 118.3 6e-37 | - | - | 57.9 5e-18 | - | - | 25.7 3e-08 |
| Afum MAT1-2-1 | - | 76.7 6e-24 | - | - | 31.3 9e-10 | - | - | 376.5 6e-115 |

Reading: every a1, alpha1 and alpha2 protein is seen by the profile of its domain (homeodomain, alpha box). Every MATa2 protein, of five lineages, is invisible to Pfam HMG_box and to the curator MAT1-2-1 profile (Geotrichum MATA excepted), while the HMG controls and a Pezizomycotina MAT1-2-1 hit strongly. So the Komagataella MATa2 result is the general MATa2 result, not a Komagataella artefact. InterPro agrees on the same proteins: the InterPro API lists a homeodomain for Kp MATa1 (F2QXQ4), Kl MATA1 (Q9HG12), Ca MTLA1 (Q9UW23), alpha box for Kp MATalpha1 (F2QXX3) and Ca MTLALPHA1, homeodomain/KN for Kp and Ca alpha2, and NO InterPro entry at all for Kp MATa2 (A0A1G4KQN4), Kl MATA2 (Q9HG13) or Ca MTLA2 (Q71U11) (the API returned no entries for these three). Note that InterPro lists a homeodomain for Yarrowia MATB2 (Q707Y0, `GN=MATB2`; CDD/SMART/PROSITE-type hits) while my Pfam PF00046 run on the same protein scores only 12.6 bits (below GA).

## 2. Iterative search

Pool: 72,343 proteins (UniProt proteomes of *K. phaffii* GS115 and CBS 7435, *S. cerevisiae*, *K. lactis*, *L. thermotolerans*, *C. albicans*, *C. glabrata*, *Y. lipolytica*, *O. parapolymorpha*, *D. hansenii*, *Z. rouxii*, *P. kudriavzevii*, *V. polyspora*, fetched 2026-10-17 from rest.uniprot.org) plus all curated proteins of `db/`. jackhmmer, inclusion E 1e-3, up to 6 iterations (pyhmmer returned only the final iteration's table, so per-iteration entry points are not available; the run had not converged at iteration 6).

- Iteration 1 (phmmer, E<=10): besides itself, Kp MATa2 hits *C. albicans* MTLA2 (38.3 bits, E 6e-9), Pichia kudriavzevii A0A1Z8JKM8 (31.1), Z. rouxii C5DUN6 (29.9), D. hansenii Q6BNT3 (29.0), L. thermotolerans C5DJY3 (27.2); not *K. lactis* MATA2, not Geotrichum MATA, not Yarrowia mata2.
- Iteration 6: 87 included proteins: the five above plus Ixr1, Abf2, Rox1, Hmo1, Rfg1, "HMG box domain-containing protein" of C. glabrata, V. polyspora, K. lactis, O. parapolymorpha, P. kudriavzevii, the Komagataella hypoxic-gene repressor, and many Pezizomycotina MAT1-2-1/matMc proteins of the db at 70-128 bits (about 20 of the top 40 hits; best non-self: C. glabrata Q6FN37 131 bits; Kp MATa2 self 260; *K. lactis* MATA2 118; *Z. rouxii* C5DQQ0 117, Ixr1 116).
- From Geotrichum MATA: Ixr1 131 bits, K. lactis Q6CMV5 136, then the same HMG-box set. From *K. lactis* MATA2: the same set, with Pezizomycotina MAT1-2-1 at 73-83 bits.
- From Kp MATa1 (converged, iteration 4, 94 included): MATalpha2 and MATa1 proteins of all yeasts and homeodomain proteins (Pho2, Yhp1), i.e. the homeodomain family, as expected.

That a MATa2 query reaches generic HMG-box proteins and MAT1-2-1 only after iteration, never in a single pass, is the signature of a weak, remote HMG-box signal.

## 3. Structure

AlphaFold DB v6 models fetched 2026-10-18 for the Kp, *K. lactis* and *C. albicans* MATa2, Kp MATa1, MATalpha1, MATalpha2 (`af_summary.txt`; secondary structure is a crude C-alpha i,i+3/i+4 distance rule, not DSSP; pLDDT digits: 9 = >=90). Foldseek web server (search.foldseek.com, mode 3Di+AA, PDB100) queried with each model (`fs_summary.txt`).

- **Kp MATa2** (178 aa, mean pLDDT 54): disordered N-terminus (1-46, pLDDT 30-45), then residues 47-122 high confidence (pLDDT 70-99) with three helical segments (about 47-64, 67-85, 87-118), then disordered C-terminus (123-178). *K. lactis* MATa2 (257 aa): structured 62-150, three helices. *C. albicans* MTLa2 (201 aa): structured 50-120, three helices. The same architecture in three lineages: one short three-helix domain in a long disordered protein.
- **Foldseek top hits, Kp MATa2 (query residues 32-126)**: 1LWM (sequence-non-specific HMGB protein, E 0.028, 15.6% id), 1J3C (C-terminal domain of HMGB2, E 0.019, 20.2%), 6JRP (CIC-HMG, E 0.095), 2RTU (HMGB1 A box, E 0.078), 5VWE (HMG domain of FACT, E 0.22), 1I11 (SOX5 HMG box, E 0.34). All six top hits are HMG boxes. *K. lactis* MATa2: LEF1 HMG (2LEF, E 0.044), SOX HMG (4S2Q, 4Y60, 6T78), HMGB1 box A; *C. albicans* MTLa2: POU/HMG/DNA (1GT0, E 0.13), SOX11, HMG1 domain A, UBF HMG box5. Top hits are all HMG boxes in all three, but E-values (0.02-0.65) are above the usual 1e-3 and the AlphaFold models are low confidence outside the core; read as consistent, not conclusive.
- Controls: Kp MATa1 (residues 160-216) hits homeodomains, best 1YRN (MATa1/MATalpha2 homeodomain heterodimer, E 0.0011, 48% id), 1LE8, 1AKH; Kp MATalpha2 (100-173) hits MATalpha2 homeodomain structures (1MNM, 1K61, E 1.5e-4 to 9e-4). Kp MATalpha1 hits only HMG-fold structures weakly (1J5N, 1WGF, E 2-9), as expected for an alpha box.

## 4. MATa2 alignment (MAFFT L-INS-i, columns of Kp MATa2 residues 30-130)

```
Kp_MATa2   AG-------KGLSGKKFNLS----------------THRPRNKFIISKSIFSKVIQRY----YPDLSFEE
Kl_MATA2   TI-----------KENT-----------------KKYTRPRNQFVLMRTLFNRRVNNHILQYYNKSKLEK
Ca_MTLA2   SI-----------DKMTVLS----------------SFRSRNSFIIARSILSKLLRK------NKADFKK
Dh_MTLA2   RK-------PTAKDKKTAIT--------------RGSFRTRNKFIIARNLLCRLVKHE-----KALSSTE
Cl_MTLA2   -----------------------------------KRYRSRNQFIMARQLISPFLKGV-------HNSAS
Geo_MATA   QA-------PPRKRKNSTED--------------KPIPRPTNSFICARSMLHKMFER------LNISTIE
Yl_mata2   DS---------RTNEPEYIE--------------PKVLKPRNSYILFRNATSRCSQS------IDPNSAF
NHP6A      ----------------------------------NAPKRALSAYMFFANENRDIVRSE----NPDITFGQ
SRY        SCNSKYQCETGENSKGNVQD---------------RVKRPMNAFIVWSRDQRRKMALE----NPRMRNSE
IXR1       S--------RIERRKQLKKQ---------------GPKRPSSAYFLFSMSIRNELLQQ----FPEAKVPE
Afum1-2-1  VS-VASHGEAVKTNKVTVTESFAPRGKPVGPLKAPKVPRPPNAFILYRQHHHPKIKEA----YPDYSNND

Kp_MATa2   -----NSKIASYLWKHRFGC--FDKYYEYLQIEESKFTSLV-TSKIDQLAFKAVPTE
Kl_MATA2   KMFTLTSKITSELWNESSPD--LKSYFSLLATLEENWHKYTHYCSWDRNSAQTLSMEPIELSQVRPRLIS
Ca_MTLA2   -----VSKSVSLLWGSVDNS--FKNYFEYLSILESQWHDKK-SSNLN--PISSTRTT
Dh_MTLA2   -----ISKVVSLVWKESNKS--LHDYFEYLSLLETYWYDKL-MSSYGNSTLQATATAAV
Cl_MTLA2   -----RSKLVSCIWKTKTSQ--FEKYFRYLAALDSNYFEEL-------IKFRSANNAFV
Geo_MATA   -----CSNLMSDIWKNYNDD-LFKHYFRYLSFLEACWHTSV-YPNYRYNPRKKLPSP
Yl_mata2   -----VSRVSSYIWGSGIPNPIVRRWFKYMGFFEALYHEVD-YPEY-------IPNK
NHP6A      -----VGKKLGEKWKALTPE--EKQPYEAKAQADKKRYESE-KELYN
SRY        -----ISKQLGYQWKMLTEA--EKWPFFQEAQKLQAMHREK-YPNYKYRPRR
IXR1       -----LSKLASARWKELTDD--QKKPFYEEFRTNWEKYRVV-RDAYE
Afum1-2-1  -----ISVMLGKQWKDENEE--IKTQFRNLAEELKKKHAED-HPDYHYTPRKPSERKR
```

The MATa2 core has the HMG-box start (`RP(R/T)N(K/S)(F/Y)I`, cf. SRY `KRPMNAFIVW`, Ixr1 `KRPSSAYFLF`, MAT1-2-1 `RPPNAFILY`), a conserved `S K x x S x x W` helix-2 motif (SRY `ISKQLGYQWK`, Ixr1 `LSKLASARWK`) and a conserved `x x Y F x Y L` stretch in helix 2/3 shared by the MATa2 proteins and Ixr1 (`PFYEEF`). What is missing relative to a canonical box is the helix-3 `Y..HRP` motif and the usual basic-rich linker: the box is short (about 75 residues) and lacks the conserved residues Pfam HMG_box is built from.

Core identity (columns of Kp MATa2 44-122, aligned positions in brackets):

| | Kp | Kl | Ca | Dh | Cl | Geo | Yl |
|---|---|---|---|---|---|---|---|
| Kp | 100 | 28 | 35 | 30 | 26 | 30 | 18 |
| Kl | 28 | 100 | 35 | 32 | 29 | 22 | 21 |
| Ca | 35 | 35 | 100 | 43 | 33 | 35 | 24 |
| Dh | 30 | 32 | 43 | 100 | 33 | 30 | 26 |
| Cl | 26 | 29 | 33 | 33 | 100 | 27 | 26 |
| Geo | 30 | 22 | 35 | 30 | 27 | 100 | 21 |
| Yl | 18 | 21 | 24 | 26 | 26 | 21 | 100 |

## 5. Is the Komagataella MATa2 truncated or split?

No (`an2.out`, part 3 and 4).
- CDS: 178 aa, ATG start (MTVTTFLSS...), TAG stop, single exon. The in-frame sequence upstream of the ATG has a stop 42 codons before it with no Met in between, so the N-terminus cannot extend; no intron is needed to read through. Length 178 sits within the range of the orthologues (176 *C. lusitaniae*, 201 *C. albicans*, 251 *D. hansenii*, 257 *K. lactis*, 281 Geotrichum, 291 Yarrowia); the extra length of the others is in the disordered termini, and the structured core (about 79 aligned columns) is present in full in Komagataella.
- Same protein everywhere: tblastn of MATa2 (and MATa1) gives a single full-length alignment (178/178 and 225/225 aa, no gap) in all 9 test genomes (GCA_023013455.1, GCA_023013825.1, GCA_023013855.1, GCA_023014715.1, GCA_023013895.1, GCA_023014845.1: 100% identity; GCA_030557195.1 *K. kurtzmanii* 97.8% MATa2, 99.6% MATa1; GCA_030567415.1 and GCA_000182785.1 *K. pastoris* 89.3% / 92.9%), in GCA_900010155.1 (100%), GCA_001708105.1 (89.3%), GCA_023213065.1 (89.9%), and both CBS 7435 assemblies (100%). A single HSP cannot span a stop codon, so each copy is an intact ORF. Not examined: the 9 test contigs are 2.1-2.3 kb; the HSPs show full-length ORFs but not their upstream context.
- GS115: no MATa genes (the assembly has only the MATalpha cassette; see the curation note).

## 6. A Saccharomycotina MATa2 profile

Training sets (annotated MATa2 proteins; the lineage labels are the groups left out): Komagataella Kp_MATa2; Saccharomycetaceae *K. lactis* MATA2; Serinales *C. albicans* MTLA2, *D. hansenii* MTLA2, *C. lusitaniae* MTLA2; Dipodascales Geotrichum MATA, *Y. lipolytica* mata2. Profiles built with MAFFT L-INS-i and pyhmmer `Builder` (hmmbuild defaults), searched against the 72,343-protein pool (plus the db copies of the proteins that are not in UniProt, named `ext|`). "Non-positives above" counts non-training proteins scoring higher than the held-out protein; the four that appear in every Dipodascales/Serinales/Komagataella/Saccharomycetaceae hold-out are D. hansenii DEHA2E19140p (a duplicate of the training protein: UniProt copy of Dh MTLA2), *Z. rouxii* ZYRO0C18326p, *L. thermotolerans* KLTH0F00264p and *P. kudriavzevii* A0A1Z8JKM8; ZYRO0C18326p sits next to ZYRO0C18348p and DEHA2E19140p next to DEHA2E19118p (the MATa1-type homeobox protein of the same genome), so the first three are very probably MATa2 orthologues, but that is from locus-tag adjacency only and was not checked on the genomes.

| Held-out lineage | Held-out protein | LOLO profile score (E) | rank in pool | non-positives above | 5th-best non-positive (bits) | Pfam HMG_box | curator MAT1-2-1 |
|---|---|---|---|---|---|---|---|
| Komagataella | Kp_MATa2 | 46.9 (1e-11) | 11 | 4 | 39.9 | 13.6 (rank 86) | none |
| Saccharomycetaceae | Kl_MATA2 | 27.6 (1e-05) | 11 | 4 | 25.4 | none | none |
| Serinales | Ca_MTLA2 / Dh_MTLA2 / Cl_MTLA2 | 49.8 (1e-12) / 52.3 (3e-13) / 38.8 (3e-09) | 8 / 7 / 12 | 2 / 2 / 5 | 40.2 | none | none |
| Dipodascales | Geo_MATA / Yl_mata2 | 43.6 (1e-10) / 22.4 (3e-04) | 10 / 11 | 4 / 4 | 21.3 | none / none | Geo_MATA 39.9 (rank 21) / none |
| all 7 (no hold-out) | the 7 | 139-262 | 1-8 | 0 | 33 (best non-orthologue HMG protein, A7TKB3 32.9) | | |

Caveat on the table, stated plainly: with a lineage held out, the held-out MATa2 scores 22-52 bits, and the best genuine HMG-box background protein scores 21-40 bits (5th non-positive), so the margin is clear for *C. albicans*/*D. hansenii* and Geotrichum (E <= 1e-10), thin for Komagataella (46.9 vs 39.9), *K. lactis* (27.6 vs 25.4) and *Y. lipolytica* (22.4 vs 21.3). The profile is specific only when several lineages are in training, and the training set is just 7 proteins whose MATa2 identity is annotation-based. The same exercise for MATa1 (homeodomain, 10 proteins in 5 lineages): held-out MATa1 scores 39-88 bits (E 1e-9 to 1e-24) and ranks 9-19 in the pool; the non-positives above are other MATa1 orthologues (Ogataea W1QBL1, D. hansenii DEHA2E19118p) and the best generic homeodomain protein scores 45-55 bits, so Saccharomycotina MATa1 is separable but the Pfam homeodomain profile also works (34-56 bits).

## 7. What this means for gene_class and idiomorph_class (recommendation only, nothing changed)

- **`gene_class`**. MATa2 is an HMG-box protein by three consistent lines (iterative links to HMG-box proteins and MAT1-2-1; three-helix AlphaFold core with HMG-box Foldseek hits in three lineages; alignment motifs) but by none of the standing profiles. The `HMG_box` class for MATkp `MATa2` and for the A-idiomorph MATA gene of `MATdc` is therefore biologically defensible, yet cannot be confirmed by a rule that asks for a Pfam or MAT1-2-1 hit; `homothallic_candidate` in `classify_locus` uses declared classes, not profile hits, so declaring it is mechanically safe. Recommend: declare `HMG_box` on MATa2 (MATkp) and on Geotrichum MATA (MATdc; it already has the curator MAT1-2-1 hit at 39.9 bits) and on Yarrowia `mata2` (MATyl; HMG per the literature quoted in the db), but only after the curator accepts the structural and iterative evidence here as the bar; keep MATa1/`mata1` unset (homeodomain has no class in the schema).
- **A class that fits both MATkp and MATdc would be a profile, not a rule**: build the Saccharomycotina MATa2 profile (section 6) into `classifiers/` the way `MAT_core` was built, and use it only as supporting evidence, because its margin over background HMG proteins is thin for single divergent lineages.
- **`idiomorph_class`**. MATalpha -> MAT1-1 is supported (alpha box above GA in Komagataella and every other yeast tested). MATa -> MAT1-2 is supported by the MATa2 HMG-box evidence above, not by Pfam. By the same logic Yarrowia A (mata1 homeodomain + mata2 HMG) is the MATa-type idiomorph and B is alpha-type (matb1 weak alpha-box hit, 15.6 bits; matb2 homeodomain per InterPro), and Geotrichum A (MATA, HMG, an MATa2 orthologue by the LOLO profile, 43.6 bits held out) / B (MATB, alpha box 35.5) line up the same way. This agrees with the earlier MATdc/MATyl A/B mapping and sharpens it: Yarrowia mata1 is homeodomain (34.0 bits) and mata2 is an a2-family protein. Leave both families `unassigned` until the curator accepts a structural/iterative evidence bar; if accepted, A -> MAT1-2 and B -> MAT1-1 for all three families.
- **Detection**: the lack of a profile match does not hurt detection of MATa. `MATkp` finds MATa1/MATa2 at 100% in the test contigs by protein homology to the CBS 7435 record (see the curation note); the profile matters only for naming the class.

## Not verified

- Foldseek E-values (0.02-0.65) are borderline and the AlphaFold models for MATa2 are low confidence outside the 75-90 residue core; the "three helices" is from a crude C-alpha distance rule, not DSSP (no DSSP in the environment).
- Identity of the *Z. rouxii*, *L. thermotolerans* and *P. kudriavzevii* proteins as MATa2 orthologues (locus-tag adjacency and scores only; not added to the training set).
- Per-iteration entry points of jackhmmer (only the final iteration was returned).
- That the 7-protein MATa2 training set is complete; Yarrowia `mata2`/Geotrichum MATA membership in the MATa2 family is by this analysis, not by a published domain annotation. I did not locate a source for the HMG/homeodomain content of K. lactis a1/a2 beyond the InterPro and UniProt records fetched here.
- Whether Pfam HMG_box's sensitivity could be recovered with a lower threshold: HMG_box E<=1000 returns only 89 proteins in the pool, none of them Saccharomycotina MATa2 except Kp at 13.6 bits.
- InterPro and Foldseek/AlphaFold results are cited as services consulted on 2026-10-17/18; no PMID was opened for them.

## 8. The profile as built into the database (deterministic build), and a correction to section 6

Built 2026-10-18 into `db/Ascomycota/classifiers/MATa2_saccharomycotina/` (`training.faa`, `MATa2_saccharomycotina.hmm`, `manifest.yaml`) by `scripts/build_mata2_saccharomycotina_hmm.py` with the `classifier_build` functions used for `MAT_core` (single-threaded MAFFT L-INS-i on id-sorted input, Builder seed 42, fixed HMM header). It is supporting evidence only: no `order.yml` entry, vote, gate or classification step references it, and a test checks both that it builds reproducibly and that nothing references it.

The deterministic build differs slightly from the exploratory build of section 6 (different MAFFT thread/order handling, and background defined as the best protein that is not a training protein and not one of three probable unlabelled MATa2 orthologues). Its numbers replace those of section 6 and of the `order.yml` comment: held-out lineage vs best background HMG-box protein (bits): Komagataella 43.9 vs 37.7 (+6.3); *K. lactis* 26.7 vs 25.4 (+1.3); *Y. lipolytica* mata2 22.4 vs 21.3 (+1.1); *C. albicans* 50.2 vs 41.3 (+8.9); *D. hansenii* 56.7 vs 41.3 (+15.4); *C. lusitaniae* 39.5 vs 41.3 (-1.8, below the background); Geotrichum MATA 43.9 vs 21.3 (+22.6). So the margins are thin for Komagataella, *K. lactis* and Yarrowia, negative for *Clavispora* and clear only for Geotrichum and *Debaryomyces*; the claim in section 6 that the profile "separates" is true only for some lineages.

## 9. Effect of gene_class on both MATa2 (HMG_box) and MATalpha1 (alpha_box): what `classify_locus` does for Komagataella

Detection with the three records, `gene_class` on MATa2 and MATalpha1, 15 genomes (the 9 test genomes, CBS 7435 on both assemblies, GS115, GCA_001708105.1, GCA_023213065.1, GCA_900010155.1), `--taxid` own, `MATPREDICT_HTML=0` (HPCC, `work/det3`). Result, exactly as reported:

- **No genome gets `homothallic_candidate`.** Every call is `mat_locus` (core gene plus flank: CBS 7435 both assemblies, GS115, K. pastoris GCA_001708105.1, K. ulmi GCA_023213065.1, and the test contigs that reach SLA2) or `idiomorph_gene_only` (MATa1 + MATa2 on contigs with no flank gene: GCA_023013895.1, GCA_023014845.1, GCA_000182785.1, GCA_900010155.1).
- **Both cassettes of one genome are reported as two separate calls**: CBS 7435 LT962479.2 `MATkp MATa` 1771-7997 and `MATkp MATalpha` 136012-141828 (both high, `mat_locus`); FR839631.1 1666-7892 and 135891-141707; GCA_001708105.1 CP014587.1 1569316-1575518 and 1438280-1444097; GCA_023213065.1 JAKWAZ010000004.1 1569573-1575776 and 1439297-1445114; GS115 only the alpha cassette.
- **Why**: the different-domains path in `classify_locus` compares core genes of ONE call (one cluster, `max_homothallic_separation_bp`, same contig). The two cassettes are ~135 kb apart, so each call holds the genes of a single idiomorph and no pair of different-class genes ever meets in one call; the declared classes (alpha_box on MATalpha1, HMG_box on MATa2) therefore change nothing for Komagataella today. MATa1/MATalpha2 are unclassed, so even within one cassette there is only one classed gene.
- **`idiomorph_class` is `unassigned` on every call; `two_idiomorphs` is `[]` on every report**: `two_idiomorphs_report` is off for MATkp (asserted by `tests/detect/test_matkp_family.py`; not enabled, still an open question). Enabling it would be the mechanism that states "both idiomorphs called in this genome" for K. phaffii; a `homothallic_candidate` label for two cassettes 135 kb apart would need either that report or a larger `max_homothallic_separation_bp` (which would also admit both cassettes into one cluster; not tried, and not wanted without a ruling).

## 10. Experiment: do the MATkp references find the MAT locus in Ogataea? (throwaway db; nothing committed to `db/`)

Setup. NCBI Taxonomy (efetch, 2026-10-18): *Ogataea* genus 461281 (rank genus, Pichiales; Pichiaceae 1156497). Two throwaway copies of the branch db on the HPCC (`/bigdata/stajichlab/jstajich/projects/kom/ogataea_scope_test/dbA`, `dbB`) differ from the branch db only in the MATkp `taxonomic_scope`: A = `[460517, 461281]` (genus Ogataea added), B = `[460517, 1156497]` (Pichiaceae). `detect --taxid <own taxid>` (`MATPREDICT_HTML=0`, v0.6.2 code of this worktree) on all 112 Ogataea genomes of the v0.6.0 campaign (Pichiales; species Ogataea; no *Hansenula*/*Pichia angusta* names occur in the list), with the unchanged branch db ("base", scope genus Komagataella, so routing is `phylum_fallback`) and with scope A (routing `lineage`). Scope B was run only on 6 non-Ogataea Pichiaceae as controls (for an Ogataea genome B routes exactly like A). No Ogataea record exists in the db; the MATkp references are the three Komagataella/GS115 records.

Independent reference. Deposit `AJ617305.1` (*O. angusta* CBS 4732, 6,239 bp, PMID 14745027, Butler et al. 2004; abstract: in *P. angusta* "MATalpha2, MATalpha1, and MATa1 genes adjacent to each other on the same chromosome"; genes yol077C partial, dic1, matalpha2, matalpha1, mata1, sla2 partial). Maekawa & Kaneko 2014 (PMID 25412462, abstract): *H. polymorpha* (*O. polymorpha*) has two MAT loci about 18 kb apart on one chromosome, MAT1 (alpha1, alpha2) and MAT2 (a1), flanked by SLA2 and SUI1 where active, and switches by inversion of the intervening region; Hanson et al. 2014 (PMID 25349420, abstract): 19 kb invertible region. tblastn of the five deposit proteins and the Komagataella proteins against four reference assemblies (`ogref.out` on the HPCC):

| Assembly | Deposit layout (all deposit proteins 85-100% identical) | Komagataella proteins on the same layout |
|---|---|---|
| *O. polymorpha* GCF_001664045.1 (NW_017264700.1) | SLA2 899,712-901,031; a1 901,050-901,592; alpha1 901,577-902,161; alpha2 902,170-902,688; DIC1 903,810-904,670 (all within 5 kb, 100% to the deposit); a second a1 fragment 920,426-920,593 and SLA2 copy 920,791-922,557 about 19 kb away | alpha1 35% (130/190 aa), alpha2 39% (82/173 aa), a1 55% (56/225 aa, homeodomain only), DIC1 67% (274/305), SLA2 44%; **MATa2: no hit at E 1e-5** |
| *O. parapolymorpha* DL-1 GCF_000187245.1 (NC_027862.1) | SLA2 424,263-425,582; a1 425,601-426,161; alpha1 426,146-426,724; alpha2 426,740-427,258; DIC1 428,386-429,246 (deposit proteins 85-97%); a1 fragment 445,001-445,207 (~19 kb away) | same pattern; no MATa2 |
| *O. polymorpha* GCA_000787595.2 (4329) | DIC1, alpha2, alpha1 on DF933571.1:913,308-916,401; a1 and SLA2 on DF933585.1:15-2048 (contig break) | same pattern; no MATa2 |
| *O. polymorpha* GCA_019207825.1 | alpha1/alpha2/DIC1 on JAHLVP010000020.1:2-3095; a1/SLA2 on JAHLVP010000036.1:2-1864 (contig break) | same pattern; no MATa2 |

So the deposit layout is real in the reference genomes (a1, alpha1, alpha2 contiguous between SLA2 and DIC1, with a second SLA2/a1 copy about 19 kb away that matches the "about 18-19 kb apart" of the two papers; I did not reconcile the deposit's contiguous a1-alpha1 arrangement with the two-loci model beyond this). No Komagataella-like MATa2 exists in these Ogataea genomes (not in the deposit either). Not verified: that this is a true absence rather than a divergent a2 below the tblastn cut-off.

Results, 112 Ogataea genomes (scope A, v0.6.2; full per-genome table in `ogsum.out` on the HPCC):

| | Base db (genus-only scope) | Scope A (MATkp covers Ogataea) |
|---|---|---|
| Genomes with any call | 112 | 111 (1 none) |
| Families called | MATsc MATalpha 93, MATsc MATa 38, MTL alpha 42, MATyl A 17, Pezizomycotina MAT 28 (MAT1-1/MAT1-2), a mix of medium and high confidence calls of several families per genome | only MATkp |
| Both idiomorphs | n/a | 74 genomes with a MATa and a MATalpha call; **48 calls contain both a and alpha gene sets in one call** (cassettes within the 25 kb cluster gap, e.g. both *O. parapolymorpha* references, 24 kb call labelled `MATa` medium) |
| MATkp MATalpha calls | | 94 calls in 93 genomes (high 68, medium 26) |
| MATkp MATa calls | | 94 calls in 92 genomes (high 6, medium 88); 77 of all MATkp calls are `partial_locus` (mostly MATa1 + a flank gene only) |
| Hit identity over genomes (median; coverage median) | | MATa1 54% (25%), MATalpha1 41% (41%), MATalpha2 42% (32%), DIC1 67% (90%), SLA2 44% (41%), MATa2 31% (24%, 51 calls) |

The 51 MATa2 "hits" are 18-38% identity over 19-72% of the protein (bit scores 20-58), i.e. not a credible MATa2; they sit in the alpha call or the merged call and are most likely an unrelated HMG-box fragment near the cassette (not examined). The MATa1 hit covers about 25% of the Komagataella MATa1: it is the homeodomain only.

Scope B on 6 non-Ogataea Pichiaceae (base calls MATsc/MATalpha, PM, etc. replaced by MATkp): *P. kudriavzevii* MATalpha high, *Ambrosiozyma kashinagicola* MATalpha high, *C. boidinii* MATalpha medium, *Saturnispora bothae* MATalpha medium, *P. manshurica* MATa medium, *Brettanomyces custersianus* MATa low. A Pichiaceae-wide scope would call all of these on 31-54% identity references; none of them was checked against a deposit.

Conclusion (recommendation only). The MATkp references do locate the cassettes in Ogataea (111 of 112 genomes, both idiomorphs separately in 74, SLA2 and DIC1 flanks at 44% and 67%), and the loci agree with the deposit and the reference-genome layout. But the evidence is not good enough for a widened scope alone: (1) identities are 31-54% with coverage of 25-41% on the core genes, so every call rests on a distant partial match (alpha high in 68 of 94 calls, MATa high in 6 of 94); (2) Ogataea has no MATa2, so a MATa call is MATa1 plus flank only (77 `partial_locus` calls overall) and the family's MATa2 requirement cannot be met; (3) in 48 calls the two cassettes (about 19 kb apart in the references) fall into one cluster and the idiomorph label is one letter, which misreports a normal two-cassette genome; (4) the call loses the other families' fragments but gains nothing in confidence over them. Ogataea therefore needs its own records first: at least O. polymorpha/O. angusta (`AJ617305.1`, a1/alpha1/alpha2/DIC1/SLA2 deposit proteins 85-100% identical in all four reference assemblies; strain-matched assembly needed under B10) and an *O. parapolymorpha* DL-1 record (GCF_000187245.1), declared with the same gene names (a1, alpha1, alpha2) so the same family could take them with MATa2 optional; widening the scope should follow the records, and the merged-cluster labelling (point 3) needs a rule decision first. Not examined: 2-3 Ogataea genomes where the MATa call spans a contig break (a1 and alpha cassettes on different contigs), and whether the weak MATsc/MTL calls of the base run sit on the same loci (not compared).

## 11. Ogataea records added (curator approval, 2026-10-18)

Two accepted records in the MATkp family directory (`db/Ascomycota/Pichiales/`), both `system: pseudohomothallic`, family scope unchanged (genus Komagataella 460517; the family comment says the Ogataea records are present but out of scope pending a curator decision):

- `460523_cbs-4732_MATkp_MATa-MATalpha`: deposit AJ617305.1 (dic1, MATalpha2, MATalpha1, MATa1, sla2 partial), taxid 460523 because NCBI Taxonomy lists CBS 4732 as an *O. polymorpha* name (the deposit says *O. angusta* 870730). **B10 not met**: the deposit is 100.000% over 4,962 bp identical to Hanpo2 NW_017264700.1:899712-904673 (and to BY4329 and GCA_019207825.1), but no source I could open says those strains are CBS 4732, so `assembly_accession` is null. All five genes validate 100%/100% (compute node, `scripts/curate_db_validate.slurm`, `taxonomy_current: true`).
- `1005962_dl-1_MATkp_MATa-MATalpha`: *O. parapolymorpha* DL-1, NC_027862.1:422526-429249 (GCF_000187245.1, `sequence_link`). SLA2, MATa1, MATalpha2 and DIC1 are RefSeq models (validated 100%/100%); MATalpha1 is not annotated in RefSeq and was derived from the deposit protein (190 aa ORF, Pfam PF04769 49.6 bits). The second MATa1/SLA2 copy about 19 kb away (HPODL_04024/04025) is described in the note, not recorded.
- **MATa2 is not recorded as present or absent.** tblastn at E 1e-5 of seven MATa2 proteins finds nothing in the Ogataea cassette regions, but at E 10 the Komagataella MATa2 aligns twice (E 0.009) at DL-1 444128-444363, and two overlapping minus-strand ORFs there (443765-444202, 444043-444387) score 14.0 and 18.1 bits with the 7-protein MATa2 profile (E 1e-3, 5e-5) and align to the HMG-box start and helix 2; a similar ORF pair occurs in Hanpo2. This may be a divergent or split MATa2 beside the second MATa1 or an unrelated HMG-like ORF. Unresolved, and dependent on the cut-off.
- Effect on detection: Komagataella calls on the 15 genomes are identical with the two extra records in the family (the extra references change nothing there); Ogataea genomes are still routed by `phylum_fallback` (MATsc/MTL/MAT calls, as in the base run of section 10) because the scope was not widened. Whole suite on the HPCC (`pytest tests`): 1265 passed, 4 skipped, 0 failed on the second run; the first run had the known intermittent `test_pipeline.py::test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window` failure (1264 passed, 1 failed).

## 12. Ogataea moved to its own family MATog (curator decision, 2026-10-18)

**UPDATE (later the same day, curator): the genus scope described below is NOT shipped.** MATog is now records plus family declaration only and routes no genome; the effects measured below with the scope [461281] are those of a throwaway db and are kept as the first measurement. See 12.1 (mechanism and verification) and 12.2 (leave-own-record-out).

Section 11 placed the two Ogataea records in MATkp; they are now in the new family `MATog` (record ids `460523_cbs-4732_MATog_MATa-MATalpha`, `1005962_dl-1_MATog_MATa-MATalpha`, regenerated by propose, validate on a compute node per `scripts/curate_db_validate.slurm` (every gene 100%/100%, `taxonomy_current: true`), accept, build-gff). MATkp is unchanged (`taxonomic_scope` [460517]).

MATog declaration: `fallback_searchable: false`, `taxonomic_scope: [461281]` (genus *Ogataea*; verified 2026-10-18 with taxonkit: `name2taxid` gives 461281 and `lineage -r` gives rank genus under Pichiales; Pichiaceae; Ogataea), genes MATa1, MATalpha1 (`gene_class: alpha_box`, Pfam PF04769 50.1 bits on the deposit protein, 49.6 on the DL-1 ORF), MATalpha2, sla2 (`sla2_homolog` by precedent), dic1, and MATa2 optional, unresolved and without class, with the E 10 / 14.0 and 18.1 bit evidence in the family comment; `idiomorph_class` unassigned, `two_idiomorphs_report` off, records pseudohomothallic. I chose a genus scope (the scope-only family routes exactly the genus Ogataea and no other genome, whereas a family with no scope at all would be searched on every phylum fallback; `fallback_searchable: false` already prevents that, so the genus scope is what makes the family searched at all) because the curator asked for a scope-only family for the genus.

Effect (HPCC, `detect --taxid <own taxid>`, `MATPREDICT_HTML=0`, v0.6.2 code, branch db):
- **Komagataella unchanged**: the 15 Komagataella genomes give identical calls with MATog in the db as without (15/15 identical; 0 MATog calls in them).
- **112 Ogataea genomes**: routing changes from `phylum_fallback` (base: 220 calls from MATsc 132, MTL 43, MAT 28, MATyl 17) to `lineage` with only MATog (188 calls). 111 of 112 genomes get a call; 73 get a MATa and a MATalpha call as separate calls, 14 calls contain both cassettes in one call; confidence is high for 82 of 83 MATa calls and 97 of 105 MATalpha calls; 185 `mat_locus`, 1 `partial_locus`, 2 `idiomorph_gene_only`. No genome is called by a non-MATog family any more.
- **Identity** (gene hits, medians; the records' own genomes are included, which raises the first row): the 37 *O. polymorpha/angusta/parapolymorpha* genomes (the species of the two records) MATa1 100%, MATalpha1 100%, MATalpha2 91.5%, SLA2 95%, DIC1 100%, 26 with both idiomorphs, no genome uncalled; the 75 other Ogataea genomes MATa1 66%, MATalpha1 45%, MATalpha2 55%, SLA2 55%, DIC1 75%, 47 with both idiomorphs, 1 uncalled, 9 calls at medium confidence. So within the polymorpha complex the family is a close match (deposit and DL-1 are the references), and for the other 75 species the match is distant (alpha1 45%) but still called high in 113 of 122 calls.
- **Named references**: *O. polymorpha* GCA_019207825.1 MATalpha (alpha1, alpha2, dic1) and MATa (a1, sla2) on two contigs, both high; Hanpo2 and *O. parapolymorpha* DL-1 one call each spanning both cassettes (a1, alpha1, alpha2, sla2, dic1; 24 kb, labelled MATalpha; the merged-cluster labelling question of section 10 point 3 remains); BY4329 two calls on two contigs.
- Not recomputed for this section: leave-own-record-out on DL-1 and the deposit-matching genomes; their 100% hits are self-hits.

Whole suite on the HPCC (`pytest tests`): 1267 passed, 4 skipped, 0 failed. Tests changed: `test_fallback_searchable.py` roster now includes MATog; `test_matkp_family.py` asserts MATkp scope, the MATog declaration (scope [461281], scope-only, MATa2 optional without class, MATalpha1 alpha_box, no two_idiomorphs report) and the two records.

### 12.1 Scope held back: mechanism and verification

MATog's `taxonomic_scope` holds the sentinel `[-461281]` (the negative of the intended genus taxid 461281). The schema requires the key (minItems 1, integer items, no minimum), no lineage or genome taxid can contain a negative number, and `fallback_searchable: false` keeps the family out of phylum fallback and explicit-phylum runs; so nothing routes to it. The family comment says the scope is held back pending 12.2 and gives the release step (replace the sentinel by 461281). `--exhaustive` still searches every family, MATog included. MATkp is untouched ([460517]). Verified on the HPCC with the branch db: on the 112 Ogataea genomes the routing (`phylum_fallback`) and the calls (`family`, idiomorph, contig, span, confidence, genes) are identical to the run with the main-behaviour db in 112 of 112 genomes, 220 calls in both (MATsc 132, MTL 43, MAT 28, MATyl 17); `tests/detect/test_matkp_family.py` also asserts that `route()` for an Ogataea lineage returns no MATog and no MATkp.

### 12.2 Leave-own-record-out analysis (throwaway dbs, scope widened to 461281 in the copies only)

Setup: three copies of the branch db under `ogataea_scope_test/loro/` with MATog scope [461281] and (dep) only the deposit record 460523, (dl1) only the DL-1 record 1005962, (both) both; `detect --taxid <own>` (v0.6.2 code, `MATPREDICT_HTML=0`) on all 112 Ogataea genomes with each, plus the earlier run with only the three Komagataella records as reference (MATkp scope widened, "A"). Groups: *O. parapolymorpha* (7 genomes, includes DL-1), *O. polymorpha/angusta* (30, includes the deposit-matching Hanpo2, BY4329, GCA_019207825.1; one hybrid), other Ogataea (75). `loro_an.py` and its output are in `ogataea_scope_test/` on the HPCC.

**What "correct" can mean.** There is no independent truth for the MAT locus of any of these genomes. Self-hits are not evidence: DL-1 with the DL-1 record, and Hanpo2, BY4329 and GCA_019207825.1 with the deposit record, find their own sequence (100%) and are shown only as "both" rows. Three weaker checks are used, each stated for what it is: (1) positional check: in four genomes whose cassette positions are known from the OTHER record's copy (DL-1: the record coordinates; the three deposit-matching assemblies: tblastn of the deposit proteins at 100%), does the foreign-only reference place each gene at the right place? (2) agreement between references: does a call made with one reference overlap a call made with a different one (deposit vs DL-1 vs Komagataella)? Agreement of the three is evidence that all find the same locus, but the deposit and DL-1 are the same clade, so only the Komagataella agreement is lineage-independent; it shows locus identity, not mating-type identity. (3) anchoring: the call contains a flank gene (SLA2/DIC1, ordinary conserved genes, adjacent to the MAT genes in all references) plus a core gene.

Per case (calls; uncalled genomes; genomes with both idiomorphs; calls that merge both cassettes; calls anchored by core + flank; confidence; median identity/coverage % of gene hits):

| Case | Genomes | Calls | Uncalled | Both idiomorphs | Merged | Anchored | Confidence | MATa1 | MATalpha1 | MATalpha2 | SLA2 | DIC1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| a: *O. parapolymorpha*, deposit only | 7 | 12 | 0 | 5 | 2 | 12 | high 12 | 91/100 | 89/100 | 93/100 | 96/100 | 97/100 |
| a': same, both records (self-hit for DL-1) | 7 | 12 | 0 | 5 | 2 | 12 | high 12 | 100/97 | 100/100 | 100/100 | 100/49 | 100/100 |
| b: *O. polymorpha/angusta*, DL-1 only | 30 | 54 | 0 | 21 | 9 | 54 | high 54 | 91/96 | 90/100 | 84/96 | 95/49 | 97/100 |
| b': same, both (self-hits in 3 genomes) | 30 | 54 | 0 | 21 | 9 | 54 | high 54 | 100/97 | 99/100 | 84/96 | 95/49 | 100/100 |
| c1: other Ogataea, deposit only | 75 | 122 | 1 | 48 | 3 | 120 | high 112, medium 10 | 67/35 | 45/79 | 55/36 | 51/99 | 75/96 |
| c2: other, DL-1 only | 75 | 122 | 1 | 47 | 3 | 120 | high 114, medium 8 | 66/37 | 44/82 | 55/27 | 55/65 | 75/95 |
| c3: other, both | 75 | 122 | 1 | 47 | 3 | 120 | high 113, medium 9 | 66/52 | 45/82 | 55/27 | 55/64 | 75/96 |
| d: other, Komagataella records only | 75 | 122 | 1 | 47 | 21 | 120 | high 43, medium 79 | 54/25 | 36/54 | 41/32 | 56/64 | 67/90 |
| d': *polymorpha/angusta*, Komagataella only | 30 | 54 | 0 | 22 | 20 | 54 | high 26, medium 28 | 55/25 | 47/35 | 45/32 | 44/41 | 67/90 |
| d'': *parapolymorpha*, Komagataella only | 7 | 12 | 0 | 5 | 7 | 12 | high 5, medium 7 | 55/25 | 42/41 | 45/32 | 43/41 | 67/90 |

Positional check (genes at the right place, foreign reference only; identity/coverage in brackets): DL-1 with the deposit only 5/5 (SLA2 96/100, MATa1 91/100, MATalpha1 89/100, MATalpha2 93/100, DIC1 97/100); BY4329 and GCA_019207825.1 with DL-1 only 5/5 each (MATa1 91/94, MATalpha1 90/100, MATalpha2 92/96, SLA2 95/49, DIC1 97/100); Hanpo2 with DL-1 only 4/5 (SLA2 assigned to the second SLA2 copy about 19 kb away, not to the cassette copy, so the call is at the right locus but that one gene is not at the position of the deposit copy). With the Komagataella records only the same four genomes: 5/5, 4/5, 5/5, 5/5 at MATa1 55/25, MATalpha1 42-47/32-41, MATalpha2 45/32, DIC1 67/90 (gene positions are right, but only the homeodomain of MATa1 and a quarter to a third of MATalpha2 are matched).

Agreement between references in the 75 other genomes (calls of reference X overlapping a call of reference Y, same idiomorph, by median core-gene identity of the call): deposit vs DL-1 <40%: 3/3; 40-55%: 62/62; 55-70%: 9/9; >=70%: 47/48. DL-1 vs deposit: 3/3, 61/61, 11/11, 46/47. Deposit vs Komagataella: 3/3, 61/62, 8/9, 47/48. DL-1 vs Komagataella: 3/3, 61/61, 9/11, 47/47. Komagataella calls overlapping a call of both Ogataea records: <40%: 68/68; 40-55%: 45/47; 55-70%: 7/7.

Confidence against identity (calls in the 75 other genomes, median core-gene identity bin; deposit alone / DL-1 alone / Komagataella alone): <40%: 1 high 2 medium / 3 high / 17 high 51 medium; 40-55%: 54 high 8 medium / 53 high 8 medium / 26 high 21 medium; 55-70%: 9 high / 11 high / 7 medium; >=70%: 48 high / 47 high / none. Confidence is therefore not a function of identity alone: the same 40-55% bin is 87% high with an Ogataea reference and 55% high with a Komagataella one; what raises it is the number of genes found and their anchoring (120 of 122 calls carry a core and a flank gene in every case).

Reading, and a recommendation (not a change):
1. **Within the polymorpha complex** (37 of the 112 genomes), a foreign Ogataea record alone finds the cassette at 84-97% identity, high confidence, all genes at the right place in the genomes where the layout is known; a record of one species is a good reference for the other.
2. **For the other 75 Ogataea genomes** (alpha1 45%, a1 66% over only 35-52% of the protein, alpha2 55% over 27-36%) the calls land on the same locus whichever reference is used: Ogataea deposit, DL-1 and the Komagataella records agree on 61-62 of 62 calls in the 40-55% bin, and 120 of 122 calls are anchored by SLA2/DIC1 at 51-75% identity. That is good evidence that the call identifies the cassette locus (position), and the only lineage-independent evidence available.
3. **It is not evidence that 'high' is justified at the gene level.** At 45% identity and 27-36% coverage of MATalpha2 and 35% of MATa1 (homeodomain only) the MAT gene models are partial; the idiomorph label rests on those partial hits; 3 calls merge both cassettes (21 with Komagataella references) and the idiomorph label of a merged call is a single letter; there is no independent truth about which idiomorph(s) the genome carries. 'High' at 45-75% identity is justified as "locus found, anchored and reproducible across three references", not as "both genes modelled, idiomorph assigned".
4. **Scope.** Supported now: the *O. polymorpha*/*angusta*/*parapolymorpha* complex (high confidence, 84-100% identity, positions checked). For the rest of the genus the evidence supports a locus-level call; a genus scope is defensible only if the report states the call is a locus detection (or if confidence is capped at medium below a coverage/identity floor, which is a rule decision for the curator). The scope stays held back as instructed; I did not choose between complex-only and genus.
5. Caveats: the 'other' group is dominated by a few species with several assemblies (23 *O. haglerorum*), so the 75 are not 75 independent points; the agreement check cannot detect a locus that all three references miss; Hanpo2/BY4329/GCA_019207825.1 are three assemblies of closely related strains of one species; the record-less alternatives were not tested (a second Ogataea clade record, e.g. from *O. haglerorum*, would show how much identity would rise for the distant genomes).

Whole suite on the HPCC (`pytest tests`, `short` partition, 2 CPUs, 8 GB) after these changes: 1268 passed, 4 skipped, 0 failed on the second run; the first run was 1267 passed, 1 failed (the known intermittent `test_pipeline.py::test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window`).

### 12.3 Complex-only scope (curator decision, 2026-10-18; supersedes the sentinel of 12.1)

MATog's `taxonomic_scope` is now `[460523, 870730, 1005962, 2053387]`, the species-level taxids of *O. polymorpha*, *O. angusta*, *O. parapolymorpha* and the natural hybrid *O. polymorpha x O. parapolymorpha*; the genus 461281 is not in it and `fallback_searchable` stays false; MATkp stays [460517]. Verified with taxonkit under Slurm (`lineage -r -t`, `reformat`): all four are rank species under *Ogataea* (461281; Pichiaceae; Pichiales); the DL-1 strain 871575 (rank strain) lies under 1005962 and routes by lineage. These are all the taxids carried by the 37 complex genomes of the 112 (460523 x18, 870730 x11, 1005962 x6, 871575 x1, 2053387 x1). The family comment says the genus scope is deferred until a record from a distant Ogataea clade exists and that the "locus detection" label versus confidence cap question is open.

Rerun of the 112 Ogataea genomes with the real scope (branch db, v0.6.2 code): the 37 complex genomes route to MATog (36 `direct`, 1 `lineage`, the DL-1 strain) and are called by MATog only (66 calls, identical in all 37 genomes to the run of 12.2 with both records, so the numbers of 12.2 rows a'/b' apply); the other 75 genomes stay `phylum_fallback` with calls identical to the main-behaviour run in 75 of 75 (144 calls, no MATog call). Total 210 calls = 144 + 66, against 220 on main (the 76 MATsc/MTL/MAT/MATyl fragment calls of the 37 complex genomes are replaced by 66 MATog calls). `tests/detect/test_matkp_family.py` asserts the scope, that the four complex taxids (including the DL-1 strain and the hybrid) route to MATog and not to MATkp, and that another Ogataea species (O. glucozyma 36023) routes to neither.

One detail of 12.2 stated explicitly: for Hanpo2 (GCF_001664045.1) the MATog call made from the DL-1 record alone assigns SLA2 to the second SLA2 copy about 19 kb away (NW_017264700.1:920791-922557) and not to the cassette copy (899712-901031) of the deposit, so the call is at the right locus but SLA2 is not at the cassette position (4 of 5 genes at the right place).
