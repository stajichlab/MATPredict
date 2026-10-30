# MATa2 support and a1 / alpha2 profiles (2026-10-30)

Question (curator): (a) is there support for the MATa2 genes of the new tier-3 records, (b) can profiles tell MATa1 from MATalpha2 (both homeodomain), and how many copies are needed.

## Data
`inventory.py`, `fetch_public.py`, `fetch_public2.py`, `build_sets.py` collect the proteins: the curated `db/` Saccharomycotina proteins plus titled NCBI Protein records (taxid 147537), **one per species per class**, mislabeled/partial/hypothetical removed (`data/set_meta.tsv`).
**a1 29 species, alpha2 25, a2 14, alpha1 12** (Saccharomycetaceae, Serinales, Pichiaceae, Dipodascales; Saccharomycetaceae dominates alpha2: 17 of 25). Profiles: mafft `--localpair --maxiterate 1000`, pyhmmer `Builder` (seed 42); `loco.py`, `a2_loo.py`, `ablation.py` run on the HPCC; `hmm/*.hmm` are the results. Pfam is `pfam_scan.py`; scoring is `score_tests.py`.

## (b) a1 vs alpha2: the data support separate profiles
- **Leave-one-clade-out** (`data/loco.tsv`): 53 of 54 correct; leave-one-species-out also 53 of 54. Typical margin (own profile minus the other) 65 bits for a1 and 75 for alpha2; the second-lowest margin is 27 bits for both.
- The one error is a labelling error, not a model error: *Hyphopichia burtonii* XP_020074339.1 "Mtla1p" scores a1 19 vs alpha2 115 (likely an NCBI annotation error). A second near-zero case (*Ogataea thermomethanolica* QFP92360.1 "ALPHA2", 0.0 vs 0.0) is probably truncated.
- **Copies needed** (`data/ablation.tsv`, 8 random draws per k): accuracy on the held-out rest is 97-98% already with **2 sequences per class** and does not rise with more; the residual 2% is those two sequences. So the classes are separable with few copies; what limits use is clade coverage (alpha2 is nearly all Saccharomycetaceae and Pichiaceae, only 2 Serinales; no Ascoideales/Phaffomycetales/Dipodascales alpha2), not copy number.
- Verdict: **build them** (`hmm/a1.hmm`, `hmm/alpha2.hmm`; a2 and alpha1 too), as supporting evidence, referenced by no gate until the curator decides, like the MATa2 profile. Report the score difference, not a hard call.

## Scoring our ORFs (`data/test_scores.tsv`, `data/test_rank.tsv`; background = the S. schoenii proteome, 3,067 proteins)
| protein | a1 | alpha2 | a2 | alpha1 | reading |
|---|---|---|---|---|---|
| P. fabianii MATa1 | 111 | 9 | | | a1 |
| B. californica MATa1 | 104 | 23 | | | a1 |
| S. schoenii right homeodomain 536,046-536,944 | 85 | 27 | | | a1 (MATa1, as recorded) |
| S. schoenii mirrored copy 529,258-529,749 | 60 | 20 | | | a1 |
| **S. schoenii left homeodomain p5 (471,961-472,489)** | 86 | 29 | | | **a1, not alpha2** |
| **S. schoenii "HMG-like" ORF p3 (470,349-470,951)** | 9 | **81** | -2 | | **alpha2, not a2** |
| B. californica 620,080-620,418 (113 aa, frame +2) | 4 | **37** | | | alpha2 (partial homeodomain) |

## (a) MATa2 support
- **a2 profile** (14 sequences; `hmm/a2.hmm`) and the existing Saccharomycotina MATa2 profile, background maximum 19 and 17 bits:
  - *P. fabianii* MATa2 **73 / 77 bits**; *B. californica* PB4207 MATa2 **81 / 82**; both rank above every background protein. The *P. kudriavzevii* control scores 167 (in training), **43.5 with Pichiaceae held out** (167 with Saccharomycetaceae held out; the candidates score 76 / 76 and 78 / 81 with Saccharomycetaceae / Pichiaceae held out, `hmm/a2_no_*`), so the profile finds a Pichiaceae a2 it has not seen.
  - Pfam HMG_box is below GA for all three (13.5, 14.5, control 18.4; GA 25.4), the same weakness as for *Komagataella*; the profile is the better evidence.
  - Synteny: MATa2 and MATa1 are divergent neighbours, 5'-5', 574 bp apart in *P. fabianii* and 625 bp in *B. californica* (275 bp in *P. kudriavzevii*).
  - Conservation: the *B. californica* a2 exon-1 segment scores 41.2 in all four assemblies; *B. discipulorum* has no a2 above 15 by a stop-to-stop scan (introns not modelled), as the earlier agent found.
- **S. schoenii: the HMG-like ORF in the left cluster is not MATa2; it is MATalpha2** (alpha2 profile 81 / 73 bits, a1 9-12, a2 negative; weak Pfam Homeodomain 16 bits). **The real MATa2 is in the right cassette and was missed** because the StringTie transcripts there are thin: whole-genome six-frame scan (`ss_a2_scan.py`, `data/ss_genome_scan.tsv`) finds only **three** a2-like ORFs in the genome, all on chromosome IV, each 800 nt open, a2 profile 77.6 / Saccharomycotina profile 66:
  - **CM143129.1:527,654-528,454 (minus)** and its mirror **537,367-538,167 (plus)**, the outer ends of the palindrome, each head-to-head with an a1 homeodomain (529,186-529,746 plus; 536,074-536,634 minus; gaps 732 and 423 bp);
  - **502,064-502,864 (minus)**, with an a1 at 503,643-504,203 (the central cluster at the centromere).
- Each of the six chromosome cassettes (CM143126-31) carries one alpha2 (71-76 bits), one alpha1 (115) and one a1-like homeodomain (59-61); a2 occurs only on chromosome IV.
- Literature naming (According to PubMed): Hanson, Byrne and Wolfe, PNAS 2014 ([10.1073/pnas.1416014111](https://doi.org/10.1073/pnas.1416014111), PMID 25349420) name MATa2 and MATalpha2 as the genes repressed beside a telomere or centromere in *Pichia pastoris* and *Hansenula polymorpha*, supporting a2/alpha2 as the second gene of each idiomorph.

## Record changes recommended (not applied; db/ untouched)
1. **S. schoenii left record**: rename the 470,349-470,951 ORF (200 aa; 227 aa with the 92 bp intron, 97 reads) **MATalpha2** and the homeodomain at 471,889-472,461 **MATa1** (single exon, 190 aa, uninterrupted; the 217 bp junction at 472,228-472,444 is a read-through isoform toward SLA2, not the gene's intron, so the 103 aa two-exon model should go). The recorded "MATalpha2 = p5" is wrong by homology. Re-check the right-cassette MATa1 intron model (536,074-536,634 is open without the 205 bp intron).
2. **S. schoenii right record**: add **MATa2** at 537,367-538,167 (plus) (mirror copy 527,654-528,454); the family declares MATa2 optional, so this is allowed.
3. **B. californica**: add a partial **MATalpha2** (620,080-620,418, plus, 113 aa, alpha2 profile 37 bits, homeodomain residues 63-155 of the profile, conserved in all four assemblies, 26.5 bits in *B. discipulorum*); keep MATa2 (81 bits).
4. **P. fabianii**: MATa2 supported (73-77 bits); no change. Still no MATalpha2 in the alpha allele: the alpha2 profile finds nothing (could be tested with `hmm/alpha2.hmm` on the TBRC 4498 window).
5. Gene names `MATa1` / `MATalpha2` need only the new profile score difference as evidence, stated as supporting.

## Limits
Public a1/alpha2 sequences are mostly Saccharomycetaceae; titles are NCBI's own (one mislabel found). The S. schoenii background is itself full of MAT paralogs (alpha2 background maximum is the p3 paralog on another chromosome). Strand and exon structures of the S. schoenii ORFs are unmodelled in `ss_a2_scan.py` (uninterrupted segments only). No record was edited; no commit was made.
