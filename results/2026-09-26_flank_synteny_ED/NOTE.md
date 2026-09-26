# Flank-ortholog synteny at the early-diverging sexM/sexP loci (draft)

2026-09-26. Scripts and outputs are in this directory. The data comes from the
early-diverging scan (`results/2026-09-26_early_diverging/`, code `634dda4`).

## Question

Mortierellomycota and Kickxellomycota have sexM/sexP-like HMG hits next to
tptA/rnhA/glrA hits. Are those flanking hits the genome's true orthologs? Does
the HMG gene sit between them, as it does in Mucorales (tptA–HMG–rnhA)?

## Method

1. **Forward search.** tblastn of all 41 curated Mucoromycota MAT-locus
   proteins against each genome. Settings: e ≤ 10, `-seg no`. This matches
   detect, which uses the tblastn default e-value. At e ≤ 1e-5,
   *Coemansia* has no tptA hit at all.
2. **Flank regions.** For each flank gene (tptA, rnhA, glrA, algA, btbA), HSPs
   on the same contig within 3 kb form one region. The region protein is the
   concatenated HSP translations (`regions.py`).
3. **Reverse search.** blastp of each region protein against three Mucorales
   proteomes: Phycomyces L51, Mucor circinelloides 1006PhL and Rhizopus
   arrhizus 97-1192 (BFD funannotate models).
   - **Confirmed:** the region's top hit is that gene's ortholog (see
     `ref/flank_orthologs.tsv`). In Mucor and Phycomyces the tptA/algA/rnhA
     orthologs are neighbouring gene models, and each curated query maps to
     exactly one protein per proteome.
4. **Genome ortholog.** The confirmed region with the best forward bitscore.
   - **Strict:** that region is also the gene's best forward region.
5. **Tests**
   - **Genome level (independent of detect):**
     - `pair`: the tptA and rnhA orthologs are ≤ 100 kb apart on one contig.
     - `pair+HMG`: a sexM/sexP HSP (e ≤ 1e-3) lies between them.
     - `HMG≤20kb`: an HMG HSP lies within 20 kb of any flank ortholog.
   - **Locus level (loci from the detection reports):**
     - Is a confirmed flank within ±5 kb of the locus?
     - Does an HMG HSP lie between confirmed regions of two different flank
       genes?

Suppressed ASMIDs (BFD `data/curation/suppress.txt`) were skipped. 608 genomes
were searched and 602 were analysed. The other 6 are the chytrids that timed out
in detect, so they have no detection report.

## Result: the Mucorales arrangement is absent in both phyla

Genome level, default thresholds (`summary_default.txt`):

| group | genomes | tptA ortholog | rnhA ortholog | pair | pair+HMG | HMG ≤ 20 kb |
|---|---:|---:|---:|---:|---:|---:|
| Mucoromycota (control) | 293 | 280 | 293 | 139 | 117 | 241 |
| Mortierellomycota | 100 | 100 | 100 | **0** | **0** | **0** |
| Kickxellomycota | 190 | 27 | 91 | **0** | **0** | **0** |
| chytrid control | 19 | 17 | 19 | 0 | 0 | 0 |

Locus level:

| group | kind | loci | confirmed flank at locus | HMG between confirmed flanks |
|---|---|---:|---:|---:|
| Mucoromycota | called | 242 | 238 | 167 |
| Mucoromycota | withheld | 393 | 28 | 0 |
| Mortierellomycota | called | 6 | 3 | **0** |
| Mortierellomycota | withheld | 388 | 43 | **0** |
| Kickxellomycota | called | 7 | 1 | **0** |
| Kickxellomycota | withheld | 744 | 19 | **0** |
| chytrid control | withheld | 163 | 2 | 0 |

- **Mortierellomycota:**
  - Both flank orthologs were confirmed in 100 of 100 genomes (99 strict).
  - tptA and rnhA are never within 100 kb of each other.
  - No HMG hit lies within 20 kb of either ortholog.
  - All 6 calls fail the test. Only 3 of them have any confirmed flank at the
    locus.
- **Kickxellomycota:**
  - No locus passes the test.
  - The reverse search confirms tptA in only 27/190 genomes and rnhA in only
    91/190. So the ortholog test is weaker in this phylum.
- **Looser thresholds** (`summary_loose.txt`: pair ≤ 300 kb, HMG ≤ 100 kb,
  HMG e ≤ 1):
  - Still no pair and no locus with HMG between confirmed flanks, in either
    phylum.
  - Some genomes gain an HMG hit within 100 kb of a flank ortholog:
    Mortierellomycota 11/100 and Kickxellomycota 13/190. I have not inspected
    these genomes.
- **The HMG hits are not rare in these phyla.** Every Mortierellomycota and
  Kickxellomycota genome has sexM/sexP HSPs at e ≤ 1e-3 (median 62 and 69 per
  genome). The chytrids have them too (24/25, median 45). This hit count is
  background from the HMG-box family. It does not by itself show that a MAT
  gene is present.

## Control shows the test works, but not in every Mucoromycota family

- **Positive control:** 167/242 called Mucoromycota loci have an HMG hit
  between confirmed flanks. For withheld loci it is 0/393.
- **Families with no tptA–rnhA pair:** Lichtheimiaceae 0/30,
  Umbelopsidaceae 0/14, Syncephalastraceae 0/9 and Endogonaceae 0/4.
  - Umbelopsidaceae is called 13/14 by detect.
  - So the tptA–rnhA pair is not conserved across all of Mucoromycota. Where
    it is conserved, it is a good positive marker.
  - Its absence alone does not prove a locus is absent.
- **Rhizopodaceae:** 45/111 genomes have a pair and 95/111 have an HMG hit
  ≤ 20 kb from a flank. This is consistent with assembly breaks or a
  different arrangement. I have not checked which.

## Reading

- **Mortierellomycota: HMG paralog next to flank paralogs.**
  - At the called loci, the "flank" hits are mostly not the tptA/rnhA
    orthologs.
  - The true orthologs are present and confirmed in every genome, but they are
    not linked to each other or to any HMG hit.
  - The calls do not have the support the scan note asked for.
  - One limit: MAT may sit at a different position in this phylum, as in
    Umbelopsidaceae. This test cannot find a locus in a new place.
- **Kickxellomycota: inconclusive, leaning paralog.**
  - No locus passes the test.
  - But the ortholog step is weak: tptA was confirmed in 14% of genomes.
- **Chytrid control:** 0 loci pass. This agrees with the 0 calls from detect.

## Limits

- The flanks are confirmed against Mucorales proteomes only. Divergent
  orthologs in Kickxellomycota can fail this check.
- tblastn HSPs are not gene models. "Between" means HSP coordinates.
- The HMG hits are not typed: sexM/sexP HSPs include other HMG-box proteins.
- The test assumes the Mucorales gene order. It cannot find a MAT locus that
  has moved to another position.
- Near-clonal *Coemansia* isolates are counted per genome. The per-species
  count (0 of 169 species) gives the same answer.

## For the curator

1. Treat the 6 Mortierellomycota and 7 Kickxellomycota calls as unsupported.
   Options: withhold them, or label them `unverified_lineage`.
2. Should we search for a relocated locus? One way: take the HMG hits in these
   phyla that pass a typing step, for example the planned sexM/sexP phylogeny,
   and look at their neighbouring genes. There is no test for this yet.
3. Lichtheimiaceae, Umbelopsidaceae and Syncephalastraceae lack the
   tptA–rnhA pair. Look at their flank arrangement before they are curated.
