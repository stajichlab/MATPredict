# Reference-gap triage: which orders lack references, and how many BFD calls may be wrong or missed

Status: measurement only (2026-10-05/06); no code or database change. Data and scripts: `results/2026-10-05_reference_gap_triage/` and `results/2026-10-05_botryosphaeriales_deposits/`. Calls are from the v0.6.0 campaigns (`results/2026-10-03_ascomycota_v060`, `results/2026-10-03_basidiomycota_v060`).

## Question
Which orders have no or poor reference members, how many existing BFD results could be wrong or missed, and would additional deposited sequences help?

## Method
1. `risk_screen.py`: join campaign `genomes.tsv` and `loci.tsv` to BFD assembly quality (`busco_genome.parquet`, `asm_stats.parquet`). Good assembly means BUSCO >= 70, N50 >= 20 kb and contigs <= 5000. A not-called genome on a good assembly is a candidate miss.
2. `risk_excess.py`: excess misses = not called on a good assembly beyond a 3% background rate per order (an assumption, not a measured miss rate).
3. Idiomorph skew screen: species with 10 or more single-idiomorph genomes and a minority share of 15% or less (binomial p < 1e-4 against 1:1).
4. `rescue_test.py`: deposited MAT proteins (`rescue_queries.faa`) searched against uncalled good-assembly genomes in four orders (core hit >= 35% identity over >= 50% of the query, plus a flank within 30 kb), with 30 called genomes per order as controls.
5. `botry_deposit_small.py`: *Diplodia sapinea* and other deposited Botryosphaeriales loci against 12 *B. dothidea* genomes and 10 controls.

## Results
### Calls and misses
- Ascomycota: 19,415 genomes; 18,166 good assemblies, 1,651 of them not called (9.1%). Basidiomycota: 3,270 genomes; 2,621 good assemblies, 217 not called (8.3%).
- Excess misses (beyond a 3% background): Ascomycota about 1,278 genomes (6.6%), nine orders carry 80% of them: Dipodascales 214, Saccharomycetales 184 (Saccharomycetaceae), Xylariales 170, Pichiales 162, Pleosporales 114 (Pleosporaceae), Orbiliales 58, Chaetothyriales 57, Amphisphaeriales 41, Microascales 37, Mycosphaerellales 37. Basidiomycota about 182 genomes (5.6%), six orders carry 80%: Trichosporonales 60, Microbotryales 29, Filobasidiales 21, Cystofilobasidiales 19, Cystobasidiales 11, Kriegeriales 8.
- Fallback-routed orders (hundreds of noise clusters per genome): Dipodascales 87%, Pichiales 99%, Orbiliales 100%, Ascoideales 100%, Trichosporonales, Microbotryales, Filobasidiales, Cystofilobasidiales 100%.
- Xylariales (71% of good assemblies not called) is a biology question: no canonical MAT locus is known in the order (`docs/superpowers/specs/2026-10-05-xylariales-models-design.md`), so many of these may be true absences of a canonical locus, not misses.
- Weak-only genomes (every call low confidence, partial or gene only): 3,171 Ascomycota and 895 Basidiomycota; called through the phylum fallback: 811 and 456.
- Idiomorph skew: 29 of 133 species screened (1,546 genomes). This screen is noisy because homothallic and clonal species legitimately show one idiomorph; do not read it as an error rate.

### Rescue test with deposited sequences
| Order | Uncalled (good) | Rescue candidates | Median identity | Controls concordant |
|---|---|---|---|---|
| Microascales | 42 | 11 (Ceratocystidaceae; KF033902/3) | 92% | 20 of 20 core hits |
| Pichiales | 178 | 9 | 43% | 12 of 12 core hits |
| Chaetothyriales | 63 | 2 | 46% | 17 of 18 core hits |
| Dipodascales | 224 | 63 | 37% | 0 of 5 core hits |

- Microascales/*Ceratocystis* is the clear quick win: high identity and every control concordant.
- Pichiales and Chaetothyriales gain little from the deposits tested (9 and 2 candidates).
- Dipodascales candidates are probably paralog noise at about 37% identity with 0 of 5 control concordance; do not curate from them.

### Botryosphaeriales
- 12 *B. dothidea* genomes: current references hit the locus at a median 75.6% identity; the *D. sapinea* deposits reach 41.3%. The 44.5% case is one 3 kb fragment on a small contig. Ten control genomes: current references median 83.7%. Deposits help only for the close relative *D. seriata* (89%). Botryosphaeriales detection is not reference-limited at order level; curating *D. sapinea* is low priority.

## What can be said
- Roughly 1,460 genomes (about 6.4% of 22,685 Ascomycota and Basidiomycota genomes) are excess misses on good assemblies, concentrated in 15 orders. This is an upper estimate of detection gaps, because absence of a canonical locus (Xylariales) or genuine loss is not separated from a miss.
- Wrong calls (as opposed to missed ones) are not measured here. The 1,360 unverified PR calls and the weak-only genomes are the main candidates (`analysis/2026-10-06_agaricomycetes-pr-arrays.md`).

## Priority list
1. Quick wins by curated records: Microascales (*Ceratocystis*; deposits exist), then Pichiales and Chaetothyriales (find better deposits first).
2. Routing and gate: Trichosporonales, Filobasidiales, Cystofilobasidiales, Microbotryales (see `analysis/2026-10-05_tremellomycetes-failures.md` on branch `tremello-failure-investigation`); needs a before/after regression.
3. Investigate despite existing records: Saccharomycetaceae (248 not called) and Pleosporaceae (114 not called).
4. Needs a lineage study before records: Xylariales, Orbiliales, Amphisphaeriales, Dipodascales outside *Yarrowia*.

## Limits
- The 3% background is an assumption. Quality filters rely on BUSCO and N50 only. Rescue tests use one deposit set per order and one identity threshold. Counts are from v0.6.0 code and predate the polish identity tier (PR #40), which recovered some Dothideomycete calls.
