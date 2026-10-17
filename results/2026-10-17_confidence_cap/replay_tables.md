21011 loci, 13740 high

## Candidate rules on the v0.6.0 `high` calls

| rule | drops | % of high | mat_locus | gene_only | homothallic |
|---|---|---|---|---|---|
| P50: best-core coverage <50 | 5615 | 40.9 | 5482 | 133 | 0 |
| P60: best-core coverage <60 | 5878 | 42.8 | 5707 | 153 | 18 |
| P70: best-core coverage <70 | 6546 | 47.6 | 6343 | 180 | 23 |
| P80: best-core coverage <80 | 6791 | 49.4 | 6550 | 216 | 25 |
| P90: best-core coverage <90 | 7121 | 51.8 | 6850 | 242 | 29 |
| P60+locus-edge: best-core coverage <60 and LOCUS segment within 500 bp of a contig end | 304 | 2.2 | 235 | 69 | 0 |
| P80+locus-edge: best-core coverage <80 and LOCUS segment within 500 bp of a contig end | 378 | 2.8 | 275 | 102 | 1 |
| W40: best-core identity <40, no flank gene | 54 | 0.4 | 0 | 54 | 0 |
| W40max: highest core identity <40, no flank gene | 27 | 0.2 | 0 | 27 | 0 |
| W45: best-core identity <45, no flank gene | 119 | 0.9 | 0 | 119 | 0 |
| W45max: highest core identity <45, no flank gene | 94 | 0.7 | 0 | 94 | 0 |
| W50: best-core identity <50, no flank gene | 162 | 1.2 | 0 | 162 | 0 |
| W50max: highest core identity <50, no flank gene | 147 | 1.1 | 0 | 147 | 0 |
| W45+edge: best-core identity <45, no flank, within 500 bp of contig end | 62 | 0.5 | 0 | 62 | 0 |
| W45+gene_only: best-core identity <45, no flank, locus_class idiomorph_gene_only | 119 | 0.9 | 0 | 119 | 0 |
| PE: partial_edge = best-core cov <80, no core gene >=90, best gene within 500 bp of a contig end | 123 | 0.9 | 63 | 59 | 1 |
| WK: weak_identity = best-core identity <45, no flank gene | 119 | 0.9 | 0 | 119 | 0 |

## Truth set (BioSample MAT1-1/MAT1-2 stated, Ascomycota:MAT/MATtub high loci)

122 high loci in truth genomes: right 103 (Botrytis 14), wrong 19 (Botrytis 16)

| rule | demoted right | demoted wrong (non-Botrytis) | demoted wrong (Botrytis, fixed by PR #95) |
|---|---|---|---|
| P50: best-core coverage <50 | 47 | 3 | 0 |
| P60: best-core coverage <60 | 49 | 3 | 0 |
| P70: best-core coverage <70 | 62 | 3 | 16 |
| P80: best-core coverage <80 | 74 | 3 | 16 |
| P90: best-core coverage <90 | 74 | 3 | 16 |
| P60+locus-edge: best-core coverage <60 and LOCUS segment within 500 bp of a contig end | 4 | 0 | 0 |
| P80+locus-edge: best-core coverage <80 and LOCUS segment within 500 bp of a contig end | 4 | 0 | 0 |
| W40: best-core identity <40, no flank gene | 0 | 0 | 0 |
| W40max: highest core identity <40, no flank gene | 0 | 0 | 0 |
| W45: best-core identity <45, no flank gene | 0 | 0 | 0 |
| W45max: highest core identity <45, no flank gene | 0 | 0 | 0 |
| W50: best-core identity <50, no flank gene | 0 | 0 | 0 |
| W50max: highest core identity <50, no flank gene | 0 | 0 | 0 |
| W45+edge: best-core identity <45, no flank, within 500 bp of contig end | 0 | 0 | 0 |
| W45+gene_only: best-core identity <45, no flank, locus_class idiomorph_gene_only | 0 | 0 | 0 |
| PE: partial_edge = best-core cov <80, no core gene >=90, best gene within 500 bp of a contig end | 0 | 0 | 0 |
| WK: weak_identity = best-core identity <45, no flank gene | 0 | 0 | 0 |

## Sensitivity of the weak-identity threshold (high, no flank gene)

| identity below | loci |
|---|---|
| 35 | 32 |
| 40 | 54 |
| 45 | 119 |
| 50 | 162 |
| 55 | 200 |
| 60 | 227 |
| 65 | 249 |
| 70 | 255 |
| 80 | 273 |

Histogram of best-core identity, high gene-only loci (5-point bins):

25:5, 30:27, 35:22, 40:65, 45:43, 50:38, 55:27, 60:22, 65:6, 70:10, 75:8, 80:37, 85:10, 90:26, 95:552, 100:1041

Sensitivity of the edge distance (partial_edge: best cov <80, no full core gene, gene-level distance):

| edge below (bp) | loci |
|---|---|
| 100 | 57 |
| 500 | 123 |
| 1000 | 170 |
| 2500 | 326 |
| 5000 | 744 |
| 10000 | 1075 |

## Implemented rule (partial_weak = PE or WK): 208 of 13740 high (1.51%); PE 123, WK 119, both 34

**all** (208)

- family: MTL 113, MAT 58, MATsc 37
- locus_class: idiomorph_gene_only 144, mat_locus 63, homothallic_candidate 1
- idiomorph: alpha 75, A 38, MATalpha 35, MAT1-1 29, MAT1-2 29, MATa 2
- class: Pichiomycetes 102, Saccharomycetes 45, Dothideomycetes 26, Eurotiomycetes 13, Sordariomycetes 12, Lecanoromycetes 4, Leotiomycetes 2, Sporopachydermiomycetes 2, Lichinomycetes 1, NA 1
- order: Pichiales 89, Saccharomycetales 25, Pleosporales 15, Ascoideales 13, Serinales 13, Eurotiales 12, Hypocreales 9, Phaffomycetales 7, Dothideales 4, Mycosphaerellales 3, Lecanorales 3, Helotiales 2, Sporopachydermiales 2, Botryosphaeriales 2, Magnaporthales 2
- genome with another call of the same family at `high`: 30

**partial_edge** (123)

- family: MAT 58, MTL 40, MATsc 25
- locus_class: mat_locus 63, idiomorph_gene_only 59, homothallic_candidate 1
- idiomorph: alpha 32, MAT1-1 29, MAT1-2 29, MATalpha 23, A 8, MATa 2
- class: Pichiomycetes 42, Dothideomycetes 26, Saccharomycetes 23, Eurotiomycetes 13, Sordariomycetes 12, Lecanoromycetes 4, Leotiomycetes 2, Lichinomycetes 1
- order: Pichiales 36, Saccharomycetales 18, Pleosporales 15, Eurotiales 12, Hypocreales 9, Serinales 6, Dothideales 4, Ascoideales 3, Mycosphaerellales 3, Lecanorales 3, Helotiales 2, Botryosphaeriales 2, Phaffomycetales 2, Magnaporthales 2, Caliciales 1
- genome with another call of the same family at `high`: 14

**weak_identity** (119)

- family: MTL 98, MATsc 21
- locus_class: idiomorph_gene_only 119
- idiomorph: alpha 64, A 34, MATalpha 20, MATa 1
- class: Pichiomycetes 90, Saccharomycetes 26, Sporopachydermiomycetes 2, NA 1
- order: Pichiales 82, Ascoideales 12, Saccharomycetales 9, Serinales 8, Phaffomycetales 5, Sporopachydermiales 2, (none) 1
- genome with another call of the same family at `high`: 17


## High calls per family (denominator) and share dropped

| family | high | dropped | % |
|---|---|---|---|
| Ascomycota:MAT | 9301 | 58 | 0.6 |
| Ascomycota:MTL | 2255 | 113 | 5.0 |
| Ascomycota:MATsc | 1920 | 37 | 1.9 |
| Ascomycota:mat1 | 96 | 0 | 0.0 |
| Ascomycota:mat2 | 69 | 0 | 0.0 |
| Ascomycota:MATyl | 53 | 0 | 0.0 |
| Ascomycota:MATtub | 42 | 0 | 0.0 |
| Ascomycota:PM | 4 | 0 | 0.0 |

## High calls per order (top, by drops)

| order | high | dropped | % |
|---|---|---|---|
| Pichiales | 228 | 89 | 39.0 |
| Saccharomycetales | 1752 | 25 | 1.4 |
| Pleosporales | 480 | 15 | 3.1 |
| Ascoideales | 25 | 13 | 52.0 |
| Serinales | 2099 | 13 | 0.6 |
| Eurotiales | 1680 | 12 | 0.7 |
| Hypocreales | 3099 | 9 | 0.3 |
| Phaffomycetales | 85 | 7 | 8.2 |
| Dothideales | 191 | 4 | 2.1 |
| Mycosphaerellales | 222 | 3 | 1.4 |
| Lecanorales | 100 | 3 | 3.0 |
| Helotiales | 275 | 2 | 0.7 |
| Sporopachydermiales | 4 | 2 | 50.0 |
| Botryosphaeriales | 714 | 2 | 0.3 |
| Magnaporthales | 443 | 2 | 0.5 |

## locus_class of dropped

| locus_class | high | dropped |
|---|---|---|
| mat_locus | 11574 | 63 |
| idiomorph_gene_only | 1848 | 144 |
| homothallic_candidate | 318 | 1 |

Gene-only `high` loci: 1848; dropped 144 (7.8%); kept 1704.
