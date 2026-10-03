# Receptor queue, step 1: STE3-like receptors in Basidiomycota proteomes

2026-09-27. Read-only exploration. Scripts: proteome_map.py, select_sample.py,
collect.py, analyze.py, draw_tree.py. Tree PDF: receptor_tree.pdf.

## Proteome coverage

- BFD funannotate proteomes: genome_annotation/<Species_Strain>/predict_results/<Species_Strain>.proteins.fa
  (directory name = SPECIES_IN + STRAIN, spaces to "_").
- 2,941 of 3,271 non-suppressed Basidiomycota genomes (90%) have one
  (proteomes.tsv). Lowest: Pucciniales 59/129.
- tables/pfam.parquet covers only 745 genomes (34 Basidiomycota), so it cannot
  replace a search.

## Search and tree

- 238 proteomes sampled across ~45 orders (sample.tsv). hmmsearch Pfam PF02076
  (Pfam 38.2, --cut_ga): 1,060 STE3-like proteins.
- 34 curated mating receptors from db records (branch curation-puccinio):
  Coprinopsis B43 (4), Schizophyllum bar3/bbr2, Ustilaginales pra1 (11),
  Cryptococcus STE3 a/alpha, Sporidiobolales STE3a1 (6) / STE3a2 (7),
  Wallemia STE3 v1/v2. Outgroup: S. cerevisiae STE3 (P06783).
- hmmalign to PF02076 (296 match columns); 983 sequences kept (>=40% of the
  model). FastTree -lg -gamma. Local support values, not bootstrap.

## Copies per genome (median, sampled genomes)

Agaricales 7, Polyporales 6, Cantharellales 6, Russulales 5, Boletales 3,
Hymenochaetales 3, Pucciniales 3, Tremellales 2, Sporidiobolales 1,
Ustilaginales 1, Sebacinales 48, Thelephorales 10, Atheliales 13.

## Do mating receptors separate?

- Sporidiobolales: yes. STE3a1 records form one clade (support 0.94) with only
  Sporidiobolales queries; STE3a2 records form another (0.88). One copy per
  genome, falling in one of the two clades.
- Ustilaginales/Exobasidiomycetes pra1: three clades (support 0.91-0.996),
  each with queries from Ustilaginales, Microstromatales, Exobasidiales and
  others; not one clean clade.
- Agaricomycetes: no. The Coprinopsis/Schizophyllum receptors sit in several
  small Agaricales-only subclades (12-31 tips, support 0.57-0.99), interleaved
  with non-mating copies; the known set is not monophyletic. 3-7 copies per
  genome (48 in Sebacinales), so tree placement alone cannot pick the mating
  receptors.
- Tremellales: Cryptococcus STE3alpha sits in a small Tremellomycetes clade
  (0.95); STE3a sits on a long branch deep in the tree.
- Wallemia: v1 STE3 groups with Tremellomycetes copies; v2 sits on a long
  branch.
- Pucciniales: 0 of 10 genomes have a copy in any known mating clade (no rust
  receptor is curated).

## Limits

- One domain alignment; FastTree local supports only.
- Sample of 238 of 2,941 proteomes.
- No positional evidence used yet (pheromone precursor with CAAX, B-locus
  neighbourhood).
