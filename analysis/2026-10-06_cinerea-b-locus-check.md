# C. cinerea B locus: curated record vs literature vs the three BFD genomes (2026-10-06)

Branch `cinerea-b-locus-check` (from `origin/main` 2d1860a). `db/` was read only; the changes below are proposals with evidence.
Question: the curated record `5346_a43-b43-okayama-7_PR_B43` lists 4 receptors, but the literature (Casselton and Olesnicky 1998) describes a ~17 kb locus of three receptor-plus-two-pheromone cassettes. The array study (`origin/b-locus-clustering-assessment`) left this unresolved.

## Summary
- Of the 4 curated receptors, 3 are the literature group 1, 2 and 3 receptors, matched by sequence to published B43 genes: Rcb1 (allele "3", 95% identity), Rcb2 B43 (96-97%) and Rcb3 B43 (99-100%).
- The 4th receptor (gene 4, KAG2006186.1, 397 aa) is a separate STE3/7-TM ORF 336 bp upstream of Rcb2 on the same strand. It is not a split copy of Rcb2 (40% identity to Rcb2). It is not described in any of the papers I could read. It is not a curation artefact of the record's coordinate selection, because it is in the genome. Whether it is a real extra receptor, a pseudogene or a mis-model is not established.
- The pheromone list has real errors. Of 4 listed pheromones, 2 are correct (Phb2.2 B43, Phb3.2 B43). One (gene 6, "fungal mating-type pheromone", 123 aa) is a wrong-frame ORF overlapping the true Phb3.1 B43. One (gene 1, named "pheromone_B44") is a B43-locus gene most similar to the published Phb2.2 B44 (70% over 60 aa); it is not evidence of a B44 haplotype. Four published B43 pheromone genes lie inside the span and are missing from the record (Phb2.1, Phb2.3, Phb3.3 and the group 1 phb1).
- The record sequence is not from the Okayama 7 #130 reference but from AmutBmut (GCA_016772295.1, strain "A43mut B43mut pab1-1 #326", CUHK 2021), yet the metadata says `differs_from_sequenced: false`. Across 29.5 kb the AmutBmut B region and the Okayama 7 #130 reference (GCF_000182895.1) differ at one site: a 3 bp (ACG) deletion in AmutBmut inside the Rcb2 coding region, a one-codon deletion in an exon of the Rcb2 gene.
- BFD holds exactly 3 C. cinerea genomes, and NCBI holds the same 3 assemblies (n = 3 assemblies, 2 independent strains). The two B43-derived genomes have an identical B layout. The third (T48-F, unknown strain/B type) has a different receptor layout and different alleles.

## 1. The curated record (read from `db/`)
`db/Basidiomycota/Agaricales/5346_a43-b43-okayama-7_PR_B43/metadata.yaml`, locus `PR`, idiomorph `B43`, `completeness: partial`, `coordinate_provenance: curator_derived`.
Source genome: GCA_016772295.1, contig JAAGWA010000010.1, 1,806,154-1,826,859 (20,706 bp). The record's `strain` is "A43 B43 / Okayama-7", `differs_from_sequenced: false`. Evidence citations: PMIDs 9539426 (O'Shea 1998), 10757757 (Halsall 2000), 15879506 (Riquelme 2005). `proposed_by: claude-literature-search`, not reviewed.

| idx | name in record | protein | locus tag | coordinates (strand) | aa | NCBI product name |
|---|---|---|---|---|---|---|
| 0 | pheromone_receptor | KAG2006176.1 | CC2G_002514 | 1806154-1808765 (-) | 576 | pheromone receptor |
| 1 | pheromone_B44 | KAG2006179.1 | CC2G_002517 | 1811951-1812676 (-) | 61 | pheromone Phb2.2 B44 |
| 2 | pheromone_B43 | KAG2006182.1 | CC2G_002519 | 1814892-1815086 (+) | 64 | pheromone Phb2.2 B43 |
| 3 | pheromone_receptor | KAG2006185.1 | CC2G_002521 | 1817324-1819437 (-) | 484 | pheromone receptor |
| 4 | pheromone_receptor | KAG2006186.1 | CC2G_002522 | 1819643-1821277 (-) | 397 | pheromone receptor |
| 5 | pheromone_receptor | KAG2006187.1 | CC2G_002523 | 1821954-1823946 (+) | 417 | pheromone receptor |
| 6 | fungal_mating_type_pheromone | KAG2006188.1 | CC2G_002524 | 1824970-1825393 (+) | 123 | fungal mating-type pheromone |
| 7 | fungal_mating_type_pheromone | KAG2006189.1 | CC2G_002525 | 1826186-1826859 (+) | 69 | fungal mating-type pheromone |

Product names are from NCBI protein records (efetch). KAG2006176, 185, 186 and 187 carry a STE3 domain (Pfam 02076) and 7 predicted TM helices in their NCBI annotation. The record's `definition_note` says four contigs (JAAGWA010000001/4/8/10/13) carried pheromone/receptor annotations and only the contig-10 cluster was used because it carries "Phb2.2 B43". The inventory shows contigs 4, 8, 9 and 12 each carry an unrelated single STE3 locus (30-50% identity to the B receptors), consistent with that exclusion.

## 2. Literature (what I retrieved and what each says)
The PubMed MCP connector was not connected in this session. Metadata came from NCBI E-utilities; full text from PMC HTML where PMC serves it.

| Citation | Retrieved | What it says (relevant) |
|---|---|---|
| Casselton LA, Olesnicky NS. Mol Microbiol Mol Biol Rev 1998;62:55-70 (PMID 9529887, PMC98906) | full text | C. cinereus B locus spans ~17 kb of locus-specific DNA with three subfamilies of functionally independent genes; "each allele within a subfamily consists of a 'cassette' of one receptor and two pheromone genes". B6 pheromones are Phb1.1, 1.2, 2.1, 2.2, 3.1, 3.2. Estimated 79 B specificities, 4-5 alleles of each subfamily would suffice. Contrast: S. commune Balpha/Bbeta are separate loci up to 3.5 map units apart. |
| O'Shea SF, Chaure PT, Halsall JR, ... Casselton LA. Genetics 1998;148:1081-90 (PMID 9539426, PMC1460031) | abstract only (PMC serves the PDF) | B6 locus: nine genes, three 7-TM receptors and six pheromone precursors in 17 kb; three functionally independent subfamilies of two pheromone genes plus one receptor; B6 and B42 share alleles of one subfamily; ~79 B specificities estimated. |
| Halsall JR, Milner MJ, Casselton LA. Genetics 2000;154:1115-23 (PMID 10757757, PMC1460978) | abstract only | B42 locus likewise nine genes in three groups of one receptor plus two pheromones; different B loci may share alleles of one or two groups; B42 carries an extra flanking gene mfs1 (major facilitator), which in other B loci lies in a shared flanking region. |
| Riquelme M, Challen MP, Casselton LA, Brown AJ. Genetics 2005;170:1105-19 (PMID 15879506, PMC1451185) | full text (text; Figure 1 image not read) | 13 B specificities analysed. Three groups of genes per B locus; 2 alleles of group 1, 5 of group 2, 7 of group 3; 14 receptors and 29 pheromones across alleles; 70 possible specificities. Nomenclature: receptors rcb1/rcb2/rcb3 = group 1/2/3, pheromones phb1/phb2/phb3 + number; "each group comprises a receptor gene and one to three pheromone genes". B43 (strain OK130 = A43B43): group 1 = the rcb1^3 allele (PCR/sequence from B3, B7, B40, B43, B45), new group 2 allele rcb2^43 and group 3 allele rcb3^43 sequenced from cosmids; B43 and B5 are homoallelic for group 2. GenBank: rcb2 B43 AY393905, rcb3 B43 AY393906, phb2.1/2.2/2.3 B43 AY393914-16, phb3.1/3.2/3.3 B43 AY393917-19. Group 3 receptors of 6 alleles are 62-81% identical, but the B3 group 3 receptor is only ~20% identical to them (68% to Rcb2 B43). Group 1 alleles are 18% identical. phb3.2 B47 is identical to phb3.3 B43. No receptor beyond rcb1/2/3 per locus is mentioned in the text I read. |
| Brown AJ, Casselton LA. Trends Genet 2001;17:393-400 (PMID 11418220) | abstract only | General review of mushroom mating-type genes; the abstract gives no gene counts. Nothing citable for the B locus structure. |
| Stajich JE et al. PNAS 2010;107:11889-94 (PMID 20547848, PMC2900686) | full text | Okayama 7 #130 genome, haploid, 10x WGS, 13 chromosomes; text discusses A and B only in general terms (B-regulated genes, STE3 receptors up-regulated). No B locus gene catalogue in the text. |
| Srivilai P et al. Pak J Biol Sci 2009;12:110-8 (PMID 19579930) | abstract | AmutBmut is a self-compatible homokaryon "whose nuclei carry mutations in both the A and B loci". The abstract does not say what the mutations are. |

Not retrieved: full text of O'Shea 1998 and Halsall 2000, the Swamy et al. AmutBmut derivation paper, and any paper describing receptors outside rcb1/2/3 near the B locus. "Not mentioned in the text I read" is therefore not evidence of absence.

## 3. Reconciling the record with the literature (genome evidence)
Method: published B43 sequences (GenBank AY393905, AY393906, AY172107, AY393914-19) were mapped onto the AmutBmut contig with 18-20-mer anchors (`scratchpad`, not committed), and proteins were compared by local BLOSUM62 alignment. Receptor identities below use the record's own proteins.

| Literature gene (B43 haplotype) | Maps to (JAAGWA010000010.1) | Record gene | Identity |
|---|---|---|---|
| Rcb1 (allele 3, group 1; AY172107, partial) | 1806650-1808634 (-) | gene 0 (1806154-1808765) | 95% over 552 aa |
| phb1.1 (group 1; AY172109, allele from B3) | 1809177-1809316 (+) | not in record | about 49/51 aa in 6-frame scan |
| Phb2.3 B43 (AY393916) | 1813942-1814104 (+) | not in record | 53/53 aa |
| Phb2.2 B43 (AY393915) | 1814892-1815087 (+) | gene 2 (exact) | 64/64 aa |
| Phb2.1 B43 (AY393914) | 1816393-1816579 (-) | not in record | 61/61 aa |
| Rcb2 B43 (AY393905) | 1817624-1819307 (-) | gene 3 (1817324-1819437) | 96% over 480 aa |
| (no match in the papers I read) | 1819643-1821277 (-) | gene 4, 397 aa | 40% to Rcb2 B43, 32% to Rcb1, 31% to Rcb3 |
| Rcb3 B43 (AY393906) | 1822313-1823827 (+) | gene 5 (1821954-1823946) | 99% over 415 aa |
| Phb3.1 B43 (AY393917) | 1824918-1825137 (+) | gene 6 is a different reading frame (1824970-1825393, 123 aa); the correct frame is the 72-aa precursor ending ...CTIA | 71/72 aa to the 6-frame hit |
| Phb3.2 B43 (AY393918) | 1826254-1826464 (+) | gene 7 | 69/69 aa |
| Phb3.3 B43 (AY393919) | 1827276-1827703 (+) | not in record | 47/47 aa (intron-split) |
| Phb2.2 B44-like (gene 1; 1811951-1812676 -) | not a published B43 gene | gene 1 | 42/60 (70%) to Phb2.2 B44, 21/58 to Phb2.1 B44 |

Reading:
1. **The cassette structure is present and matches the literature for B43.** Group 1 (rcb1 + phb1), group 2 (rcb2 + phb2.1, 2.2, 2.3), group 3 (rcb3 + phb3.1, 3.2, 3.3). B43 therefore has up to 3 pheromones per group, which Riquelme allows ("one to three"); "two per group" is the B6/B42 pattern.
2. **Span.** Rcb1 start (1,806,154) to Phb3.3 end (1,827,703) is 21.5 kb. The record spans 20.7 kb because it omits Phb3.3, and the article's ~17 kb is for B6/B42, whose groups carry fewer pheromones. The span disagreement with "17 kb" is explained by allele content, not by error.
3. **Gene 4 (4th receptor).** Not explained by any of the cassettes. It is a separate, distinct ORF with its own STE3 domain and 7 TM helices per the NCBI annotation. The PubMed-derived text does not mention it. The miniprot model for Rcb2 in the array inventory (521-522 aa, locus 1,817,626-1,821,179) extends 40 aa upstream into gene 4's 5' region; it is a model artefact, because the published gene (AY393905) lies wholly inside gene 3.
4. **Gene 1 naming.** NCBI names it Phb2.2 B44 because it is most similar to that deposit; B43 already carries its own Phb2.1/2.2/2.3 B43 in the same group. The `definition_note` inference that a B43 and B44 pheromone co-occurring is itself consistent with B-locus complexity is therefore weaker than stated: it is one more B43-locus pheromone-like gene, not a B44 allele. It is not shown that this gene is functional.
5. **Gene 6** is a mis-framed ORF (NCBI model CC2G_002524, 123 aa). Frame 0 at 1,824,918 encodes MSDLFASLDLFLSSTEDNG...SWFCTIA, a 72-aa precursor with a CAAX end, 71/72 identical to Phb3.1 B43. The 123-aa ORF is in frame 1 and shares no sequence with it.

Verdict on the question: the 4-receptor listing is not a curation error in the sense of the wrong locus or wrong haplotype. The three expected receptors are right. The fourth receptor is a real gene in the B region that is not described in the literature I could read. The pheromone side has real errors (gene 6, gene 1 label, 4 missing genes), and the span and strain metadata need fixing.

## 4. Strain and genome identity (the record vs the sequenced genome)
- NCBI (esummary/BioSample): GCA_016772295.1 = "Coprinopsis cinerea AmutBmut pab1-1", strain "A43mut B43mut pab1-1 #326", Lab derived, CUHK (Hong Kong), 31 contigs, N50 2.82 Mb (BioSample title "Genome resequencing of C. cinerea #326"). GCF_000182895.1 = "okayama7#130", Broad Institute, 2010, 68 contigs, N50 3.47 Mb (13 chromosomes in the RefSeq report). GCA_982397435.1 = "T48-F", Indian Biological Data Centre, 2026, 605 contigs (458 in the BFD copy), N50 0.30 Mb; strain not given in the BioSample record I could fetch (SAMEA118140661 gave an error on efetch).
- These are the only 3 assemblies NCBI returns for "Coprinopsis cinerea"/"Coprinus cinereus" in the assembly database (3 of 3 are in BFD v0.6.0, `samples.parquet`).
- Okayama 7 #130 vs AmutBmut, B region (Okayama 7 #130 NW_003307535.1:1,713,000-1,742,500 vs AmutBmut JAAGWA010000010.1:1,801,000-1,830,500): global alignment of 29.5 kb gives zero mismatches and a single indel, a 3 bp insertion (ACG) in Okayama 7 #130 at AmutBmut 1,818,588. That site lies in the Rcb2 gene (published B43 gene AY393905 maps to 1,817,624-1,819,307, minus strand) and the published B43 sequence carries the ACG, so AmutBmut is the one with the deletion. The CDS model lengths agree (521 vs 522 aa). A one-codon deletion in Rcb2 might be connected to the "Bmut" self-compatible phenotype, but nothing I retrieved states the molecular basis of Bmut, and this is a hypothesis only.
- Outside the B region the two assemblies differ by a 7 kb insertion in AmutBmut at about 1,854-1,861 kb (k-mer anchors; between the B locus and the unrelated 450-aa STE3 paralog), which explains why the B array is 77.4 kb in AmutBmut and 70.4 kb in Okayama 7 #130.

## 5. Per-genome comparison (pipeline v0.6.0 calls, array inventory `ste3_loci_table.tsv.gz`, n = 3 genomes)
| Genome (strain) | STE3 loci (arrays) | B array (loci) | B array span | Loci in array | Pipeline PR call | PR call status | Group assignment (by identity to published Rcb) |
|---|---|---|---|---|---|---|---|
| GCA_016772295.1 AmutBmut #326 | 8 (5) | JAAGWA010000010.1 (4) | 1,801,307-1,878,728 (77.4 kb) | 3 B receptors at 1801307-1808657 (628 aa), 1817626-1821179 (521 aa), 1822312-1829569 (417 aa); plus a 450-aa partial paralog at 1873437-1878728 (65% identity to the T48-F query) | PR, 5 genes, 1,806,650-1,826,463 (19.8 kb), high confidence | unverified-free: `in_pipeline_call`; idiomorph undetermined (not typed as B43) | Rcb1 allele 3 (95%), Rcb2 B43 (99%), Rcb3 B43 (99-100%) |
| GCF_000182895.1 Okayama 7 #130 | 8 (5) | NW_003307535.1 (4) | 1,713,771-1,784,143 (70.4 kb) | 3 B receptors at 1713771-1721121, 1730090-1733646, 1734779-1742036; same 450-aa paralog 36.9 kb away | PR, 5 genes, 1,719,114-1,738,930 (19.8 kb), high | `in_pipeline_call`; idiomorph undetermined | identical to AmutBmut (B region identical bar the 3 bp) |
| GCA_982397435.1 T48-F | 6 (5) | CEVXIV010000011.1 (2) | 74,889-108,669 (33.8 kb) | 537 aa partial (no start) at 74889-81847, 250 aa partial (no stop) at 99314-108669 | PR, 2 genes, 80,060-107,633 (27.6 kb), medium | unverified: admitted only through a strict-CAAX precursor | 537 aa: 56% to Rcb2 B43 (also 57% to Rcb3 B3: the Ia subfamily); 250 aa: 63% to Rcb3 B43 (group 3 type, IIb cluster). No group 1 (Rcb1-type) receptor found in the array |

Other STE3 loci in each genome (outside the B array) are single-locus arrays with 30-50% identity to the B receptors and are not called: AmutBmut contigs 4 (787 aa, partial, in a withheld cluster), 8 (328 aa), 9 (651 aa), 12 (523 aa); Okayama 7 #130 has the same four; T48-F has contigs 5 (345 aa partial), 12 (467 aa complete), 52 (649 aa partial), 170 (633 aa complete).

Variation among genomes (n = 3):
- AmutBmut and Okayama 7 #130 carry the same B43 haplotype (same genes, same order, same alleles; one codon indel). These two are not independent; they share a strain lineage.
- T48-F differs: 2 receptor loci in the B array instead of 3, a 17.5 kb gap between them (4.7 kb between Rcb2 and Rcb3 in B43, with gene 4 in between), and no Rcb1-type receptor detected. Its group 2 and group 3 receptors are 56% and 63% identical to the B43 alleles, so T48-F is not B43 at groups 2 and 3 (an allele match would be about 100%). It is therefore a different B haplotype. The identity numbers are consistent with the allele diversity Riquelme reports (group 3 alleles 62-81% identical, group 1 alleles 18%).

What cannot be concluded:
- Whether gene 4 is a functional receptor, a pseudogene or an annotation artefact (no expression, no transformation test retrieved).
- Whether T48-F truly lacks a group 1 receptor: both array models are partial (no_start / no_stop), the assembly is fragmented (contig 11 is 637 kb, so the loci are not at a contig end), and group 1 alleles are only ~18% identical to each other, so a divergent Rcb1 would be missed by the identity search.
- Which B specificity T48-F has, or its relationship to the literature alleles: only pairwise protein identity to the published Rcb set was done, with partial models. No pheromone-level analysis was done for T48-F.
- Anything about strain-level variation beyond n = 3 assemblies, 2 independent lineages. Other C. cinerea isolates exist in labs (the literature lists 13 B specificities) but have no assemblies in NCBI or BFD.
- The molecular basis of Amut/Bmut.
- Group 1 pheromone count for B43 (only the B3-allele phb1.1 was mapped; I did not look for more).

## 6. Proposed curation changes (not applied; for the curator)
| # | Change | Evidence | Confidence |
|---|---|---|---|
| 1 | Source genome and strain: either keep GCA_016772295.1 and set `differs_from_sequenced: true` with a note that it is AmutBmut (A43mut B43mut pab1-1 #326), or rebuild from GCF_000182895.1 (Okayama 7 #130, the actual A43 B43 reference). The two B regions are identical except the 3 bp. | Sections 1, 4 | High |
| 2 | Fix gene 6: replace the 123-aa wrong-frame ORF with Phb3.1 B43 (JAAGWA010000010.1:1,824,918-1,825,137 + , 72 aa, AY393917); or drop it and note the NCBI model is wrong | Section 3, point 5 | High |
| 3 | Rename gene 1 from `pheromone_B44` to an unambiguous label (for example "phb2.x B44-like"), and revise the `definition_note` sentence that B43/B44 co-occurrence is consistent with B-locus complexity | Section 3, point 4 | High |
| 4 | Add the missing published B43 genes: Phb1.1 (1,809,177-1,809,316 +), Phb2.3 (1,813,942-1,814,104 +), Phb2.1 (1,816,393-1,816,579 -), Phb3.3 (1,827,276-1,827,703 +) and extend the span to 1,827,703 (21.5 kb) | Section 3 table | High for genes 2.1/2.3/3.3 (exact), medium for phb1.1 (about 49/51 aa) |
| 5 | Name the receptors by group: gene 0 = rcb1 (group 1), gene 3 = rcb2 (group 2), gene 5 = rcb3 (group 3), with the Riquelme 2005 allele names | Section 3 table | High |
| 6 | Gene 4: keep it but label it as an extra STE3-type receptor not described in the literature (not group-assigned), or move it to a flank/`non-core` role until it is supported; do not describe the locus as "4 receptors, 3 groups" without the caveat. Do not merge it into Rcb2 | Section 3, point 3 | Medium |
| 7 | Replace "three gene groups, each with two pheromones" with "three groups, each with one receptor and one to three pheromones (B43: 1 + 3 + 3)" and add Casselton and Olesnicky 1998 (PMID 9529887) as a source; the note that "~14 receptors and ~29 pheromones" occur genome-wide is a misreading: Riquelme counts them across alleles of 13 B haplotypes | Section 2 | High |
| 8 | Do not add T48-F as a B43 record. If a record is wanted, make it a separate record with `idiomorph` unassigned (or by allele after mapping) and `status: provisional`; both array models are partial | Section 5 | Medium |
| 9 | `elements` for the array study: the 17 kb locus is smaller than its array because the array rule extends 45 kb to the unrelated 450-aa paralog; no curation change, only a note to report array span separately from B-locus span | Section 5 | High |

Files: this note only. No change to `db/`, `src/` or `results/`.
