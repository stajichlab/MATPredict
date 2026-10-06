# C. cinerea B43 record: curation changelog (2026-10-06)

Branch `curate-cinerea-b43` (from `origin/main` 2d1860a). Applies changes 1-7 of the proposals in
`analysis/2026-10-06_cinerea-b-locus-check.md` (branch `cinerea-b-locus-check`) to
`db/Basidiomycota/Agaricales/5346_a43-b43-okayama-7_PR_B43`, record_version 1 -> 2. The curator approved grouping
them in one branch. Changes 8 and 9 are not applied (follow-ups below). The same text is in the record's
`curation.notes`; this file adds the verification table, the preserved previous content and the validation results.
Sign-off (`reviewed_by`) is still pending.

Nothing was removed without a trace: genes 0-5 and 7 keep their protein accession, locus tag and coordinates; the replaced gene 6
is preserved below; the previous `definition_note` is preserved below; all PMIDs stay and one is added.

## Changes

| # | Change | Result in the record |
|---|---|---|
| 1 | Source strain/genome | `strain.differs_from_sequenced: false -> true`. The sequence is AmutBmut pab1-1 (A43mut B43mut #326, GCA_016772295.1), not the Okayama 7 #130 reference (GCF_000182895.1). The B regions differ only by a 3 bp ACG deletion in the AmutBmut rcb2 coding region (observation, no inference about the mutant phenotype). |
| 2 | Wrong-frame gene 6 | Replaced by `phb3.1`, 1,824,918-1,825,136 (+), 72 aa. |
| 3 | Gene 1 name | `pheromone_B44` -> `phb2_B44-like`; `definition_note` B43/B44 sentence withdrawn. |
| 4 | Missing B43 pheromones | Added `phb1.1` (gene 8), `phb2.3` (9), `phb2.1` (10), `phb3.3` (11); segment end 1,826,859 -> 1,827,702. |
| 5 | Receptor names | gene 0 `rcb1`, gene 3 `rcb2`, gene 5 `rcb3`. |
| 6 | Gene 4 | Kept (KAG2006186.1), renamed `ste3_extra_unassigned`, described as an STE3-type receptor not described in the literature, function unknown. |
| 7 | Group structure, counts, sources | One receptor plus one to three pheromones per group (B43 1+3+3); 14 receptors and 29 pheromones are alleles across 13 haplotypes; PMID 9529887 added to `locus_existence` and `boundaries` citations. |

`order_in_locus` was renumbered to genomic order (12 genes). `gene_index` values were not changed; new genes are 8-11.

## Verification against the genome (done independently of the note)

Region JAAGWA010000010.1:1,800,001-1,830,000 fetched from NCBI; each published B43 deposit searched in it (both strands) and
translated.

| Gene | Result | Note's coordinates | Record coordinates | Check |
|---|---|---|---|---|
| phb1.1 (new) | ORF ATG at 1,809,160, 156 bp with stop; protein 49/51 aa identical to AY172109 (strain B3) | 1,809,177-1,809,316 | 1,809,160-1,809,315 (+) | note started 17 bp inside the ORF (lacks MDDLTV) and ended 1 bp after the stop; the 2 differences are D6V, V21M |
| phb2.3 (new) | exact nucleotide match to AY393916 (162 bp) | 1,813,942-1,814,104 | 1,813,942-1,814,103 (+) | note 1 bp long |
| phb2.2 (gene 2) | exact match to AY393915 | | unchanged 1,814,892-1,815,086 | confirmed |
| phb2.1 (new) | exact match to AY393914 (186 bp, minus strand) | 1,816,393-1,816,579 | 1,816,393-1,816,578 (-) | note 1 bp long |
| phb3.1 (replaces gene 6) | frame 0 from 1,824,918 gives the 72-aa precursor; 1 nt differs from AY393917 (V63G); the NCBI 123-aa model is frame 1 and unrelated | 1,824,918-1,825,137 | 1,824,918-1,825,136 (+) | note 1 bp long (219 bp = 72 codons plus stop) |
| phb3.2 (gene 7) | exact match to AY393918 at 1,826,254-1,826,463 | | unchanged (NCBI model 1,826,186-1,826,859) | confirmed; model longer than the gene |
| phb3.3 (new) | exact match to AY393919 (427 bp, 3 exons: 1,827,276-1,827,338, 1,827,419-1,827,469, 1,827,559-1,827,702), 85 aa | 1,827,276-1,827,703 | 1,827,276-1,827,702 (+) | note 1 bp long; the proposed span end 1,827,703 would include one base past the gene |
| rcb1, rcb2, rcb3 | gene 0 vs AY172107: 527/552 aa (95.5%); gene 3 vs AY393905: 466/480 (97.1%); gene 5 vs AY393906: 414/415 (99.8%) | | unchanged | confirmed |
| gene 4 | 40% to rcb2, 32% to rcb1, 31% to rcb3 (global BLOSUM62) | | unchanged | confirmed |
| gene 1 | 42/61 aa (69%) to Phb2.2 B44 (AY393913); at most 26 identities over 53-61 aligned columns to any B43 pheromone | | unchanged | confirmed (note: 42/60, 70%) |
| ACG deletion | MAFFT of NW_003307535.1:1,713,000-1,742,500 vs JAAGWA010000010.1:1,801,000-1,830,500: 0 mismatches, one 3 bp gap at 1,818,588 (AmutBmut lacks ACG; Okayama 7 #130 and AY393905 carry it) | | | confirmed |

Discrepancies with the note: the five coordinate off-by-one differences above, and the 17 bp start error for phb1.1. The note's
phb1.1 identity (about 49/51 aa) is right; the 2-aa differences from the B3 deposit make phb1.1 the B43 (rcb1 allele 3) copy of a
shared group-1 pheromone, which agrees with Riquelme et al. 2005 but the deposit is from strain B3, so the allele identity is inferred.

## Detection side effects (db/Basidiomycota/order.yml, src)

Gene names are used by detection (`detect.search._attribute` drops a hit whose gene name is not in the family roster). The new
names are declared as `aliases` of existing roster genes, so the roster, the required core and scoring are unchanged:
`pheromone_receptor` aliases rcb1, rcb2, rcb3, ste3_extra_unassigned; `pheromone_B44` aliases phb2_B44-like;
`fungal_mating_type_pheromone` aliases phb1.1, phb2.1, phb2.3, phb3.1, phb3.3. The legacy roster label `pheromone_B44` is kept so
the reported gene names do not change; renaming it is a follow-up. `build_gene_class_index` now resolves aliases (so the GenBank
`/gene_class` qualifier is written for aliased names); test added. What the changed reference proteins do to calls is in the
regression section.

## Preserved: replaced gene 6

NCBI model KAG2006188.1, locus tag CC2G_002524, JAAGWA010000010.1:1,824,970-1,825,393 (+), 123 aa, frame 1 of the true Phb3.1:

    MDASVSTPIYPRPLKVKVARFFQSKRPSPPKSLTESLPTSRDEVGEARLGSALSRSADRIPHMAVASSYTTEYSVQRSTDGEPETSLYRFDSLHLRSRRHDTPFVISILFMSALFIVSSIPSL

## Preserved: previous definition_note (record_version 1)

> The B (PR) locus in C. cinerea is a large multi-group pheromone/receptor complex (O'Shea et al. 1998; Halsall et al. 
> 2000; Riquelme et al. 2005): three gene groups per locus, each with two pheromones and one receptor. A 
> literature-derived gene list (from a separate research agent's synthesis) reported pheromone/receptor annotations for 
> this species scattered across FOUR different contigs of the fragmented assembly GCA_016772295.1 
> (JAAGWA010000001/4/8/10/13). Given C. cinerea has many paralogous pheromone/receptor gene clusters genome-wide (the 
> source literature reports ~14 receptors and ~29 pheromones total across the species), only ONE contig cluster 
> (JAAGWA010000010.1, ~1806154-1826859) was used for this record: it is the only one containing a gene explicitly 
> annotated "Pheromone Phb2.2 B43", directly matching this strain's known B43 allele. Genes on the other three contigs 
> were deliberately excluded as likely unrelated paralogous loci elsewhere in the genome, not confirmed as part of this 
> specific B43 specificity module. This cluster also contains a gene annotated "Pheromone Phb2.2 B44" immediately 
> adjacent to the B43 gene. Per user review (J. Stajich, 2026-09-16, domain expert judgment, pending a book-chapter 
> reference to be added later): PR/B loci in basidiomycetes routinely carry duplicated pheromone/receptor gene copies, 
> and finding one-to-many gene copies clustered in adjacent/nearby contig space is sufficient to define a B locus for 
> this database's purposes -- the co-occurrence of both B43 and B44 pheromone genes in one genome is itself consistent 
> with known B-locus complexity in fungi, not a sign of misassignment. This judgment upgraded boundaries evidence below 
> from tier 2 to tier 1. completeness remains "partial" only because the individual gene roster is not claimed to be the 
> full B-locus gene complement, not because the cluster's validity as a B locus is in doubt.