# Parastagonospora MAT1-1 / MAT1-2 calls: deeper check (2026-10-17)

Branch `parastagonospora-check`. No change to `src/` or `db/`. Data and scripts: `results/2026-10-17_parastagonospora_check/`.
Inputs: v0.6.0 reports (`core-hmms/asco`), records-run reports (`core-hmms/dothideo`), 186 genomes with species starting
*Parastagonospora* (175 *P. nodorum*, 11 *P. avenae*). Profiles from branch `ascomycota-core-hmms` (alpha box, HMG box),
run unchanged.

## 1. Are the second-idiomorph hits in the records-run loci real genes?

No. `gene_audit.py` translated every MAT-named gene entry of every Ascomycota:MAT locus (both runs, 1,168 entries) and scored
it on both profiles (`gene_audit.tsv`). "Strong" = modelled, at least 300 aa, own-profile score at least 150 bits.

Records run, 182 loci (one locus per genome, all `mat_locus`, all `high`):

| locus label | n | strong gene | other core gene |
|---|---|---|---|
| MAT1-1 | 91 | MAT1-1-1: modelled, 343-344 aa, alpha 225-271 bits, 71-100% identity, bitscore 434-563 | 88 weak MAT1-2-1, 3 none |
| MAT1-2 | 91 | MAT1-2-1: modelled, 340-345 aa, HMG 224-246 bits, 75-100% identity, bitscore 265-377 | 89 weak MAT1-1-1, 2 none |

- Loci with two strong, modelled, full-length genes (homothallic-like): **0** of 182.
- The weak partners are all unmodelled, 36-113 aa. Weak MAT1-2-1 (n=88): bitscore at most 35.4, HMG score at most 4.8, alpha at most 15.3.
  Weak MAT1-1-1 (n=89): bitscore at most 31.2 (106 aa, 32% identity, from the *F. graminearum* record `5518_3639_MAT_combined`), own alpha score 0,
  HMG-profile score 0-64.5 (a fragment of the real HMG gene). These are cross-hits of the shared HMG/alpha-like weakly conserved regions, not genes.
- The earlier "88 loci list both MAT1-1-1 and MAT1-2-1" is this: 88 MAT1-1 loci carry a weak MAT1-2-1 hit, and 89 MAT1-2 loci carry a weak MAT1-1-1 hit.
- MAT1-1-4 (179 loci, every label): 36 aa, 29-33 bits, alpha and HMG score 0, unmodelled, from `199306_rmscc1040_MAT_MAT1-1`. It is ubiquitous, carries no
  idiomorph information here and should not vote. MAT1-1-3 (91 loci in both runs): the HMG gene, see section 4. MAT1-2-10: absent (only on the *Sclerotinia trifoliorum* record).
- v0.6.0 for contrast: in the 87 wrong loci the MAT1-1-1 entry is an unmodelled 106-aa, 31.2-bit hit (alpha score 0), MAT1-2-1 an unmodelled 30-aa 55.1-bit hit
  (HMG 45.3), and MAT1-1-3 a modelled 76-aa HMG fragment (HMG 114.6, alpha 5.9). No full-length gene was modelled in any of them.

## 2. Independent truth

- Literature. Syme et al. 2013 G3 3:959-969 (doi 10.1534/g3.112.004994, PMID 23589517), "Mating type locus" section: SN15 is MAT1-1, Sn4 and Sn79 are MAT1-2.
  Records-run calls: SN15 (GCA_008452785.1) MAT1-1; LDN03-Sn4 (GCA_002267005.1) MAT1-2; Sn79 (GCA_002267025.1 and GCA_051013865.1) MAT1-2. 4 of 4. v0.6.0 called Sn4, Sn79 and Sn2000
  all MAT1-1 (1 of 4 right, SN15 only). Sn2000 (GCA_002267045.1) is called MAT1-2; I found no published mating type for it.
  Bennett et al. 2003 Fungal Genet Biol 40:25-37 (doi 10.1016/s1087-1845(03)00062-8, PMID 12948511) sequenced both idiomorphs, SN435PL98 MAT1-1 (AY212018.1, MAT1-1-1 AAO31740.1)
  and SN436GA98 MAT1-2 (AY212019.1, MAT1-2-1 AAO31742.1), a pair known to mate; these are the repo's own Pleosporales records and the records-run references.
  Population surveys at about 1:1: Sommerhalder et al. 2006 (doi 10.1094/PHYTO-96-0234, 1,207 isolates, 17 of 18 fields not different from 1:1), Adhikari et al. 2008
  (doi 10.1094/PHYTO-98-1-0101), Cowger and Silva-Rojas 2006 (doi 10.1094/PHYTO-96-0860), Kaur et al. 2024 (doi 10.1094/PHYTO-10-22-0392-R).
  Jones et al. 2024 Commun Biol 7:1614 (doi 10.1038/s42003-024-07327-7), the study behind 152 of our genomes (PRJNA612761): core WA population "close to 1:1 ... heterothallic meiotic potential";
  clonal sub-groups carry a single mating type; one WA and one Iranian isolate matched both loci. McDonald 2012/2013 and the Richards et al. 2019 PLoS Genet paper
  (doi 10.1371/journal.pgen.1008223) were not read for mating-type statements; I did not confirm any claim from them.
- BioSample. None of the 186 is in `results/2026-10-17_biosample_truth/` (0 Parastagonospora rows in `truth_set.tsv`, `biosample_mating_attributes.tsv`, `agreement.tsv`).
  I fetched all 186 BioSamples (`para_bs.tsv`): no mating-type or sex attribute in any (attributes: strain, isolate, geo_loc_name, host, collection_date, and similar).
- The 4 "none" genomes. tBLASTn of the two Bennett proteins (`tblastn_none_genomes.tsv`, `assembly_stats_none_and_controls.txt`):

| genome | species | assembly | tBLASTn | verdict |
|---|---|---|---|---|
| GCA_003497985.1 | nodorum | 296 contigs, N50 317 kb | MAT1-2-1 at 99.4% / 98.9% (two HSPs, 345 aa) on QVWR01000152.1:57395-58482 | gene present: software miss, MAT1-2 |
| GCA_003503065.1 | nodorum | 6,706 contigs, N50 7 kb | MAT1-2-1 at 99.4% / 98.9% on QVWJ01004591.1:318-1405 | gene present, 318 bp from the contig end: software miss on a fragmented assembly, MAT1-2 |
| GCA_022559985.1 | avenae | 408 contigs, N50 317 kb | MAT1-1-1 at 85.6% (271 aa) + 56% (N-terminal) on JAFJTR010000080.1:3562-4648 | gene present: software miss, MAT1-1 |
| GCA_003503165.1 | avenae | 5,714 contigs, N50 5 kb, 25.7 Mb | no hit at e-value 1e-3 | assembly too fragmentary (25.7 Mb against 36-41 Mb for the others); not a call |

  Three of four are misses by the software on genomes that carry a clear MAT gene; the reason is not established (all four reports have `detected: []` in both runs).
  The likeliest cause is the evidence floor (one gene at the locus) but I did not test it.

## 3. Composition and clonality

(`genomes_table.tsv`; labels = records run; clusters = mash distance below 0.001, `para_mash.tsv.gz`, k=21, s=20000)

| group | genomes | MAT1-1 | MAT1-2 | none | exact binomial p |
|---|---|---|---|---|---|
| all | 186 | 91 | 91 | 4 | 1.00 |
| *P. nodorum* | 175 | 86 | 87 | 2 | 1.00 |
| *P. nodorum*, Western Australia (Jones et al. set) | 152 | 77 | 75 | 0 | 0.94 |
| *P. nodorum*, elsewhere (USA 9, Iran 3, Sweden 2, Finland 1, Switzerland 1, Brazil 1, China 1, no location 3, South Africa 1, Australia outside WA 1) | 23 | 9 | 12 | 2 | 0.66 |
| *P. avenae* (f. sp. tritici and avenae) | 11 | 5 | 4 | 2 | 1.00 |
| one genome per clone cluster | 122 | 60 | 62 | - | 0.93 |
| WA, one per cluster | 94 | 46 | 48 | - | 0.92 |

The 186 form 126 clusters: 105 singletons and 21 multi-member clusters holding 81 genomes (largest 14, 12, 8, 8, 5). No cluster contains both labels, so clonality is
label-consistent, matches Jones et al. (single-mating-type clonal groups), and does not inflate the 1:1 (it holds after clone-correction, 60:62). Duplicates: Sn79 is in the set twice
(GCA_002267025.1 NDSU, GCA_051013865.1 ANU, same BioSample SAMN07500996, mash distance about 0, both MAT1-2); SN15 sits in a cluster of 8. Strain names are otherwise unique (185 of 186).
The 1:1 is a population-level result; the *P. avenae* and non-WA subsets are too small to test.

## 4. Root cause, generalisation, safety of dropping accessory genes

Exact cause (all 87 wrong loci are identical): the vote did not use a MAT1-1-3 bitscore win. Sequence of events in v0.6.0 (`GCA_002267005.1`, report lines ~96-190):
1. Fused references: v0.6.0 had no Pleosporales record (the Pleosporales/Mycosphaerellales records came in `7c829bf`/`d9a2fac` after the tag). The HMG gene hit `5518_3639_MAT_combined`
   (*Fusarium graminearum* 3639) gene **MAT1-1-3**, which is an HMG-box gene in the Sordariomycete MAT1-1 idiomorph, and `28548_06cwm-g22_MAT_combined`
   (*Sclerotinia trifoliorum*) gene **MAT1-2-1**.
2. `resolve_idiomorph_overlaps` (src/MATPredict/detect/idiomorph.py) treats MAT1-1-3 and MAT1-2-1 as mutually exclusive and picks the winner on `_rank` = (proteome flag, first-pass identity):
   MAT1-1-3 41.071% against MAT1-2-1 40.0%, a 1.07-point margin, identical in all 87 loci (`basis: first_pass_identity`).
3. The loser MAT1-2-1 (55.1 bits) was marked superseded, so `idiomorph_candidates` held only MAT1-1: 45.4 (MAT1-1-3), confidence `medium`; the report's genes_found lists no MAT1-2-1.
The gene named MAT1-1-3 in the record is the real HMG-box MAT1-2-1 of *Parastagonospora*; the record itself documents (order.yml lines 83-93) that MAT1-1-3 and MAT1-2-1 are HMG paralogs that cross-match.
Records run: `13684_sn436ga98` (MAT1-2-1, 100% identity, 345 aa) beats MAT1-1-3 by identity (90 loci), `13684_sn435pl98` (MAT1-1-1) beats MAT1-2-1 weak hits (88 loci), and confidence is `high` in all 182.

Important correction to the replay in `analysis/2026-10-16_ascomycota-core-profiles.md`: `replay_vote.py` reads `e.get("lost_to")`, a key that does not exist in the reports, so superseded hits still vote
in the replay. The replay therefore credits "drop MAT1-1-2/-3 from the vote" with fixing the 87 loci. In the real code the MAT1-2-1 hit has already been superseded, so dropping MAT1-1-3 from the vote leaves
MAT1-1-1 (31.2 bits) alone and the call stays MAT1-1 (wrong). The same correction applies to its "0 lost" and "186 gained" figures. I checked this on the 87: with superseded hits removed and accessory genes dropped,
the call stays MAT1-1; with superseded hits kept (as in the replay) it becomes MAT1-2 (derived from the v0.6.0 report's `idiomorph_candidates` and `genes_found`, not re-run through the pipeline).
A fix must change the resolution step too, not only the vote.

Proposed change (not applied):
```
# idiomorph.py
-    for gene_a, gene_b in _mutually_exclusive_pairs(family, hits):
+    core = {g["name"] for g in family.genes if not g.get("optional")}        # MAT1-1-1, MAT1-2-1, APN2, SLA2 ...
+    for gene_a, gene_b in _mutually_exclusive_pairs(family, hits):
+        if not (gene_a in core and gene_b in core):    # an optional gene may not supersede a core gene's hit
+            continue
...
# idiomorph_candidates: skip optional genes (MAT1-1-2, -3, -4, -5, MAT1-2-10) unless no core gene has a hit.
```
Then re-replay with hit-level superseded handling on the full campaign before any release.

Generalisation (`campaign_pattern.py`, v0.6.0 vs records-run reports; pattern = MAT1-1-3 wins a cross-match over a higher-scoring MAT1-2-1, 3,118 v0.6.0 loci overall):
- Dothideomycetes with a higher-scoring MAT1-2-1 beaten by MAT1-1-3 and no strong alpha gene, v0.6.0: *Parastagonospora* 89 (87 wrong), *Pyrenophora* 59 (label MAT1-1, unchanged by records; with records the alpha gene scores at least 100 bits, so the call is now supported),
  *Corynespora* 10 (MAT1-2, unchanged), singles in *Plenodomus*, *Sporormia*, *Coniothyrium*, *Lentithecium*, *Ramulariopsis*. In the records run: *Pyrenophora* 59 (alpha supported), *Corynespora* 10, *Plenodomus* 1, *Ramulariopsis* 1, *Lentithecium* 1 remain; *Parastagonospora* 0.
  So the records fix works where the order has its own records, and *Corynespora* (10), *Plenodomus* and *Ramulariopsis* have not been checked against any truth.
- Outside Dothideomycetes the pattern is common (e.g. *Pyricularia* 197, *Fusarium* 433, *Calonectria* 78, *Verticillium* 29; 489 loci with no alpha gene of 100 bits or more outside the records-run orders in the whole campaign), where MAT1-1-3 is a genuine MAT1-1 gene. I did not examine these; the label is right where a core MAT1-1-1 is also found.
  Orders without records (Botryosphaeriales, Dothideales and others in the full Dothideomycetes run) show no such loci in the records run.
- Is a core-only vote safe? On the 182 Parastagonospora loci, voting on MAT1-1-1 and MAT1-2-1 (unsuperseded and superseded alike) reproduces the records-run label in 182 of 182 in both runs. Beyond that I measured nothing new: the earlier
  replay's 974 scored loci are affected by the superseded-hit bug above, and a core-only vote cannot call a locus whose only evidence is MAT1-1-2/-3 (the replay leaves those unchanged). Do not enable it without the resolution change and a full-campaign replay with correct supersession.

## Files
`gene_audit.py`/`.tsv`, `summarize.py`, `analyze.py`, `locus_summary.py`/`.tsv`, `campaign_pattern.py`, `campaign_accessory_wins*.tsv`, `genomes_table.tsv`, `para_asm.tsv`, `para_bs.tsv`, `para_mash.tsv.gz`,
`tblastn_none_genomes.tsv`, `assembly_stats_none_and_controls.txt`. HPCC jobs `para-audit`, `para-mash`, `para-tblastn` (outputs in `/bigdata/stajichlab/jstajich/projects/core-hmms`).
