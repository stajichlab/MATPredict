# CTG-clade yeasts without a polished MTLA1 model: what v0.9.0 shows (#260)

Date 2026-11-21. Tables in `results/2026-11-21_ctg_mtla1_260/`. The 2026-10-23 diagnostic folder behind the issue's "55 yeasts" has been purged, and the v0.6.0 and v0.7.0 reports were not kept, so this measures the v0.9.0 campaign (commit 7d73be5, v0.9.1 preset) and does not reproduce the 55.

## Where MTLA1 stands in v0.9.0
1,661 MTLA1 gene entries in 1,533 genomes with an MTL-family locus.
- **130 are `lost_to` MTLalpha2** (127 of them `not_polish_candidate`, 3 `polished_single`): MTLA1's hit was taken by the alpha2 gene, mostly in *Candida* (108). That is a cross-gene assignment, not a lost model, and it is left out below.
- **1,531 are live:** 932 `polished_agree`, 236 `polished_single`, 143 `polished_disagree`, **220 `unpolished`** (14%). Unpolished means the gene was a polish candidate but neither exonerate nor miniprot made a model from the raw tblastn hit.

## The 220 unpolished entries
- 85 species. *Lodderomyces* 104 (96 *L. parapsilosis*, 8 *L. metapsilosis*, so many genomes of one species), *Metschnikowia* 29, *Sungouiella* 15, *Yamadazyma* 13, *Scheffersomyces* 12, *Kurtzmaniella* 7.
- Their hits are weak: median identity 31.1% and median coverage 31.0%, against 100% and 100% for `polished_agree`. Tiers: 177 medium, 35 high, 8 low. None is `defective`.
- These genera have almost no `polished_agree` at all, even where a model exists: *Lodderomyces* 0 agree / 27 single / 104 unpolished; *Metschnikowia* 0 / 84 / 29; *Kodamaea* 0 / 18 / 1 (`mtla1_status_by_genus.tsv`). Against *Candidozyma* 680 of 690 agree. So the 220 are the weak end of a lineage where the references fit badly.

## Probe: is there a model the pipeline missed?
One genome per species for 24 of those species (five genera, `probe.tsv`), plus four controls where the pipeline polished MTLA1.
- **Controls:** miniprot with all MTLA1 references of `db/` finds MTLA1 at the same locus in 4 of 4 (*C. albicans* ×2 and *C. dubliniensis* at 87-99%, *D. hansenii* at 71%). The Pfam homeodomain model (gathering cutoff) fires in 3 of 4; it misses *D. hansenii*. So a missing Pfam hit is not evidence of absence.
- **Unpolished genomes:** miniprot finds **no alignment at all** to any of the references in 24 of 24, at a 5 kb intron limit.
- **Homeodomain at the pipeline's own hit:** 15 of 24 (all the *Scheffersomyces*, most *Metschnikowia*, *Sungouiella blattae*, *L. metapsilosis*; identity 40-50% over about 20% of the protein, which is the homeodomain alone). 9 of 24 have none (*L. parapsilosis*, *S. intermedia*, 5 *Metschnikowia*, 2 *Yamadazyma*).

## What this says about the issue's causes
- **Missing or distant reference record: consistent.** The database has no a1 record from *Lodderomyces*, *Metschnikowia*, *Scheffersomyces* or *Sungouiella*; the nearest reference is a *C. albicans*, *D. hansenii* or *C. auris* record, and an independent miniprot search with all of them finds nothing. In 15 genomes the only conserved part, the homeodomain, is what the weak hit carries.
- **Subphylum guard:** not the route for these entries; they are live in the report (no `lost_to`).
- **Code change since v0.6.0:** not tested (no old reports).
- **CTG genetic code and the 5 kb intron cap:** not tested separately. The miniprot probe uses the standard code and a 5 kb limit and also finds nothing; the unpolished hits are single-exon fragments. I did not rerun the pipeline with these toggled.
- **Not shown:** that these genomes have an MTLA1 at all. For the 9 genomes with no homeodomain (including *L. parapsilosis*, 96 genomes) the hit may be noise.

## For the curator
Leaving these as they are is consistent with the evidence: they are called at medium tier on the weak hit and the pipeline does not claim a model. Making them polished would need a real MTLA1 record from one of these lineages (as for C. auris MTLA2), taken from a verified annotation or deposit; it cannot come from the raw hits. Say if you want me to look for sources in *Metschnikowia*, *Scheffersomyces* or *Lodderomyces*, or to close #260 as measured.
