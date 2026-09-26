# Basidiomycota MAT anchors, and the first Russulales, Boletales, Polyporales references

2026-09-26, branch `basidio-anchors` (from `polish-scope-cuts` f81dad1), commit
`1403399`. Results: `results/2026-09-26_basidio_anchors/` (main checkout).

Curator request: define Basidiomycota anchors before the 3,275-genome launch,
test them on diverse genomes, and curate Boletales, Polyporales and Russulales
references if realistic.

## 1. Literature: which anchors exist

* **Agaricomycetes, HD locus: MIP 5' and beta-fg 3'.** "A gene organization that
  is consistent for nearly all characterized Agaricomycetes with limited
  exceptions" (Mujic et al. 2017, G3 7:1775, PMID 28450370). MIP is "consistently
  syntenic" with MAT-HD in Polyporales and other Agaricomycetes, with no linkage
  to the P/R genes (James et al. 2013, Mycologia 105:1374, PMID 23928418).
  Exceptions: *Lentinula edodes*, MIP not closely linked (PMID 24295887);
  *Sparassis*, beta-fg absent (PMID 33213780).
* **PR locus: no conserved anchor found in the literature.** In Boletales,
  isoprenyl cysteine methyltransferase (ICMT) genes flank the B locus, but their
  copy number varies with breeding system and their orthology is not
  established (Mujic et al. 2017).
* **Ustilaginomycotina and Pucciniomycotina: no anchor found in the literature
  searched.** The *Leucosporidium scottii* paper (Maia et al. 2015, PMID 26178967)
  was not available as full text through PubMed.

## 2. Positional test on 30 genomes

`find_anchors.sh`, `analyze_anchors.py`. For each genome, tblastn (e<=1e-5) of
MIP1 and beta-fg (both from *A. bisporus* H97), *S. cerevisiae* STE14 (ICMT),
54 curated HD-class proteins and 19 curated receptor proteins. Anchor = the
top-bitscore HSP of that query. "Adjacent" = an HD HSP within 20 kb on the same
contig. Background = the fraction of the genome within 20 kb of any HD cluster
(the chance that a random point is "adjacent").

| subphylum | genomes | HD next to MIP1 | HD next to beta-fg | receptor next to ICMT | background |
|---|---:|---:|---:|---:|---:|
| Agaricomycotina | 18 | **12** | 9 | 0 | 0.1% |
| Ustilaginomycotina | 4 | 0 | 0 | 0 | 0.4% |
| Pucciniomycotina | 7 | 0 | 0 | 0 | 0.2% |
| Wallemiomycotina | 1 | 0 | 0 | 1 | 0.4% |

Stable at 10, 50 and 100 kb (Agaricomycotina 10, 12, 13 of 18; all other
subphyla 0). Within Agaricomycotina, the 6 non-adjacent genomes are
*Cryptococcus* (Tremellomycetes, different MAT architecture), *Dacryopinax*
(Dacrymycetes), *Rhizoctonia* (Cantharellales), *Auricularia*, and *Ganoderma*
and *Russula*, which had no HD HSP at e<=1e-5. So in the orders now in the HD
family's scope, 11 of 13 genomes put an HD hit next to MIP1.

**Search for other anchors outside Agaricomycetes** (`neighbors.py`,
`analyze_anchors2.py`). Neighbours within 30 kb of the best HD hit in 10
annotated genomes were compared across genomes. Candidates were then tested by
tblastn on all 30 genomes:

* A DEAD/H helicase, an IES1 (INO80) subunit and two uncharacterised proteins sit
  next to bE/bW in *M. maydis*, *S. reilianum* and *Malassezia*. They are
  adjacent in 2/2 Ustilaginales/Malasseziales genomes (one is the query's own
  source), and in 0/2 Exobasidiomycetes.
* Splicing factor SF3B5 sits next to the best HD hit in 3 Microbotryomycetes
  proteomes, but by tblastn it is adjacent only in its own source genome.
* The protein family shared next to HD in *Mixia*, *Sporobolomyces*,
  *Leucosporidium* and *Wallemia* is long (900-1,285 aa) and weakly similar to
  bW. It is most likely the divergent HD1 partner, not a flank.

**Limit:** outside Agaricomycotina the HD hits themselves are uncertain (31-42%
identity to bE). The absence of adjacency there can mean "no anchor" or "the HD
hit is not the MAT gene".

## 3. Curated records (all genes 100% identity and coverage)

| record | order | deposit | genes | citation |
|---|---|---|---|---|
| `984962_um274_HD_A1` | Russulales | KF280353.1 | MIP1 (5' partial), HD1, HD2, HD1, HD2 | PMID 23864721 |
| `90004_unknown_HD_A1` | Boletales | AB646132.2 | MIP1 (5' partial, codon_start 3), HD1, HD2 | doi:10.1007/s11557-012-0840-z (no PMID) |
| `2822231_sb25_HD_A1` | Polyporales | HQ188438.1 | MIP1 (5' partial, codon_start 3), HD2, HD1 | PMID 21131435 |
| `192523_h97_HD_A1` v2 | Agaricales | NW_006267344.1 | + beta_fg XP_006454074.1, + MIP1 XP_006454076.1 | PMID 34356095 |

The *A. bisporus* MIP1 is annotated "hypothetical protein". Its identity rests
on homology: 38.9% over 717 of 757 aa to yeast OCT1 (P35999, E=1e-160).

Not curated:
* *Suillus luteus* ON315855.1 (HD1+HD2, PMID 37070772). The HD alleles were
  reconstructed from genome reads by local assembly, with no crosses or targeted
  sequencing. This is closer to tier 2.
* *Ganoderma boninense* HD1/HD2 and STE3 (ON855036-7, ON212092...). These are mRNA
  deposits with no genomic locus.
* *Coriolopsis trogii* receptors (MF990238-41). These are single-gene records, and
  there is no PR family outside Agaricales to attach them to.

Roster: `MIP1` and `beta_fg` are optional `flanking_conserved` genes of `HD`.
The HD scope is widened from Agaricales to Agaricales, Russulales, Boletales
and Polyporales (order taxids verified against NCBI).

## 4. Pilot detection, three arms, same 30 genomes

`compare_arms.py`, `compare_arms.txt`. Arms: **base** = f81dad1 (before this
work); **noanchor** = 1403399 with MIP1/beta_fg `exclude_from_search`;
**anchor** = 1403399.

| arm | Agaricomycotina called | with an anchor gene in the call | median s | Ustilago called | Puccinio called |
|---|---:|---:|---:|---:|---:|
| base | 8/18 | 0 | 204 | 4/4 | 1/7 |
| noanchor | 13/18 | 0 | 54.5 | 4/4 | 1/7 |
| anchor | 13/18 | **13/13** | 74 | 4/4 | 2/7 |

* **The new records and the wider scope give the call gain (8 to 13) and the
  speed-up (median 204 to 54.5 s).** The speed-up is from routing: 9 genomes in
  Russulales/Boletales/Polyporales moved from `phylum_fallback` (all 9 families)
  to `lineage` (HD only). Newly called: *Serpula*, *Trametes*, *Gelatoporia*,
  *Ganoderma*, *Stereum*.
* **The anchors change no call in scope. They add positional corroboration.**
  Every one of the 13 Agaricomycotina calls carries MIP1 and/or beta_fg. HD is
  already `high` on core genes alone, so an optional flank cannot promote it.
  Cost: median 54.5 to 74 s (+36%).
* **Outside scope, the anchors carried one unverified call.** *Rhodotorula
  toruloides* (phylum_fallback) gains an HD `mat_locus`, medium, HD1+HD2+MIP1,
  over 40 kb. The positional test found no HD within 20 kb of its top MIP1 hit.
  This call has not been checked, and it should not be trusted. Under the
  2026-09-26 `not_searched` ruling this genome would not be searched by default.
* *Ganoderma boninense* G3 has two complete MIP1-HD-beta_fg blocks 276 kb apart
  on one chromosome (CM035305.1). This could be two haplotypes in one
  assembly, or a duplication. Not checked.

## 5. Defects seen, not fixed here

* **`homothallic_candidate` fires on normal heterothallic Basidiomycota loci.**
  8 of 30 genomes (anchor arm; 5 in base) get the class. The cause is that
  HD1+HD2 (and bE+bW) are two different `gene_class` values, which is what
  338ee18's rule takes as two unrelated MAT genes. In Basidiomycota the HD1/HD2
  pair *is* one heterothallic locus. The same applies to receptor+pheromone.
  This must be fixed before the phylum launch, or the homothallism screen will be
  mostly false in this phylum.
* *Cryptococcus neoformans* (GCA_002221985.1) gets 0 calls in all arms under
  lineage routing to Tremellales MAT. Not investigated here.

## Curator decisions needed before the launch

1. **Idiomorph labels for unnamed multiallelic HD alleles.** Heterobasidion and
   Phanerochaete carry a placeholder "A1". Choose a convention (for example,
   widen the pattern to allow a strain label).
2. **Should anchors be searched only in scope?** Under `phylum_fallback` they
   reached *Rhodotorula*, where they are not supported. The `not_searched` ruling
   may make this moot.
3. **Keep the anchors at +36% runtime for corroboration only?** They do not change
   any call in scope. Their value is to back each call with position, and in
   future to rescue loci the bar withholds in uncurated Agaricomycete orders
   (Cantharellales, Auriculariales, Hymenochaetales). This pilot did not measure
   that.
4. **Fix `homothallic_candidate` for Basidiomycota** (section 5) before launch.
5. **No anchors exist yet for Ustilaginomycotina, Pucciniomycotina or PR loci.**
   The Ustilaginales/Malasseziales helicase/IES1 neighbourhood is a candidate
   for the existing `bLocus` family only.
