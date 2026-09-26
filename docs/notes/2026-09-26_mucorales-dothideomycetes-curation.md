# Curation round 2026-09-26: Lichtheimiaceae, Syncephalastrum, Dothideomycetes

Branch `curation-mucor-dothideo` (from `polish-scope-cuts` at `f81dad1`).
Curator request: curate these three targets "if realistic".

## Dothideomycetes: four tier-1 records added (commit `dca6ccc`)

| record | deposit | genes | source |
|---|---|---|---|
| `5016_c5_MAT_MAT1-1` | AF029913.1 | MAT1-1-1 | Turgeon et al. 1993, PMID 8479433 |
| `5016_c4_MAT_MAT1-2` | AF027687.1 | MAT1-2-1 | Turgeon et al. 1993, PMID 8479433 |
| `1047171_ipo323_MAT_MAT1-1` | AF440399.1 | MAT1-1-1 | Waalwijk et al. 2002, PMID 11929216 |
| `1047171_ipo94269_MAT_MAT1-2` | AF440398.1 | APN2, MAT1-2-1 | Waalwijk et al. 2002, PMID 11929216 |

- Evidence: *Cochliobolus* MAT identity by complementation, RFLP and genetic
  mapping. *Zymoseptoria* by absolute F1 cosegregation with mating type.
- Validation passes: 5/5 proteins at 100% identity and coverage. 649 tests pass.
- The records are `accepted` on this branch only, so detect could use them.
  Each record says curator sign-off is pending.
- Pleosporales flanks (GAP1, BGL1, ORF1) are not roster genes and are not curated.

### Effect on the 100-genome Dothideomycetes pilot list

Same code, database differs. `results/2026-09-26_dothideo_curation/`
(`compare.py`, `compare_output.txt`).

| | before (`f81dad1`) | after (`dca6ccc`) |
|---|---|---|
| genomes called | 90 | 91 |
| high / medium loci | 46 / 45 | 47 / 45 |
| idiomorph flips | | 0 |
| median idiomorph margin, Pleosporales | 9.0 (n=19) | 18.7 (n=22) |
| median idiomorph margin, Mycosphaerellales | 20.1 (n=8) | 26.5 (n=8) |
| median wall per genome | 480 s | 468 s |

Two genomes changed: *Polyplosphaeria fusca* (Polfu1) MAT1-2 medium -> high,
and GCA_060176455.1 (Pleosporales) none -> MAT1-2 medium. The 9 uncalled
genomes did not change. The records mainly widen the margin between the two
idiomorphs. They add little recall, because recall was already 90%.

### Not added, curator decision

- *Diplodia sapinea* `KF551229.1` / `KF551228.1` (Bihon et al. 2014, PMID
  24220137). Heterothallic by genome content and a 1:1 idiomorph ratio in
  populations. No sexual state has been observed, so there is no cross. It
  is weaker than the Morchella precedent, which had single-ascospore
  segregation. Botryosphaeriales is already called 8/8.
- Other tier-1 deposits in `docs/notes/2026-09-24_mat-reference-gap-literature.md`
  (*Leptosphaeria maculans*, *Parastagonospora nodorum*, *Pseudocercospora
  fijiensis*, *Fulvia fulva*) were not added. They are the next candidates if
  more Dothideomycetes coverage is wanted.

## Lichtheimiaceae and Syncephalastrum: nothing meets the bar

- NCBI nuccore has no MAT or sex-locus deposit for *Lichtheimia*,
  *Rhizomucor*, *Thermomucor* or *Syncephalastrum* (searches for sexP, sexM,
  "mating type", "sex locus", HMG). Only WGS contigs.
- PubMed has no paper that describes a sex locus in these genera. The
  *L. corymbifera* genome paper (Schwartze et al. 2014, PMID 25121733) mentions
  sexP/sexM only in general terms.
- So no tier-1 record is possible. A tier-2 record needs a gene tree that puts
  the candidate HMG gene in the sexM/sexP clade (the 2026-09-24 ruling). The
  sexM/sexP phylogeny spec
  (`docs/superpowers/specs/2026-09-26-sexM-sexP-candidate-phylogeny-design.md`)
  would provide it.

### Flank order in the annotated genomes checked

The CDS in the scan's withheld loci were blasted against the curated
Mucorales MAT proteins (blastp, E <= 1e-5):

- *S. racemosum* NRRL 2496, MCGN01000001.1:1,493,958-1,545,990 (52 kb):
  tptA ortholog ORZ02716.1 (71% identity to curated tptA) at the left end.
  An HMG-box gene ORZ02736.1 at the right end, about 50 kb away, with 19
  genes between them. Its best MAT hit is sexP at 26.5% over 83 aa (the HMG
  box only). No rnhA in the window.
- *S. racemosum*, MCGN01000004.1:1,379,847-1,471,236: rnhA-like ORY97662.1
  (29% over 463 aa) and an HMG gene ORY97701.1 (sexP 30% over 77 aa) about
  86 kb apart. No tptA.
- *L. ornata* CBS 291.66, NW_026695375.1:118,280-206,612: HMG gene
  XP_058339816.1 (sexM 35% over 105 aa). The only tptA-like CDS nearby,
  XP_058339812.1, matches tptA at 25%, which is too low to call it the
  ortholog.
- *L. ornata*, NW_026695325.1:147,273-186,138: rnhA-like XP_058347066.1 (31%)
  and an HMG gene XP_058347053.1 (sexP 28% over 106 aa), about 37 kb apart.

This agrees with the synteny check (`results/2026-09-26_flank_synteny_ED/`):
tptA and rnhA are not within 100 kb of each other in these families. In the
windows checked, tptA and rnhA each sit next to a different HMG-box gene. None
of the HMG-box genes matches sexM/sexP outside the HMG box. The standard
tptA-sex-rnhA roster does not fit these genomes. Do not choose flank genes for
a Lichtheimiaceae or Syncephalastraceae record until the gene tree places one
of these HMG genes in the sexM/sexP clade.
