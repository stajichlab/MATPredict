# 020. Annotation errors at MAT loci that mislead homology-based detection

- **Category:** assembly-or-annotation artefact
- **Status:** verified (each entry carries its evidence code)
- **Lineage:** various

## Summary
Curation found missing, truncated, mislabelled and misnamed MAT genes in public
records. Several would silently mislead a search by gene name or annotation.
The full list is `ANNOTATION_ERRORS_FIXED_REPORT.md`; the most citeable are below.

## Evidence (IDs from `ANNOTATION_ERRORS_FIXED_REPORT.md`)
- A1 / B3.1: D. hansenii CBS767 MTLalpha1 XP_460134.1 named only
  "DEHA2E19096p"; the record wrongly said no alpha1 exists (entry 004).
- B1.6: U. maydis mfa1 (41-aa pheromone precursor) absent from the modern RefSeq
  annotation (NC_026482.1), present in the 1995 deposit U37795.
- B2.4: Morchella importuna KY782629.1 / KY782630.1 APN2 split into two CDS.
- B4.1, B4.2: C. albicans assemblies carrying one idiomorph while reads show
  both (entry 011).
- B5.2 (curation-puccinio): Melampsora larici-populina `GL883124.1`: RefSeq
  names EGG03504.1 "MlpbE1", but it is the ortholog of P. graminis bW1
  PGTG_05144 (30% over 248 aa, E=2e-31); EGG03439.1 is the bE1 ortholog (36%
  over 218 aa, E=9e-26). The bE label sits on the bW ortholog.
- B5.4 (curation-puccinio): Wallemia SXI1 not annotated (entry 001).
- C1: Yarrowia AJ617307.1 declares strain W29, but W29 lacks the matb sequence.
- C5, C8, C9 (curation-mucor-dothideo): Zymoseptoria, Leptosphaeria and
  Pseudocercospora MAT deposits carry no strain qualifier.
- D3, D4, D5, D6, D7: biological absences or non-determinant genes that look
  like detection failures (Colletotrichum MAT1-2 in both partners; CTG-clade
  MTLalpha2 absence; C. parapsilosis MTLa1 pseudogene; L. elongisporus and
  C. sojae without MAT genes; Hanseniaspora MAT genes present, 15/18 no-call
  genomes carry a core hit next to SLA2).
- D10 (curation-mucor-dothideo): Lodderomyces beijingensis GCF_963989305.1,
  MTLa1+MTLa2-like blocks at 51-53% on 7 chromosomes (entry 021).

## Method that found it
Curation against deposits and papers; reciprocal best hits; read mapping.

## Verification done / still open
Upstream corrections not filed.

## Limits
Collected opportunistically during curation, not by a systematic audit.

## Related
Entries 001, 004, 011, 021.
