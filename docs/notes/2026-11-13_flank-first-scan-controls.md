# Flank-first MAT scan: Pfam models and the Rhizophagus control set
Status: in progress (controls only; no genome scan has run)

## Question
Can adjacent-gene (flank) conservation find MAT loci in Mortierellomycotina and Glomeromycotina, which have no curated MAT locus?

## Inputs
- Pfam 38.2 (`/srv/projects/db/pfam/current/Pfam-A.hmm`): HMG_box (PF00505.26), Homeodomain (PF00046.36). I also fetched HMG_box_2 (PF09011.17); it was not requested.
- Control set: 93 proteins linked to PMID 27572831 (Ropars et al. 2016, Nat Microbiol, doi 10.1038/nmicrobiol.2016.33) in NCBI Protein. Fetched on 2026-11-13 with the elink and efetch calls.
  - HD2: 41 (7 titled "hd2", 34 "HD2 AMF mating type protein").
  - HD1-like: 39.
  - Flank genes: choline transporter-like 6, phosphoglycerate mutase family 7.
- The full paper text was read on 2026-11-13 after the control set was built (see 2026-11-13_amf-hd-locus.md).

## Method
`hmmsearch --cut_ga` of the 93 proteins against the three models. Second run with E <= 1.

## Results
| Protein group | n | Pfam hit (GA cutoff) |
|---|---|---|
| HD2 | 41 | Homeodomain, 41/41 |
| HD1-like | 39 | none, 0/39 (also none at E <= 1) |
| choline transporter-like | 6 | none |
| phosphoglycerate mutase family | 7 | none |

## Limits
- "HD1-like" proteins have no Pfam homeodomain or HMG hit. A homeodomain-pair search keyed on Pfam would find HD2 only. The HD1-like protein needs a sequence search (diamond) instead.
- The 93 records are alleles from several isolates. They are not 93 independent loci.
- Whether this locus controls mating is not established by these data. The authors call it putative.

## Files
`results/2026-11-13_flank_first_scan/pfam/hmg_homeo.hmm`, `.../ropars/` (ropars_proteins.faa, domains.tbl, domains_E1.tbl, prot_summary.json).
