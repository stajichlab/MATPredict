# GenBank input check (2026-11-04)

Question: does `matpredict detect` work on a GenBank file, both **annotated** (features and CDS translations) and **unannotated** (sequence only)?
Test genome: *S. cerevisiae* S288C RefSeq (GCF_000146045.2, 17 records, 12.3 Mb; MAT, HML and HMR on NC_001135.5). Annotated = the NCBI GBFF as is (33 MB);
unannotated = the same file with every feature except `source` removed (15 MB, 0 CDS). MATPredict `main` (afcd7ea), GmG v0.3.1, HPCC job 31661658
(`GmG_runs/matpredict_phylum/gb_out/`). Each input was run with `--phylum Ascomycota` and with `--phylum auto`, next to the same genome as plain FASTA.

| input | route | CDS proteins used | loci found | report |
|---|---|---|---|---|
| annotated GBFF | `--phylum Ascomycota` | 6,021 (proteome route) | 3 MATsc | report, GFF3, `tracks.yaml` (gene-model track) |
| annotated GBFF | `--phylum auto` | 6,021 | 3 MATsc | same |
| unannotated GenBank | `--phylum Ascomycota` | 0: "searching the genome sequence (no gene-model track)" | 3 MATsc | report, GFF3 |
| unannotated GenBank | `--phylum auto` | 0 | 3 MATsc | same |
| FASTA `.fna` / `.seq` | both routes | n/a | 3 MATsc | same |

- **Every input finds the same three loci**: HML (MATalpha) and MAT (MATalpha, medium confidence) and HMR (MATa, high), with the same idiomorph and confidence. Spans differ by a few
  bp between annotated (gene models from the file) and sequence-only runs (for example 199,546-200,969 against 199,549-200,966); the 8-family `--phylum` run reports a wider HML span
  (12,389-36,302) than the 1-family `auto` run, for FASTA and GenBank alike, so that is routing, not GenBank.
- `--phylum auto` works on GenBank: the file is converted to FASTA before routing, so GmG sees the sequence (`routed by family taxid 4893`, 1 family searched, 12-30 s against 42-61 s).
- A FASTA named `.seq` is accepted by `detect` (it has no extension filter) and by `auto` (GmG gets it through an `.fna` link, #248).
- The annotated file reports 6,021 CDS proteins; the page's converter counts 6,027 CDS features (the difference is features without a `/translation`).
- **Browser page converter** (`genbank_input.mjs`, Node, same files): annotated gives 17 contigs and 6,027 CDS, unannotated 17 contigs and 0; all 17 sequences equal NCBI's FASTA.
  The full browser run (WebAssembly tools) was not run here; it needs the CI bundle.
