# Evidence for docs/notes/2026-10-21_gene-model-extension-design.md
Prototype and measurements from a read-only review (Fable agent, 2026-10-21). All run against the patched run `runs190r` of PR #117 (3,181 reports).
- `survey.tsv`: per core-gene model (all reports): exonerate/miniprot lengths, internal stops, in-frame continuation past the model end (`*_ext_aa`, `*_ext_state`), upstream-to-stop and first-M offsets.
- `modelkit.py` (survey), `modelkit2.py` (chooser prototype; MAFFT identical-columns surrogate score, NOT the pipeline's BLOSUM62 score), `nulltest.py` (shuffle test of the 3' extensions), `peek.py` (Botrytis).
- `test2_out.tsv`: chooser output on 35 genomes (the cases the thresholds were fitted on, so a demonstration, not a validation).
- `refs_MAT1-1-1.faa`: the 26 curated MAT1-1-1 proteins used as references.
- `fable_review_report.md`: the full review report as delivered. Its claims are labelled measured/read/inferred; see the design note for what is carried forward.
Scripts expect the HPCC layout (genomes in BFD `input_clean_genomes`, the pixi `test` env); paths are not portable.
