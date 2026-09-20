# MAT Locus Visualization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give MATPredict two real, scriptable diagram outputs: a single-locus gene-structure figure (via pyGenomeViz) and a multi-locus synteny/comparison figure (via clinker), for both curated records and real detection results — closing the real "all-N placeholder sequence, no CDS features" gap in the existing GenBank/GFF3 writers that currently blocks either tool from working.

**Architecture:** Fix the two existing writers to emit real sequence + CDS features first (Tasks 1-2), since both downstream tools depend on that data existing at all. Then a small investigation task installs and hand-trials both tools against the newly-fixed real output to resolve the real open questions the prior research left unconfirmed (Task 3), before committing to the exact CLI wiring the trial's findings inform (Tasks 4-5). A final task produces and commits one real example of each diagram type from real curated data, as both a smoke test and a durable example.

**Tech Stack:** Python, Biopython (`Bio.SeqRecord`/`Bio.SeqFeature`, already a dependency), `clinker` (new, `bioconda::clinker-py`), `pyGenomeViz` (new, `conda-forge::pygenomeviz`), this project's existing `NcbiClient`/`_independent_translation`/`gff_export`/`report.py` modules.

**Spec:** `docs/notes/2026-09-18_mat-locus-visualization-research.md` — read this in full before Task 1. It already resolves the tool choice (clinker for group/synteny, pyGenomeViz for single-locus) and identifies the exact writer gap this plan starts by fixing. It also leaves 2 real open questions for Task 3 to resolve empirically rather than by further reading: clinker's exact CDS/translation feature requirement, and whether pyGenomeViz's GFF3 parser tolerates this project's fragmented-multi-contig-locus GFF3 shape.

## Global Constraints

- Never build a second NCBI client — reuse `NcbiClient.fetch_nucleotide_sequence` (already used by `db/validate.py`'s `_independent_translation`) for any real segment sequence fetch this plan needs.
- Never fabricate a sequence — if a segment's real nucleotide sequence can't be fetched (e.g. an assembly-level `GCA_`/`GCF_` accession `NcbiClient` can't resolve yet, a known, documented limitation elsewhere in this codebase), fall back to today's `N`-placeholder behavior for that segment rather than inventing data, and say so in the output (e.g. a feature qualifier or log warning), matching this project's "never fabricate, disclose the gap" convention used throughout its curation work.
- Add `clinker`/`pygenomeviz` to `pixi.toml` under the existing `[dependencies]` table (both channels — `conda-forge`, `bioconda` — are already configured; no new channel needed).
- No change to `run_pipeline`'s detection behavior — this plan only adds NEW output alongside existing GFF3/YAML reports, never modifies detection logic itself.
- `db/candidates/` records are out of scope for any new CLI command this plan adds (matching every other DB-wide tool in this project).

---

### Task 1: Fix `gff_export.write_genbank` to emit real sequence and real CDS features

**Files:**
- Modify: `src/MATPredict/db/gff_export.py`
- Modify: `src/MATPredict/db/cli.py` (`build_gff_for_record`'s call site, if `write_genbank`'s signature changes)
- Test: `tests/db/test_gff_export.py`

**Interfaces:**
- Produces: `write_genbank(record: dict, sequences: dict[int, str], out_path: Path, ncbi: NcbiClient | None = None) -> None` — same signature as today plus one new optional `ncbi` parameter (defaulting to `None`, preserving today's placeholder-`N` behavior for any caller that doesn't pass one, e.g. existing tests). When `ncbi` is given, each segment's REAL nucleotide sequence is fetched via `NcbiClient.fetch_nucleotide_sequence` (falling back to the `N`-placeholder for that segment only if the fetch fails or the segment's `sequence_source.type` isn't `insdc_nucleotide` — never for the whole file), and each present gene with a real sequence available in the `sequences` dict gets a real `CDS` feature (not just `gene`) with a `translation` qualifier, alongside `role`/`gene_class`/`present_in_idiomorphs` qualifiers copied from the gene's own schema fields.

- [ ] **Step 1: Write the failing test for real-sequence fetching**

```python
# append to tests/db/test_gff_export.py
from unittest.mock import MagicMock

from MATPredict.db.gff_export import write_genbank


def test_write_genbank_uses_real_sequence_when_ncbi_client_given(tmp_path):
    record = {
        "record_id": "111_a_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 130,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC1.1", "seq_region": "ACC1.1"}},
        ]}},
        "genes": [
            {"gene_index": 0, "name": "G1", "role": "core_MAT", "present": True,
             "segment_index": 0, "start": 100, "end": 130, "strand": "+"},
        ],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.return_value = "ATG" * 10 + "TAA"

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={0: "M" * 10}, out_path=out_path, ncbi=fake_ncbi)

    fake_ncbi.fetch_nucleotide_sequence.assert_called_once_with("ACC1.1", 100, 130, None)
    text = out_path.read_text()
    assert "N" * 31 not in text  # the old placeholder is gone
    assert "ATGATGATG" in text.replace("\n", "").replace(" ", "")  # real sequence is present
    assert "CDS" in text
    assert "/translation=" in text.replace("\n", "").replace(" ", "")


def test_write_genbank_falls_back_to_placeholder_when_fetch_fails(tmp_path):
    record = {
        "record_id": "222_b_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 100, "end": 130,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC2.1", "seq_region": "ACC2.1"}},
        ]}},
        "genes": [],
    }
    fake_ncbi = MagicMock()
    fake_ncbi.fetch_nucleotide_sequence.side_effect = Exception("simulated NCBI outage")

    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={}, out_path=out_path, ncbi=fake_ncbi)  # must not raise

    text = out_path.read_text().replace("\n", "").replace(" ", "")
    assert "N" * 31 in text  # placeholder fallback, not a crash


def test_write_genbank_preserves_old_placeholder_behavior_when_no_ncbi_client_given(tmp_path):
    record = {
        "record_id": "333_c_MAT_combined",
        "locus": {"core": {"segments": [
            {"segment_index": 0, "start": 1, "end": 20,
             "sequence_source": {"type": "insdc_nucleotide", "accession": "ACC3.1", "seq_region": "ACC3.1"}},
        ]}},
        "genes": [],
    }
    out_path = tmp_path / "locus.gbk"
    write_genbank(record, sequences={}, out_path=out_path)  # no ncbi= at all -- default behavior

    text = out_path.read_text().replace("\n", "").replace(" ", "")
    assert "N" * 20 in text
```

- [ ] **Step 2: Run to verify it fails**

```bash
pixi run pytest tests/db/test_gff_export.py -v
```

Expected: FAIL (`write_genbank() got an unexpected keyword argument 'ncbi'`).

- [ ] **Step 3: Implement the fix**

In `src/MATPredict/db/gff_export.py`, add the `ncbi` parameter and rewrite the per-segment sequence-fetch and per-gene feature logic:

```python
def write_genbank(
    record: dict, sequences: dict[int, str], out_path: Path, ncbi: "NcbiClient | None" = None
) -> None:
    """Write a GenBank record for the core locus, from the same segments/genes data as write_gff3.

    Builds one Bio.SeqRecord per segment. When `ncbi` is given, each segment's REAL
    nucleotide sequence is fetched via NcbiClient.fetch_nucleotide_sequence -- the
    same mechanism db/validate.py's _independent_translation already uses -- and
    falls back to an all-"N" placeholder ONLY for that segment, on a fetch failure
    or an unfetchable sequence_source.type (e.g. an assembly-level GCA_/GCF_
    accession NcbiClient can't resolve yet), never fabricating a sequence. Passing
    no `ncbi` (the default) preserves the prior all-placeholder behavior exactly,
    for any caller/test that doesn't need real sequence.

    Each present gene with a sequence available in `sequences` gets a real `CDS`
    feature (not just `gene`) carrying a `translation` qualifier and `role`/
    `gene_class`/`present_in_idiomorphs` qualifiers copied from the gene's own
    schema fields -- the per-gene attributes clinker's `--colour_map`/
    `--gene_functions` (or pyGenomeViz's `--feature_type2color`) need to color/label
    by MAT-domain biology directly, without a separate manual mapping step.
    """
    segments = record["locus"]["core"]["segments"]
    genes_by_segment: dict[int, list[dict]] = {}
    for gene in record["genes"]:
        if not gene.get("present", True):
            continue
        genes_by_segment.setdefault(gene["segment_index"], []).append(gene)

    seq_records = []
    for segment in segments:
        segment_index = segment["segment_index"]
        segment_length = segment["end"] - segment["start"] + 1
        seq_region = segment["sequence_source"]["seq_region"]
        source = segment.get("sequence_source", {})

        nucleotide_sequence = None
        if ncbi is not None and source.get("type") == "insdc_nucleotide" and source.get("accession"):
            try:
                nucleotide_sequence = ncbi.fetch_nucleotide_sequence(
                    source["accession"], segment["start"], segment["end"], None
                )
            except Exception:
                nucleotide_sequence = None
        if not nucleotide_sequence:
            nucleotide_sequence = "N" * segment_length

        seq_record = SeqRecord(
            Seq(nucleotide_sequence),
            id=f"{record['record_id']}.segment{segment_index}",
            name=seq_region[:16] if seq_region else f"segment{segment_index}",
            description=f"{record['record_id']} core locus segment {segment_index} ({seq_region})",
        )
        seq_record.annotations["molecule_type"] = "DNA"

        for gene in genes_by_segment.get(segment_index, []):
            local_start = gene["start"] - segment["start"]
            local_end = gene["end"] - segment["start"] + 1
            strand = 1 if gene.get("strand") != "-" else -1
            location = FeatureLocation(local_start, local_end, strand=strand)

            gene_feature = SeqFeature(location, type="gene", qualifiers={
                "gene": [gene["name"]], "role": [gene["role"]],
            })
            seq_record.features.append(gene_feature)

            translation = sequences.get(gene["gene_index"])
            if translation:
                qualifiers = {
                    "gene": [gene["name"]], "role": [gene["role"]],
                    "translation": [translation],
                }
                if gene.get("gene_class"):
                    qualifiers["gene_class"] = [gene["gene_class"]]
                if gene.get("present_in_idiomorphs"):
                    qualifiers["present_in_idiomorphs"] = [",".join(gene["present_in_idiomorphs"])]
                cds_feature = SeqFeature(location, type="CDS", qualifiers=qualifiers)
                seq_record.features.append(cds_feature)

        seq_records.append(seq_record)

    SeqIO.write(seq_records, out_path, "genbank")
```

(`FeatureLocation` is already imported at the top of this file per its existing `write_gff3`/`write_genbank` code — confirm the exact existing import line and reuse it rather than adding a duplicate.)

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run pytest tests/db/test_gff_export.py -v
```

Expected: PASS, all tests including the 3 new ones.

- [ ] **Step 5: Update `build_gff_for_record`'s call site**

In `src/MATPredict/db/cli.py`, pass `ncbi=ncbi` into the existing `gff_export.write_genbank(record, sequences, out_path=...)` call (the `ncbi` client is already constructed and in scope there).

- [ ] **Step 6: Run the full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/db/gff_export.py src/MATPredict/db/cli.py tests/db/test_gff_export.py
git commit -m "feat: write real segment sequence and CDS/translation features in curated-record GenBank output"
```

---

### Task 2: Add a companion FASTA + CDS features to detection-result GFF3 output

**Files:**
- Modify: `src/MATPredict/detect/report.py`
- Test: `tests/detect/test_report.py`

**Interfaces:**
- Produces: `write_detection_gff3(outcome: DetectionOutcome, out_path: Path, genome_fasta: Path | None = None) -> None` — same signature plus one new optional `genome_fasta` parameter. When given, for every segment the function slices the REAL nucleotide sequence directly out of the genome FASTA (already on local disk at detection time — no network fetch needed, unlike Task 1's curated-record case) and writes it to a companion FASTA at `out_path.with_suffix('.fasta')`, using the SAME contig/seqid naming the GFF3 itself uses so clinker's GFF3+FASTA input convention (same base filename) is satisfied. Also adds a `CDS` feature (with a `translation` qualifier, from `GeneEvidence`'s data — check whether it already carries a translatable sequence after this session's earlier `GeneEvidence.exons` work, or whether a translation needs deriving the same way `benchmark.py`'s `_extract_translated_gene` already does for its own purposes; reuse that logic rather than writing new splicing/translation code if it's usable) alongside the existing `gene` feature.

- [ ] **Step 1: Read `report.py`'s real, current `write_detection_gff3` in full**, plus `benchmark.py`'s `_extract_translated_gene`/`_splice_transcript` (built in an earlier plan this session, already tested and reviewed for correctness including the minus-strand exon-order subtlety) — reuse that translation logic rather than re-deriving it, since `GeneEvidence.exons` (added in that same earlier plan) is exactly the data this needs.

- [ ] **Step 2: Write the failing test**

Construct a small, realistic `DetectionOutcome` fixture (a single `DetectionResult` with one `GeneEvidence` carrying real `exons`) plus a real small `tmp_path`-based genome FASTA (a short synthetic sequence with a real ORF at the recorded coordinates, matching this project's existing test-fixture conventions in `tests/detect/test_report.py` — read that file's existing fixtures first) and assert: (a) a companion FASTA is written at the expected path with the real sliced sequence, (b) the GFF3 gains a `CDS` feature line with a `translation=` attribute for the gene, (c) calling `write_detection_gff3` WITHOUT `genome_fasta` (the default) produces byte-identical output to before this task (no regression for existing callers).

- [ ] **Step 3: Run to verify it fails, then implement**, following the same "additive, default-preserving, real-fetch-with-graceful-fallback" discipline as Task 1 (a segment whose sequence can't be sliced — e.g. contig name mismatch — falls back to no CDS/no companion sequence for that segment, logged, never fabricated or crashing the whole write).

- [ ] **Step 4: Run tests to verify they pass; run the full suite; commit**

```bash
pixi run pytest -v
git add src/MATPredict/detect/report.py tests/detect/test_report.py
git commit -m "feat: emit a companion FASTA and CDS/translation features alongside detection-result GFF3"
```

---

### Task 3: Install and hand-trial both tools against real, newly-fixed output (investigation)

**Files:**
- Modify: `pixi.toml` (add `clinker`/`pygenomeviz` dependencies)
- No other files — this task resolves real open questions empirically, informing Tasks 4-5's exact CLI design; it does not write MATPredict application code itself.

- [ ] **Step 1: Add the dependencies and install**

```toml
# pixi.toml [dependencies], alongside the existing diamond/blast/miniprot/exonerate entries
clinker = "*"
pygenomeviz = "*"
```

```bash
pixi install
pixi run clinker --version
pixi run python -c "import pygenomeviz; print(pygenomeviz.__version__)"
```

- [ ] **Step 2: Generate real Task-1 output for 2-3 real curated records** (pick small, fast ones already exercised elsewhere in this session's work, e.g. the Onygenales `Coccidioides` records or the Teloschistales Xanthoria records whose real translated sequences were just derived) via `build_gff_for_record`, producing real `locus.gbk` files with real sequence and real CDS/translation features.

- [ ] **Step 3: Run a real clinker trial**

```bash
mkdir -p /tmp/clinker_trial
cp db/Ascomycota/Onygenales/*/locus.gbk /tmp/clinker_trial/  # or whichever 2-3 records you picked, renamed uniquely if filenames collide
pixi run clinker /tmp/clinker_trial/*.gbk -o /tmp/clinker_trial/output.html
```

Resolve the real open question from the research doc: does clinker accept these real `CDS`-with-`translation` GenBank files cleanly? Does it need anything Task 1's output doesn't yet provide (check clinker's actual error output, not just its docs, if it fails)? Document the real, concrete finding (works as-is / needs X additional feature or qualifier) — do not guess if the trial itself is ambiguous; read clinker's actual source (`pixi run python -c "import clinker; print(clinker.__file__)"` to locate it, then read its real parser code) if the CLI's own error messages aren't clear enough.

- [ ] **Step 4: Run a real pyGenomeViz trial** — both single-locus (one real curated record's `locus.gbk` or `locus.gff3`+companion FASTA) and, if feasible, its synteny mode against 2 real Task-2 detection-result GFF3+FASTA outputs (or 2 curated records). Resolve the research doc's second open question: does pyGenomeViz's GFF3 parser tolerate this project's fragmented-multi-contig-locus shape (`locus_group` attribute, per-segment `MAT_locus` parent features)? Use a real fragmented record if one exists in `db/` (check the genome-scale-detection-rollout findings/curated records for a real fragmented example — the Xanthoria `2903220_liq80xsp_MAT_combined` record is a real, already-curated fragmented-across-2-scaffolds example) as the real test case, not a synthetic one, since this is exactly the shape the research flagged as unconfirmed.

- [ ] **Step 5: Write up the real findings** to a new dated section appended to `docs/notes/2026-09-18_mat-locus-visualization-research.md` (do not create a new file — this is a direct follow-up to that research), resolving both of its own open questions with real, concrete evidence (exact error messages if something failed, exact command lines that worked, screenshots/description of real output if useful). This is what Tasks 4-5 will build against — if either tool needs an adaptation beyond what Tasks 1-2 already built, name it precisely here rather than discovering it mid-implementation in a later task.

- [ ] **Step 6: Commit**

```bash
git add pixi.toml pixi.lock docs/notes/2026-09-18_mat-locus-visualization-research.md
git commit -m "feat: add clinker/pygenomeviz dependencies; resolve real-data trial open questions"
```

---

### Task 4: `matpredict draw locus` — single-locus diagram via pyGenomeViz

**Files:**
- Create: `src/MATPredict/detect/draw.py` (or `src/MATPredict/db/draw.py` if Task 3's trial shows this is primarily a curated-record-side feature — decide based on Task 3's real findings about which input shape works better, and place it accordingly; document the choice)
- Modify: whichever `cli.py` registers the new subcommand (`db/cli.py` if this is a curation-workflow command like `matpredict curate-db draw-locus`, or `detect/cli.py` if it's `matpredict detect draw-locus` for a real detection result — again, base this on Task 3's real findings about the more natural real-world use case, likely BOTH: one entry point per input shape)
- Test: a new test file matching wherever the module lands

**Interfaces:** to be finalized once Task 3's real trial confirms pyGenomeViz's actual working API/CLI shape against this project's real GenBank/GFF3 output — do not guess pyGenomeViz's exact function signatures here; read its real, installed API (`pixi run python -c "import pygenomeviz; help(pygenomeviz)"` or its real source) as this task's own first step, informed by what Task 3 already confirmed works.

- [ ] **Step 1: Read Task 3's real findings** (the dated section it appended to the research doc) for the exact working command/API shape against this project's real output.

- [ ] **Step 2: Write a real, complete function** (`draw_locus(record: dict, gbk_path: Path, out_path: Path) -> Path`, or the GFF3+FASTA equivalent, whichever Task 3 confirmed works better) that calls pyGenomeViz's real API to produce one static image (SVG or PNG — pyGenomeViz's own default, per Task 3's trial) per curated record or detection result, with genes colored by `role` (a fixed, small color palette: e.g. `core_MAT` one color, `flanking_conserved` another, `flanking_variable` a third — check this project's `dataviz` skill/palette conventions if one applies to a scientific figure like this, otherwise pick a simple, colorblind-safe 3-color set and document the choice) and labeled with `gene_class`/idiomorph where present.

- [ ] **Step 3: Write a real test** using a real small curated record fixture (or the actual `tmp_path`-based synthetic fixture pattern this project's other tests already use) confirming a real output file is produced and is non-empty/well-formed (parseable as the expected image/SVG format — don't just check `.exists()`).

- [ ] **Step 4: Wire the CLI subcommand(s)**, following this project's existing `action.add_parser(...)`/`set_defaults(func=...)` registration style.

- [ ] **Step 5: Run the full suite, commit**

```bash
pixi run pytest -v
git add <files touched>
git commit -m "feat: add matpredict draw-locus command for single-locus gene-structure diagrams"
```

---

### Task 5: `matpredict draw synteny` — group/comparison diagram via clinker

**Files:**
- Create: `src/MATPredict/db/synteny.py` (or wherever fits best per Task 3's findings)
- Modify: whichever `cli.py` registers the new subcommand
- Test: a new test file matching wherever the module lands

**Interfaces:** same caveat as Task 4 — finalize the exact function signature against Task 3's real, confirmed clinker command shape, not a guess.

- [ ] **Step 1: Design the real color/function-label CSV generation**: a function that takes a list of record IDs (or a mix of curated record IDs and a real detection result's GFF3+FASTA) and produces clinker's real `--colour_map`/`--gene_functions` CSV inputs directly from each gene's `role`/`gene_class` schema fields — this is the concrete, real value MATPredict adds on top of raw clinker (auto-generating exactly the per-gene annotation clinker needs from data this project already curates, rather than a curator hand-building those CSVs).

- [ ] **Step 2: Write a real function** (`draw_synteny(record_ids: list[str], db_root: Path, out_path: Path) -> Path`, or the real shape Task 3's trial confirmed) that gathers each record's real `locus.gbk` (from Task 1's now-real output), generates the real colour-map/gene-functions CSVs from Step 1, and shells out to the real `clinker` CLI (reusing this project's existing subprocess-invocation conventions from `search.py`'s `_run_checked` pattern if one generalizes cleanly, or a similar direct `subprocess.run` call otherwise).

- [ ] **Step 3: Write a real test** using 2-3 real small curated record fixtures, confirming clinker is actually invoked with the right file list and the right generated CSVs (mock the subprocess call itself per this project's existing search.py-adjacent test conventions — read `tests/detect/test_search.py`'s subprocess-mocking pattern for the style to match), and confirming the generated CSVs' real content (role→color mapping, gene_class labels) matches the input records' real schema fields.

- [ ] **Step 4: Wire the CLI subcommand**, following this project's existing registration style.

- [ ] **Step 5: Run the full suite, commit**

```bash
pixi run pytest -v
git add <files touched>
git commit -m "feat: add matpredict draw-synteny command for multi-locus comparison diagrams via clinker"
```

---

### Task 6: Real end-to-end verification (data step)

**Files:** none new — this task runs Tasks 4-5's real tooling against real curated data and verifies genuine, correct output.

- [ ] **Step 1: Generate one real single-locus diagram** for a real, already-curated Onygenales or Xanthoria record (pick one with real, verified CDS/translation data from this session's earlier work) via `matpredict draw-locus`, and visually/structurally confirm it's correct: right number of genes, right colors by role, right relative gene order/orientation matching the record's own real coordinates.

- [ ] **Step 2: Generate one real synteny/group diagram** comparing 2-3 real, related curated records (e.g. the 2 real Coccidioides Onygenales records, or the 2 real Xanthoria records whose `MAT1-1-1` genes now have real, comparable translated sequences from this session's earlier curation work) via `matpredict draw-synteny`, and confirm the real identity links clinker draws between homologous genes look sensible (e.g. `MAT1-1-1` in one record should link to `MAT1-1-1` in the other with a real, non-trivial percent identity — cross-check against this session's own earlier finding that the 2 real Xanthoria `MAT1-1-1` sequences differ by exactly one residue, so clinker's own computed identity for that pair should come out very high, near 100%, as a sanity check that the whole pipeline is wired correctly end to end).

- [ ] **Step 3: Report the real results** — do not silently close this plan; present both real diagrams' locations and a short honest assessment (did they look correct? any surprises?) to the human reviewer, matching this project's established practice for every rollout/verification task.

## Self-review notes (controller, at plan-writing time)

- **Spec coverage:** the research doc's 2 concrete recommendations (pyGenomeViz for single-locus, clinker for group/synteny) map to Tasks 4-5; its identified writer gap (no real sequence/CDS features) maps to Tasks 1-2; its 2 unresolved open questions map to Task 3, deliberately sequenced BEFORE Tasks 4-5 commit to exact interfaces, so this plan does not guess at either tool's real behavior.
- **No placeholders in Tasks 1-2** (the parts with enough existing, read codebase context to write real, complete code): both are grounded in this session's actual current `gff_export.py`/`report.py` source, reusing already-built, already-reviewed machinery (`NcbiClient.fetch_nucleotide_sequence`, `benchmark.py`'s exon-splicing/translation logic) rather than inventing new mechanisms. Tasks 4-5 are deliberately left to finalize their exact interfaces against Task 3's real trial results rather than the plan guessing pyGenomeViz/clinker's precise real API/CLI shape sight-unseen — this is not a placeholder in the "TBD, figure it out later" sense the writing-plans skill forbids; it is a genuine external-tool-behavior dependency the plan's own Task 3 is structured to resolve with real evidence before Tasks 4-5 execute, the same pattern the very first plan of this session used for its own genuinely-unresolved "is `datasets` CLI available" question.
- **Type/signature consistency:** Task 1's `write_genbank`'s new `ncbi` parameter is additive and optional, preserving every existing caller's behavior (verified by an explicit backward-compatibility test). Task 2's `write_detection_gff3`'s new `genome_fasta` parameter follows the identical additive pattern.
