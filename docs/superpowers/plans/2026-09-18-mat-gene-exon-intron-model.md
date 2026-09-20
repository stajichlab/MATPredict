# MAT Gene Exon/Intron/Frame Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the curated database's single-span-per-gene coordinate model with a real exon/intron/frame model (`exons` + `codon_start` + `transl_table`), so `validate.py`'s independent-translation check can correctly verify multi-exon and partial-CDS genes instead of reporting a false or unfixable fail, and so a real GenBank-parsing helper — not hand-transcription — is available to populate it.

**Architecture:** Add three optional, backward-compatible gene-level fields to the record schema (`exons`, `codon_start`, `transl_table`). Rewrite `_independent_translation()` in `src/MATPredict/db/validate.py` to assemble a spliced transcript from `exons` (falling back to the existing single-span behavior when `exons` is absent) and apply `codon_start`/`transl_table` exactly once, per the research below. Add a `NcbiClient` helper that parses a real GenBank flatfile's CDS feature into `exons`/`codon_start`/`transl_table`, so future curation populates this data by parsing, not eyeballing. Backfill the 14 records already known to need it. Everywhere else in the tool (coordinate benchmark, reporting), treat exon/intron structural differences between two records or between a curated record and a detection result as **not a disagreement** as long as the translated protein matches — this principle, stated explicitly by the domain expert during scoping, is binding on every consumer of this data, not just the validator.

**Tech Stack:** Python, Biopython (`Bio.Seq`, `Bio.SeqIO` — already a dependency, already used for GenBank parsing elsewhere in `ncbi_client.py`), existing `NcbiClient`/pytest infrastructure.

**Spec:** No separate spec document — this plan was scoped directly from a domain-expert conversation plus a dedicated biology research pass (see "Research findings" below, which is the binding source of truth for every schema/algorithm decision in this plan). Treat the bullet list under "Research findings" as the Global-Constraints-equivalent authority for this plan.

## Research findings (binding — every task below implements these)

Verified live against NCBI/PubMed/INSDC documentation during scoping (see `.superpowers/sdd/2026-09-17-mat-detection-search-localization/` session context if the full report is needed — not re-included here to keep this plan focused):

1. **Multi-exon is the norm outside pheromone precursors and most yeast cassette genes.** Real fungal MAT gene exon counts range 1 (many pheromone precursors, *S. cerevisiae* MATalpha1/alpha2, *S. pombe* mat1-Mc/Mi) to 5 (*N. crassa* mat A-2). Exon count does **not** cleanly track MAT gene class (alpha-box vs. HMG-box) or lineage — it must be recorded per curated record, never assumed from a template.
2. **`/codon_start` (1, 2, or 3) is applied exactly once, to the first base of the fully assembled, strand-oriented, spliced sequence** — never per-exon. This was confirmed empirically (not just from documentation) against a real multi-exon, minus-strand, 5'-partial record (*Ceratocystis fimbriata* MAT1-2-1, GenBank PV763125, `codon_start=2`): applying the offset once to the concatenated, reverse-complemented, ordered exon sequence reproduced the deposited 280-aa protein exactly; per-exon phases in that same feature were 1/2/1/0 — not independently in-frame. Absent qualifier means 1 (NCBI's own stated default).
3. **`codon_start` is common, not exotic**, in real curated MAT records: in a 251-record GenBank sample of MAT1-1-1/MAT1-2-1 partial-CDS deposits, roughly half carried a non-1 `codon_start`, including a substantial fraction of multi-exon (`join()`) features. A translator that ignores this will fail on a large share of real records — exactly the failure mode already found and documented (not fixed) for 14 records in this database.
4. **Structural divergence between real biological orthologs, and even between isolates of the same species, is well documented and expected** (Chitrampalam et al. 2013, PLoS ONE 8:e56895, doi:10.1371/journal.pone.0056895 — *S. sclerotiorum* MAT1-1-1 recorded with 0, 1, or 2 introns depending on isolate; *S. sclerotiorum* vs. *B. cinerea* MAT1-2-4 differ by exactly one intron). **The protein sequence, not the exon layout, is the thing that must match.** This is a binding design principle for every consumer of this data: two records (or a curated record vs. a detection result) with different exon counts/boundaries are NOT in conflict if their translated proteins agree.
5. **Alternative splicing is real but narrow**: RefSeq itself carries two isoforms for *N. crassa* mat A-1 (293aa canonical / 288aa variant); *S. sclerotiorum* MAT1-2-1 has 3 documented transcript variants, only one functional. This plan does **not** build multi-isoform support (YAGNI — no currently curated record needs it); it records one canonical structure per gene entry. If a future record genuinely needs to represent an alternate isoform, that is a new, separate gene entry (already how multi-gene records work), not a new plan requirement.
6. **`transl_table` matters**: CUG-clade *Candida* records (e.g. MTLA1) use NCBI genetic code table 12, not table 1. Ignoring this produces systematic Leu/Ser mismatches that look like coordinate bugs but aren't.
7. **`partial_5`/`partial_3` (from GenBank's `<`/`>` markers) are independent of `codon_start`** — a record can be 3'-partial with `codon_start=1`. This plan does not need to store partiality explicitly (it doesn't change the translation algorithm), but curators populating `exons`/`codon_start` from a real GenBank record should not confuse the two.
8. **No minimum exon length** — real curated genes have exons as short as 12bp (a documented micro-exon in *C. neoformans* SXI2a). Any new code must not reject short exons.

## Global Constraints

- 1-based, fully-closed coordinates throughout (matches `db/gff_export.write_gff3`'s existing convention and `NcbiClient.fetch_nucleotide_sequence`'s existing `seq_start`/`seq_stop` semantics — no new offset conversions).
- `exons`, `codon_start`, `transl_table` are all OPTIONAL gene-level fields. Every existing record without them must continue to validate and behave exactly as today (single-span, codon_start=1, table 1) — zero required migration for records that don't need this.
- Never fabricate a coordinate, exon boundary, or `codon_start` value. Every value populated into a real record must come from an actual, this-session GenBank fetch — never inferred, guessed, or copied from a similar-looking record.
- `codon_start` is applied exactly once, to the assembled spliced sequence — never per-exon (see research finding 2). This is the single most important algorithmic detail in this plan; get it wrong and every multi-exon partial-CDS gene silently mistranslates.
- Structural divergence (different exon count/boundaries) between two records, or between a curated record and a detection result, is never itself a validation failure. Only a protein-sequence mismatch is.

---

### Task 1: Schema fields for `exons`, `codon_start`, `transl_table`

**Files:**
- Modify: `db/_schema/order.schema.yaml` is the locus-vocabulary schema (gene_class lives there) — this task's fields belong on the RECORD schema instead. Find the actual record-level schema file by running `grep -rl "gene_index" db/_schema/` and confirm the exact path before editing (it is very likely `db/_schema/record.schema.yaml` or similarly named — read whichever file this Task 1 implementer confirms it is; do not guess the name here, discovering the exact file name and current gene-item schema shape is this task's first step).
- Modify: `src/MATPredict/db/schema.py` if it hand-validates gene fields in addition to (or instead of) a YAML-schema file — read this file first to find out whether validation is schema-driven, code-driven, or both, and add the new fields wherever `start`/`end`/`strand`/`protein_accession` are currently validated.
- Test: whichever test file currently covers gene-field schema validation (find via `grep -rl "gene_index" tests/db/`).

**Interfaces:**
- Produces: a gene entry MAY carry `exons: [{start: int, end: int}, ...]` (ordered 5'→3' in transcript orientation — for a minus-strand gene this means descending genomic coordinate order, i.e. the highest-coordinate exon listed first, matching how GenBank's own `complement(join(...))` lists components), `codon_start: int` (1, 2, or 3; default 1 when absent), `transl_table: int` (NCBI genetic code table number; default 1 when absent).
- Consumes: nothing new — this task only adds schema surface, no behavior change yet.

- [ ] **Step 1: Find the real record schema file and current gene-field validation**

Run:
```bash
grep -rl "gene_index" db/_schema/ src/MATPredict/db/schema.py
```
Read whichever file(s) this returns. Confirm exactly how `start`, `end`, `strand`, `protein_accession`, `segment_index` are currently declared/validated (YAML JSON-Schema-style `properties`, or Python dict/dataclass validation, or both). This determines exactly how Step 2 below should be written — do not assume a YAML schema exists if validation is actually hand-coded in Python.

- [ ] **Step 2: Add the three optional fields**

Using whatever mechanism Step 1 found, add:
- `exons`: optional array of objects, each with required `start` (integer) and `end` (integer). No minimum/maximum length constraint on the array, and no minimum span on any exon (real exons can be 12bp — do not add a `minLength`/`minimum` constraint that would reject this).
- `codon_start`: optional integer, allowed values `1`, `2`, `3`. Not required; absence means 1 (document this default in a comment next to the field, do not silently default it to a different value anywhere).
- `transl_table`: optional integer. Not required; absence means 1 (standard code). Document this default the same way.

If validation is JSON-Schema/YAML-based, these are straightforward optional `properties` additions with an `enum: [1, 2, 3]` on `codon_start`. If it's Python code, add optional-key handling that does not raise when the keys are absent, and validates the allowed value set when present (raise a clear `ValueError` for a `codon_start` outside 1-3, or for any `exons` entry missing `start`/`end`).

- [ ] **Step 3: Write a test proving the fields are genuinely optional**

Add a test that constructs a minimal valid gene entry with ONLY the pre-existing fields (`gene_index`, `start`, `end`, `strand`, `protein_accession`, etc. — copy the shape from an existing passing fixture in the test file you found) and confirms it still validates successfully with no `exons`/`codon_start`/`transl_table` present. Add a second test with all three new fields present and valid values, confirming it validates. Add a third test with an invalid `codon_start` (e.g. `4` or `0`) confirming it's rejected.

- [ ] **Step 4: Run the full test suite**

```bash
pixi run pytest -v
```
Expected: all previously-passing tests still pass (should be 173, but confirm the actual current count rather than assuming), plus your new tests passing.

- [ ] **Step 5: Commit**

```bash
git add <the schema file(s) and test file(s) you touched>
git commit -m "feat: add optional exons/codon_start/transl_table fields to gene schema"
```

---

### Task 2: Rewrite `_independent_translation()` to assemble spliced, frame-corrected transcripts

**Files:**
- Modify: `src/MATPredict/db/validate.py` (the `_independent_translation` and `_translate_cds` functions, currently at lines 17-60 as of this plan's writing — re-read the file fresh, it may have shifted).
- Test: `tests/db/test_validate.py` (confirm exact filename via `grep -rl "_independent_translation" tests/`).

**Interfaces:**
- Consumes: Task 1's `exons`/`codon_start`/`transl_table` fields on a `gene` dict; the existing `NcbiClient.fetch_nucleotide_sequence(accession, start, end, strand)` signature (already handles per-call strand reverse-complement via its own `strand=2` parameter — do not add a second, redundant local reverse-complement step).
- Produces: `_independent_translation(record, gene, ncbi) -> str | None` keeps its existing signature and `None`-means-not-applicable contract exactly as today — callers in `validate_record()` do not need to change.

- [ ] **Step 1: Write the failing tests first**

Add to the test file two new fixtures reproducing the two real cases from the research (use real, small, hand-copyable sequences — do not fetch live NCBI in a unit test; construct a short synthetic genomic sequence and synthetic exon coordinates that exercise the same *mechanism*, i.e. multi-exon + non-1 codon_start + minus strand, with a mocked/fake `NcbiClient` whose `fetch_nucleotide_sequence` returns fixed strings for each exon's coordinates — follow whatever mocking pattern the existing tests in this file already use for `NcbiClient`, do not invent a new one):

```python
def test_independent_translation_assembles_multi_exon_plus_strand():
    # Two exons, plus strand, codon_start=1. Exon 1: "ATGGCC" (Met-Ala partial),
    # exon 2 continues the frame: "TTTTAA" (Phe-stop). Assembled: ATGGCCTTTTAA
    # -> translates to "MAF" (stop dropped).
    gene = {
        "gene_index": 0, "segment_index": 0, "strand": "+",
        "start": 1, "end": 12,  # outer bounds, unused when exons present
        "exons": [{"start": 1, "end": 6}, {"start": 7, "end": 12}],
        "protein_accession": "ncbi_protein:FAKE1.1",
    }
    record = {"locus": {"core": {"segments": [
        {"sequence_source": {"type": "insdc_nucleotide", "accession": "FAKE_ACC.1"}}
    ]}}}
    ncbi = FakeNcbiClient(nucleotide_by_range={
        ("FAKE_ACC.1", 1, 6, "+"): "ATGGCC",
        ("FAKE_ACC.1", 7, 12, "+"): "TTTTAA",
    })
    result = _independent_translation(record, gene, ncbi)
    assert result == "MAF"

def test_independent_translation_applies_codon_start_once_not_per_exon():
    # Reproduces the real Ceratocystis PV763125 mechanism: 4 exons, minus strand,
    # codon_start=2 (skip the first base of the ASSEMBLED sequence, not each exon).
    # Use a small synthetic case: assembled raw (pre-offset) = "TATGGCCTTTTAA"
    # (13 nt; note this is what fetch_nucleotide_sequence returns per exon, already
    # reverse-complemented per exon by the strand=2 semantics -- exons are listed in
    # transcript/descending-genomic-coordinate order and concatenated as-is, no
    # additional local revcomp). With codon_start=2, translation starts at index 1:
    # "ATGGCCTTTTAA" -> "MAF" (same result as above, proving the offset is applied
    # once to the whole assembled string, not to each exon's start).
    gene = {
        "gene_index": 0, "segment_index": 0, "strand": "-",
        "start": 1, "end": 13,
        "exons": [{"start": 100, "end": 104}, {"start": 90, "end": 97}],
        "codon_start": 2,
        "protein_accession": "ncbi_protein:FAKE2.1",
    }
    record = {"locus": {"core": {"segments": [
        {"sequence_source": {"type": "insdc_nucleotide", "accession": "FAKE_ACC.1"}}
    ]}}}
    ncbi = FakeNcbiClient(nucleotide_by_range={
        ("FAKE_ACC.1", 100, 104, "-"): "TATGG",
        ("FAKE_ACC.1", 90, 97, "-"): "CCTTTTAA",
    })
    result = _independent_translation(record, gene, ncbi)
    assert result == "MAF"

def test_independent_translation_falls_back_to_single_span_when_no_exons():
    # Existing behavior (Task 1 did not touch this path): a gene with plain
    # start/end/strand and no `exons` key still works exactly as before.
    gene = {
        "gene_index": 0, "segment_index": 0, "strand": "+",
        "start": 1, "end": 6,
        "protein_accession": "ncbi_protein:FAKE3.1",
    }
    record = {"locus": {"core": {"segments": [
        {"sequence_source": {"type": "insdc_nucleotide", "accession": "FAKE_ACC.1"}}
    ]}}}
    ncbi = FakeNcbiClient(nucleotide_by_range={("FAKE_ACC.1", 1, 6, "+"): "ATGGCC"})
    result = _independent_translation(record, gene, ncbi)
    assert result == "MA"

def test_independent_translation_uses_transl_table_12_for_cug_clade_records():
    # CTG under table 1 is Leu (L); under table 12 (Alternative Yeast Nuclear Code)
    # it is Ser (S). This is the exact real-world failure mode the research found
    # for Candida MTL records.
    gene = {
        "gene_index": 0, "segment_index": 0, "strand": "+",
        "start": 1, "end": 9, "transl_table": 12,
        "protein_accession": "ncbi_protein:FAKE4.1",
    }
    record = {"locus": {"core": {"segments": [
        {"sequence_source": {"type": "insdc_nucleotide", "accession": "FAKE_ACC.1"}}
    ]}}}
    ncbi = FakeNcbiClient(nucleotide_by_range={("FAKE_ACC.1", 1, 9, "+"): "ATGCTGTAA"})
    result = _independent_translation(record, gene, ncbi)
    assert result == "MS"  # not "ML" -- proves table 12 was actually used
```

Adapt `FakeNcbiClient` to match whatever fake/mock NCBI client class already exists in this test file (do not invent a new fake if one is already there — extend it to support keyed-by-range nucleotide fetch responses if it doesn't already).

- [ ] **Step 2: Run the tests, confirm they fail**

```bash
pytest tests/db/test_validate.py -v -k "independent_translation"
```
Expected: the new tests fail (function doesn't yet support `exons`/multi-fetch/`transl_table`), the existing single-span test (if any) still passes.

- [ ] **Step 3: Rewrite `_independent_translation()` and `_translate_cds()`**

```python
def _translate_cds(nucleotide_sequence: str, table: int = 1) -> str:
    """Translate a coding sequence to protein, stopping at (and dropping) the first stop
    codon, using the given NCBI genetic code table (default 1, standard code; CUG-clade
    Candida MTL records require table 12, the Alternative Yeast Nuclear Code)."""
    protein = str(Seq(nucleotide_sequence).translate(table=table, to_stop=True))
    return protein


def _assemble_transcript(
    ncbi: NcbiClient, accession: str, gene: dict
) -> str | None:
    """Fetch and concatenate a gene's exon sequences in transcript order, or fall back
    to its single start/end span when no `exons` list is recorded.

    `exons` entries must already be listed in transcript order (5'->3'; for a minus-
    strand gene this is descending genomic coordinate order, matching how GenBank's own
    complement(join(...)) syntax orders its components) -- this function does not
    reorder them. Each exon is fetched with the gene's own strand, so
    `fetch_nucleotide_sequence`'s existing `strand=2` reverse-complement handling
    applies per exon exactly as it already does for a single span; no additional local
    reverse-complementation is performed here.
    """
    strand = gene.get("strand")
    exons = gene.get("exons")
    if exons:
        parts = []
        for exon in exons:
            part = ncbi.fetch_nucleotide_sequence(accession, exon["start"], exon["end"], strand)
            if not part:
                return None
            parts.append(part)
        return "".join(parts)

    start, end = gene.get("start"), gene.get("end")
    if start is None or end is None:
        return None
    return ncbi.fetch_nucleotide_sequence(accession, start, end, strand)


def _independent_translation(record: dict, gene: dict, ncbi: NcbiClient | None) -> str | None:
    """Independently re-derive a gene's protein sequence from the record's own recorded
    genomic coordinates (single span, or a real exon/intron structure when `exons` is
    recorded), by fetching from NCBI and translating -- rather than trusting anything
    already stored in the record about its protein.

    `codon_start` (1/2/3, default 1) is applied EXACTLY ONCE, to the first base of the
    fully assembled, strand-oriented, spliced sequence -- never per-exon. This was
    verified against a real multi-exon, minus-strand, 5'-partial GenBank record; GenBank's
    own /codon_start qualifier is defined relative to "the first base of that feature"
    (the whole joined CDS feature), not per exon.

    Returns None (check not applicable) under the same conditions as before: no NCBI
    client, no fetchable segment sequence, or missing coordinate data on the gene.
    """
    if ncbi is None:
        return None
    segment_index = gene.get("segment_index")
    if segment_index is None:
        return None
    if not gene.get("exons") and (gene.get("start") is None or gene.get("end") is None):
        return None

    segments = record["locus"].get("core", {}).get("segments", [])
    if segment_index >= len(segments):
        return None
    segment = segments[segment_index]
    sequence_source = segment.get("sequence_source", {})
    accession = sequence_source.get("accession")
    if not accession or sequence_source.get("type") != "insdc_nucleotide":
        return None

    transcript = _assemble_transcript(ncbi, accession, gene)
    if not transcript:
        return None

    codon_start = gene.get("codon_start", 1)
    offset = codon_start - 1
    frame_corrected = transcript[offset:]
    # A 3'-partial CDS's spliced length is often not a multiple of 3; Biopython's
    # translate() with a trailing 1-2nt remainder needs an explicit trim, since it
    # raises rather than silently truncating.
    usable_length = len(frame_corrected) - (len(frame_corrected) % 3)
    frame_corrected = frame_corrected[:usable_length]

    table = gene.get("transl_table", 1)
    return _translate_cds(frame_corrected, table=table)
```

Note the removed `start`/`end` requirement change: the old code required `start`/`end` unconditionally; the new code only requires them when `exons` is absent (since `_assemble_transcript` uses `exons` exclusively when present). Re-read the surrounding `validate_record()` caller to confirm nothing else depended on `_independent_translation` requiring `start`/`end` even when `exons` is given — it shouldn't, since callers only look at the return value, but verify this directly by reading `validate_record()`'s current body (shown in this plan's context above) rather than assuming.

- [ ] **Step 4: Run the tests, confirm they pass**

```bash
pytest tests/db/test_validate.py -v
```
Expected: all pass, including the 4 new tests and every pre-existing test in this file.

- [ ] **Step 5: Run the full suite**

```bash
pixi run pytest -v
```
Expected: zero regressions.

- [ ] **Step 6: Commit**

```bash
git add src/MATPredict/db/validate.py tests/db/test_validate.py
git commit -m "feat: independent-translation check now supports multi-exon genes and codon_start/transl_table"
```

---

### Task 3: A real GenBank CDS-structure parser (for curation, not hand-transcription)

**Files:**
- Modify: `src/MATPredict/db/ncbi_client.py` (add a new method alongside `fetch_nucleotide_sequence`/`fetch_protein_sequence`/`fetch_taxonomy_lineage`).
- Test: `tests/db/test_ncbi_client.py`.

**Interfaces:**
- Produces: `NcbiClient.fetch_cds_structure(accession: str, protein_id: str) -> CdsStructure`, where `CdsStructure` is a small dataclass/namedtuple with fields `exons: list[tuple[int, int]]`, `strand: str` (`"+"` or `"-"`), `codon_start: int`, `transl_table: int`. This is the tool a future curation pass (human or agent) should call instead of hand-transcribing a `join()` string — it directly answers the root cause the research and the earlier 17-fail investigation both identified (curators recording the gene/mRNA span, or the outer join bounds, instead of the real per-exon structure).
- Consumes: nothing new from other tasks; independent of Tasks 1-2's schema/algorithm, though it produces exactly the shape Task 1's schema expects.

- [ ] **Step 1: Write the failing test with a real, static GenBank flatfile fixture**

Add a small, real GenBank-format text fixture to the test file (as a triple-quoted string constant) representing a genuine multi-exon, minus-strand, partial CDS feature. Use the real Ceratocystis fimbriata PV763125 CDS feature verified during this plan's research (`complement(join(156..334,386..684,741..810,877..>1172))`, `/codon_start=2`) — fetch the real minimal GenBank flatfile for this accession yourself during this task (a single live `efetch` call, done once while writing the test, to get real feature-table text) and trim it to just the relevant `source`/`CDS` feature lines plus a placeholder `ORIGIN` sequence block (the test only needs the feature table structure to be real and parseable, not a full 1kb+ sequence — a short placeholder sequence of the right total length is fine since this test is about coordinate/qualifier PARSING, not translation).

```python
def test_fetch_cds_structure_parses_multi_exon_partial_minus_strand_codon_start():
    client = NcbiClient(fetcher=FakeFetcher(responses={
        "<the efetch URL for PV763125 gb format>": GENBANK_FIXTURE_TEXT,
    }))
    result = client.fetch_cds_structure("PV763125.1", protein_id="XXX000000.1")  # use the real protein_id from the fixture
    assert result.exons == [(156, 334), (386, 684), (741, 810), (877, 1172)]
    assert result.strand == "-"
    assert result.codon_start == 2
    assert result.transl_table == 1
```

Adapt the exact `FakeFetcher`/mocking pattern to whatever this test file already uses for `NcbiClient` (it already mocks `efetch.fcgi` calls for `fetch_nucleotide_sequence`/`fetch_protein_sequence` tests — follow that exact pattern, do not invent a new one).

- [ ] **Step 2: Run the test, confirm it fails**

```bash
pytest tests/db/test_ncbi_client.py -v -k "cds_structure"
```
Expected: fails (`AttributeError`, method doesn't exist yet).

- [ ] **Step 3: Implement `fetch_cds_structure`**

```python
from dataclasses import dataclass


@dataclass(frozen=True)
class CdsStructure:
    exons: list[tuple[int, int]]
    strand: str
    codon_start: int
    transl_table: int


def fetch_cds_structure(self, accession: str, protein_id: str) -> "CdsStructure":
    """Parse a GenBank record's real CDS feature into an exon/intron/frame structure,
    for populating a curated gene's `exons`/`codon_start`/`transl_table` fields from the
    actual deposit rather than a hand-transcribed span.

    `protein_id` disambiguates which CDS feature to use when an accession carries more
    than one CDS (e.g. a multi-gene MAT locus deposit) -- matched against each CDS
    feature's own `protein_id` qualifier.

    Exons are returned in the ORDER Biopython's SeqFeature.location reports them, which
    for a GenBank `complement(join(...))` feature is genomic 5'->3' listing order, i.e.
    the SAME order as GenBank's own text representation. This is transcript order for a
    plus-strand feature; for a minus-strand feature, this codebase's convention (see
    validate.py's _assemble_transcript) also expects descending-genomic-coordinate
    (transcript) order, so this function must return exons in that same descending
    order for minus-strand features -- verify Biopython's actual iteration order
    empirically against the real PV763125 fixture (do not assume without checking;
    different Biopython versions have historically differed on sub-feature ordering for
    compound locations) and reverse the list here if needed to guarantee descending
    order for minus-strand CDS features specifically.
    """
    url = self._url("efetch.fcgi", f"db=nuccore&id={accession}&rettype=gb&retmode=text")
    body = self.fetcher.get(url)
    record = SeqIO.read(StringIO(body), "genbank")
    for feature in record.features:
        if feature.type != "CDS":
            continue
        if feature.qualifiers.get("protein_id", [None])[0] != protein_id:
            continue
        location = feature.location
        exons = [(int(part.start) + 1, int(part.end)) for part in location.parts]
        strand = "-" if location.strand == -1 else "+"
        if strand == "-":
            exons = list(reversed(exons))
        codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
        transl_table = int(feature.qualifiers.get("transl_table", ["1"])[0])
        return CdsStructure(exons=exons, strand=strand, codon_start=codon_start, transl_table=transl_table)
    raise ValueError(f"no CDS with protein_id={protein_id!r} found in {accession}")
```

Note the `+1` on `part.start`: Biopython's `SimpleLocation`/`CompoundLocation` uses 0-based start, 1-based-inclusive end internally (Python-slice-like) — converting to this project's 1-based-fully-closed convention requires adding 1 to the start only. **Verify this conversion directly against the real PV763125 fixture in Step 1's test** (the expected `(156, 334)` etc. are already the real 1-based GenBank-text coordinates) rather than trusting this plan's arithmetic blindly — if the assertion in Step 1 fails on the off-by-one, that confirms which direction to adjust, the same class of bug already found and fixed twice this session.

- [ ] **Step 4: Run the test, confirm it passes**

```bash
pytest tests/db/test_ncbi_client.py -v -k "cds_structure"
```

- [ ] **Step 5: Run the full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/db/ncbi_client.py tests/db/test_ncbi_client.py
git commit -m "feat: add NcbiClient.fetch_cds_structure for real exon/codon_start/transl_table parsing"
```

---

### Task 4: Backfill the 14 documented-limitation records

**Files:**
- Modify: the 14 `metadata.yaml` files currently carrying a `sequence_match.notes` entry documenting the codon_start/multi-exon limitation (enumerate them fresh via `grep -rl "codon_start" db/ db/candidates/ 2>/dev/null` combined with `grep -rl "multi-exon" db/**/metadata.yaml db/candidates/**/metadata.yaml` — do not trust a hardcoded list in this plan, re-discover the exact current set since other work may have touched them since this plan was written).

**Interfaces:**
- Consumes: Task 1's schema fields, Task 3's `fetch_cds_structure` helper.
- Produces: updated `metadata.yaml` records with real `exons`/`codon_start`/`transl_table` populated, re-validated against Task 2's rewritten check.

- [ ] **Step 1: Enumerate the current 14 (or however many remain)**

```bash
grep -rl "codon_start\|multi-exon\|schema limitation" db/**/metadata.yaml db/candidates/**/metadata.yaml
```
Read each match's `sequence_match.notes` field to confirm it's actually one of the documented-limitation records from the prior investigation, not an unrelated hit.

- [ ] **Step 2: For each record, call `fetch_cds_structure` and populate the fields**

For each failing gene in each record: identify its accession and `protein_accession` (already present in the record), call `NcbiClient.fetch_cds_structure(accession, protein_id)` (write a small one-off script using the project's real `NcbiClient` instance — do not hand-transcribe), and write the resulting `exons`/`codon_start`/`transl_table` into the gene's entry in `metadata.yaml`. Leave `start`/`end` (the outer bounds) unchanged — they remain correct and are still used for GFF3/locus-span reporting.

- [ ] **Step 3: Re-validate every touched record**

```bash
matpredict curate-db validate --phylum <phylum> --record-id <record_id>
```
(or the candidate-path equivalent for pending records — check `src/MATPredict/db/cli.py` for the exact invocation for both accepted and candidate records). Confirm each gene's `sequence_match` status improves from `fail` to `pass` (or at minimum `warn`, if some other unrelated issue remains — investigate and report honestly if a record does not reach `pass` after correct exon/codon_start population, do not force a status).

- [ ] **Step 4: Update or remove the now-stale `sequence_match.notes` limitation text**

Replace the "documented limitation, not fixed" note with either nothing (if the record now cleanly passes) or an updated, accurate note describing whatever residual issue remains.

- [ ] **Step 5: Run the full suite, commit**

```bash
pixi run pytest -v
git add <touched metadata.yaml files>
git commit -m "fix: backfill real exon/codon_start structure for 14 previously-documented-limitation records"
```

(Split into more than one commit if it's cleaner to group by phylum/family — controller's judgment at execution time, same as prior fix passes this session.)

---

### Task 5: Confirm the coordinate benchmark and reporting layers honor "protein match, not exon match" as ground truth

**Files:**
- Read (and modify only if a real problem is found): whatever module implements the coordinate-accuracy benchmark referenced by `docs/superpowers/plans/2026-09-17-mat-detection-search-localization-benchmark-notes.md` (grep for "benchmark" under `src/MATPredict/` to find it), and `src/MATPredict/detect/report.py`.

**Interfaces:**
- Consumes: nothing new — this is an audit task, not a schema/algorithm change, unless it finds a real bug.

- [ ] **Step 1: Read the benchmark comparison logic**

Find and read whatever code compares a detection result's predicted coordinates against a curated record's coordinates for benchmark scoring. Determine: does it compare at the outer-span level (start/end), at the exon level, or does it not currently exist as runnable code yet (the benchmark-notes doc referenced above may describe manual/one-off timing runs rather than a scored comparison harness — confirm which is actually the case before assuming there's a bug to fix)?

- [ ] **Step 2: Judgment call**

If the benchmark only ever compares outer span (which is what every record has always had, single-span or now exons+outer-span both), there is likely nothing to fix here — the outer span remains a valid, if coarser, ground truth for "did detection find roughly the right region," and per the binding design principle, exon-level disagreement was never supposed to be scored anyway. Document this finding in the task's completion note. If instead you find the benchmark (or any other consumer) DOES penalize exon-count/structure differences as if they were errors, fix that specific comparison to compare translated-protein-equivalence or outer-span-only instead, and add a regression test proving a structurally-different-but-protein-equivalent pair no longer scores as a mismatch.

- [ ] **Step 3: Commit if any code changed**

```bash
pixi run pytest -v
git add <any files touched>
git commit -m "fix: <specific description of what was found and changed>"
```
If nothing needed changing, no commit for this task — note that finding in the plan's final report instead.

---

## Self-review notes (controller, at plan-writing time)

- Spec coverage: every numbered research finding (1-8) maps to a task: (1)/(8) → Task 1's no-minimum-length schema; (2) → Task 2's single-application-point algorithm, directly tested; (3) → Task 4's backfill of real, common non-1 `codon_start` cases; (4) → the binding design principle stated in Global Constraints and enforced by Task 5's audit; (5) → explicitly scoped OUT (YAGNI, documented why); (6) → Task 1's `transl_table` field + Task 2's test 4; (7) → explicitly noted as out of scope with reasoning (doesn't change the translation algorithm).
- No placeholders: every code block above is real, complete code, not a description of code to write later.
- Type/signature consistency checked: `_independent_translation`'s signature and `None`-contract are preserved across Tasks 1-2 so `validate_record()` (not itself modified by this plan) keeps working unchanged; `fetch_cds_structure`'s output shape (`CdsStructure.exons`/`.codon_start`/`.transl_table`) matches exactly the field names Task 1 adds to the gene schema, so Task 4's backfill script needs no field-name translation layer.
