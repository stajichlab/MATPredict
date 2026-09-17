# Search Stage Revision: Localize-then-Polish Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `matpredict detect`'s unrestricted whole-genome `exonerate protein2genome` fallback with a `tblastn` genome-wide localization pass, followed by a two-tool cross-validated polish (`miniprot` + `exonerate --refine region`, both restricted to the small localized window), reporting per-gene tool-agreement status without feeding it into confidence tiering in v1.

**Architecture:** New `search.py` functions for `tblastn` localization (batched across all routed families) and `miniprot` polishing; extend `search.py`'s existing exonerate wrapper to support `--refine region` against a sliced window (reusing `_extract_window`); a new `polish.py` module that classifies each gene's two-tool outcome into one of four statuses (`polished_agree` / `polished_disagree` / `polished_single` / `unpolished`) and picks canonical coordinates; `pipeline.py`'s `run_pipeline` gains an explicit Stage-0 split between the genome-only path (needs fresh localization) and the fast-path per-family rescue (already has a cluster to anchor a window on, skips localization); `tiering.assign_tier`'s `second_pass_used` parameter is retired in favor of an `unpolished` flag with the identical Medium-tier-cap effect; `GeneEvidence`/`DetectionResult`/`report.py` are extended to carry per-gene status and both tools' models when they disagree.

**Tech Stack:** Python 3.11+, `blast` (bioconda, provides `tblastn`) and `miniprot` (bioconda) as new dependencies alongside the existing `diamond`/`exonerate`, pytest with injected `runner` callables (no unit test invokes a live binary, matching every existing wrapper in this project).

**Spec:** `docs/superpowers/specs/2026-09-17-mat-detection-search-localization-design.md` (revises the "Search" and "Boundary calling" sections of `docs/superpowers/specs/2026-09-16-mat-detection-pipeline-design.md`)

## Global Constraints

- No test invokes a live external binary (`tblastn`, `miniprot`, `exonerate`) — every wrapper takes an injected `runner: Callable`, defaulting to `subprocess.run`, exactly like every existing wrapper in `search.py`.
- Coordinates are 1-based, fully-closed everywhere.
- Attribution (which `(phylum, locus_name)` family a hit belongs to) is keyed on `record_id` via `family_registry.load_record_families`/`_attribute` — never on the bare gene name. Every new wrapper must reuse `_attribute`/`_roles_by_family`, not reimplement gene-name-based lookup.
- Minus-strand HSPs from any new tool must be normalized to `start <= end` with `strand` derived correctly — this project has already hit this exact bug class once (window-offset re-basing for exonerate) and once for a related attribution bug (Task 6's cross-family gene-name pooling); do not repeat either shape of mistake with the new tools.
- `tblastn`'s query set for a genome-only run is batched once across every routed family's genes (both `core_MAT` and `flanking_conserved`) — never one `tblastn` invocation per family.
- Tool agreement between `miniprot` and `exonerate --refine` is reported per gene but must NOT be consulted by `tiering.assign_tier` in this revision — `polished_agree` and `polished_disagree` must produce identical tiering outcomes.
- The fast-path per-family rescue (existing cluster, one missing core gene) skips `tblastn` localization entirely and polishes directly against a window padded around the existing cluster's span.

---

## Task 1: Add `blast` and `miniprot` dependencies

**Files:**
- Modify: `pixi.toml`

**Interfaces:**
- Produces: `tblastn` (from the `blast` bioconda package) and `miniprot` available on PATH in the pixi `default`/`test` environments.

- [ ] **Step 1: Add the dependencies**

In `pixi.toml`'s `[dependencies]` table, add alongside the existing `diamond`/`exonerate`/`taxonkit` entries:

```toml
blast = "*"
miniprot = "*"
```

Update the existing comment above the search-tool dependencies to mention the new tools and what they're for:

```toml
# Search binaries `matpredict detect` shells out to (see detect/search.py):
# diamond for the proteome fast path; tblastn (blast) for genome-wide
# localization; miniprot and exonerate (--refine region, on a sliced
# window) for precision polishing within a localized candidate region.
```

- [ ] **Step 2: Verify the environment resolves and the binaries are present**

Run: `pixi install` (or the project's equivalent env-sync command), then:
`pixi run which tblastn miniprot` (or `.pixi/envs/default/bin/tblastn --version` / `.pixi/envs/default/bin/miniprot --version` if `pixi run which` isn't available in this environment).
Expected: both binaries resolve to a path inside the pixi environment.

- [ ] **Step 3: Commit**

```bash
git add pixi.toml
git commit -m "$(cat <<'EOF'
build: add blast (tblastn) and miniprot as detect-pipeline dependencies

Needed for the localize-then-polish search revision: tblastn for
genome-wide localization, miniprot for precision polishing alongside
exonerate --refine region.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: `tblastn` localization wrapper

**Files:**
- Modify: `src/MATPredict/detect/search.py`
- Test: `tests/detect/test_search.py`

**Interfaces:**
- Consumes: `Family`/`FamilyKey` (`family_registry.py`), `_attribute`/`_roles_by_family`/`_run_checked`/`SearchToolError` (existing, in `search.py`), `SearchHit` (existing).
- Produces:
  ```python
  METHOD_TBLASTN = "tblastn_genome"

  def search_localize(
      genome_fasta: Path,
      families: list[Family],
      reference_fasta: Path,
      record_families: dict[str, FamilyKey],
      runner: Callable = subprocess.run,
  ) -> list[SearchHit]:
      """Genome-wide tblastn localization, batched across every routed
      family's genes (core_MAT and flanking_conserved). One makeblastdb
      build plus one tblastn invocation covers every family -- never one
      call per family. Returns SearchHit(method="tblastn_genome", ...)."""
  ```
  Low-complexity filtering is disabled (`-seg no`) so short, simple
  pheromone-precursor queries are not suppressed (spec section on
  short-ORF genes, restated for tblastn).

- [ ] **Step 1: Write the failing test**

Add to `tests/detect/test_search.py`:

```python
from MATPredict.detect.search import METHOD_TBLASTN, search_localize

# tblastn's query is the curated reference protein set, its database is the
# genome -- the OPPOSITE role assignment from diamond's fast path (where the
# predicted proteome is the query). So qseqid carries the reference header
# (record_id|geneN|gene_name) and sseqid is the genome's own contig name.
TBLASTN_TSV = (
    # qseqid                sseqid  pident length sstart send sframe
    "rec1|gene0|mfa1\tcontigA\t95.0\t40\t400\t100\t-1\n"  # minus strand: sstart > send
    "rec1|gene1|pra1\tcontigA\t90.0\t300\t3600\t4914\t1\n"  # plus strand
)


def fake_tblastn_runner(cmd, **kwargs):
    assert "-seg" in cmd and cmd[cmd.index("-seg") + 1] == "no"
    class Result:
        returncode = 0
        stdout = TBLASTN_TSV
        stderr = ""
    return Result()


def test_search_localize_normalizes_minus_strand_and_batches_one_call(tmp_path):
    calls = []

    def counting_runner(cmd, **kwargs):
        calls.append(cmd)
        return fake_tblastn_runner(cmd, **kwargs)

    hits = search_localize(
        genome_fasta=tmp_path / "genome.fa",
        families=[FAMILY],  # reuse the existing FAMILY fixture (aLocus, mfa1/pra1)
        reference_fasta=tmp_path / "reference.faa",
        record_families={"rec1": FAMILY.key},
        runner=counting_runner,
    )
    by_gene = {h.gene_name: h for h in hits}
    assert by_gene["mfa1"].start == 100 and by_gene["mfa1"].end == 400 and by_gene["mfa1"].strand == "-"
    assert by_gene["pra1"].start == 3600 and by_gene["pra1"].end == 4914 and by_gene["pra1"].strand == "+"
    assert all(h.method == METHOD_TBLASTN for h in hits)
    # exactly one tblastn invocation regardless of how many families/genes
    tblastn_calls = [c for c in calls if c[0] == "tblastn"]
    assert len(tblastn_calls) == 1
```

(Adjust the exact `-outfmt` column list and fake TSV shape to whatever real
column layout you choose to request — the test's binding contract is the
*resulting* `SearchHit` values, not a specific outfmt string, matching
this project's established pattern from the original diamond/exonerate
wrappers.)

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/detect/test_search.py -v -k search_localize`
Expected: FAIL (`ImportError`/`AttributeError`, `search_localize` doesn't exist yet).

- [ ] **Step 3: Implement**

Add to `search.py` (near the other search functions), reusing
`_roles_by_family`, `_attribute`, `_parse_reference_header`, and
`_run_checked` unchanged:

```python
METHOD_TBLASTN = "tblastn_genome"

# qseqid is the reference protein header (parsed via _parse_reference_header);
# sseqid is the genome's own contig name -- tblastn's query/db roles are the
# OPPOSITE of diamond's fast path (search_fast_path's query is the predicted
# proteome; here the query is the curated reference set and the genome is
# the database), so do not copy search_fast_path's qseqid/sseqid roles.
_TBLASTN_OUTFMT = "6 qseqid sseqid pident length sstart send sframe"


def search_localize(
    genome_fasta: Path,
    families: list[Family],
    reference_fasta: Path,
    record_families: dict[str, FamilyKey],
    runner: Callable = subprocess.run,
) -> list[SearchHit]:
    """Genome-wide tblastn localization -- one blastdb build and one
    tblastn call cover every routed family's genes (core_MAT and
    flanking_conserved), never one call per family. Coordinates are
    approximate (no splice awareness); Stage 2 polishing refines them.
    """
    roles_by_family = _roles_by_family(families)

    with tempfile.TemporaryDirectory() as tmp_dir_name:
        db_prefix = Path(tmp_dir_name) / "genome_db"
        _run_checked(runner, [
            "makeblastdb", "-in", str(genome_fasta), "-dbtype", "nucl", "-out", str(db_prefix),
        ])
        cmd = [
            "tblastn", "-query", str(reference_fasta), "-db", str(db_prefix),
            "-seg", "no", "-outfmt", _TBLASTN_OUTFMT,
        ]
        result = _run_checked(runner, cmd)

        hits: list[SearchHit] = []
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            qseqid, contig, pident, _length, sstart, send, sframe = line.split("\t")
            record_id, gene_name = _parse_reference_header(qseqid)
            attribution = _attribute(record_id, gene_name, record_families, roles_by_family)
            if attribution is None:
                continue
            family_key, role = attribution
            start, end = sorted((int(sstart), int(send)))
            strand = "+" if int(sframe) > 0 else "-"
            hits.append(SearchHit(
                family_key=family_key, gene_name=gene_name, role=role,
                contig=contig, start=start, end=end, strand=strand,
                identity=float(pident), reference_record_id=record_id,
                method=METHOD_TBLASTN, coverage=None,
            ))
        return hits
```

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/detect/test_search.py -v -k search_localize`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/detect/search.py tests/detect/test_search.py
git commit -m "$(cat <<'EOF'
feat: add tblastn genome-wide localization wrapper

Replaces exonerate's role as a genome-wide search tool. One
makeblastdb + one tblastn call covers every routed family's genes
(core_MAT and flanking_conserved) in a single invocation, with -seg no
so short pheromone-precursor queries aren't filtered as low-complexity.
Minus-strand HSPs are normalized (start<=end, strand from sframe).

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: `PolishModel` — shared exon-aware gene-model type

**Files:**
- Create: `src/MATPredict/detect/polish.py`
- Test: `tests/detect/test_polish.py`

**Interfaces:**
- Produces:
  ```python
  @dataclass(frozen=True)
  class ExonSpan:
      start: int  # 1-based, inclusive, genome coordinates
      end: int

  @dataclass(frozen=True)
  class PolishModel:
      gene_name: str
      family_key: FamilyKey
      role: str
      contig: str
      start: int
      end: int
      strand: str
      exons: list[ExonSpan]
      identity: float
      reference_record_id: str
      method: str  # "miniprot" | "exonerate_refine"

  def boundaries_agree(a: PolishModel, b: PolishModel, tolerance_bp: int = 10) -> bool:
      """True when a and b describe the same gene's exon structure within
      tolerance_bp at every boundary (same exon count, each corresponding
      exon's start/end within tolerance). Different exon counts never agree."""
  ```

This is a pure data/comparison module — no tool invocation here (Tasks
4/5 build the two tools' wrappers that *produce* `PolishModel`s).

- [ ] **Step 1: Write the failing test**

```python
# tests/detect/test_polish.py
from __future__ import annotations

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.polish import ExonSpan, PolishModel, boundaries_agree

KEY = FamilyKey("Basidiomycota", "aLocus")


def _model(exons, method="miniprot"):
    return PolishModel(
        gene_name="mfa1", family_key=KEY, role="core_MAT", contig="c1",
        start=exons[0].start, end=exons[-1].end, strand="+", exons=exons,
        identity=95.0, reference_record_id="rec1", method=method,
    )


def test_boundaries_agree_within_tolerance():
    a = _model([ExonSpan(100, 200), ExonSpan(250, 300)], method="miniprot")
    b = _model([ExonSpan(103, 198), ExonSpan(252, 297)], method="exonerate_refine")
    assert boundaries_agree(a, b, tolerance_bp=10) is True


def test_boundaries_disagree_beyond_tolerance():
    a = _model([ExonSpan(100, 200)], method="miniprot")
    b = _model([ExonSpan(150, 260)], method="exonerate_refine")
    assert boundaries_agree(a, b, tolerance_bp=10) is False


def test_boundaries_disagree_on_different_exon_count():
    a = _model([ExonSpan(100, 300)], method="miniprot")
    b = _model([ExonSpan(100, 200), ExonSpan(250, 300)], method="exonerate_refine")
    assert boundaries_agree(a, b) is False
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/detect/test_polish.py -v`
Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Implement**

```python
# src/MATPredict/detect/polish.py
"""Cross-tool gene-model comparison for the localize-then-polish search
revision -- see docs/superpowers/specs/2026-09-17-mat-detection-search-localization-design.md."""
from __future__ import annotations

from dataclasses import dataclass

from MATPredict.detect.family_registry import FamilyKey


@dataclass(frozen=True)
class ExonSpan:
    start: int
    end: int


@dataclass(frozen=True)
class PolishModel:
    gene_name: str
    family_key: FamilyKey
    role: str
    contig: str
    start: int
    end: int
    strand: str
    exons: list[ExonSpan]
    identity: float
    reference_record_id: str
    method: str


def boundaries_agree(a: PolishModel, b: PolishModel, tolerance_bp: int = 10) -> bool:
    """Same exon count, and each corresponding exon's start/end within tolerance_bp."""
    if len(a.exons) != len(b.exons):
        return False
    return all(
        abs(ea.start - eb.start) <= tolerance_bp and abs(ea.end - eb.end) <= tolerance_bp
        for ea, eb in zip(a.exons, b.exons)
    )
```

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/detect/test_polish.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/detect/polish.py tests/detect/test_polish.py
git commit -m "$(cat <<'EOF'
feat: add PolishModel and exon-boundary agreement comparison

Pure data/comparison layer for the localize-then-polish revision --
no tool invocation here. Agreement requires identical exon count and
every corresponding exon's boundaries within tolerance; different
exon counts never agree regardless of overall span overlap.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: `exonerate --refine region` polishing wrapper (exon-aware)

**Files:**
- Modify: `src/MATPredict/detect/search.py`
- Test: `tests/detect/test_search.py`

**Interfaces:**
- Consumes: `_extract_window` (existing, unchanged), `PolishModel`/`ExonSpan` (Task 3).
- Produces:
  ```python
  def polish_with_exonerate(
      genome_fasta: Path,
      family: Family,
      gene_name: str,
      reference_fasta: Path,
      record_families: dict[str, FamilyKey],
      window: tuple[str, int, int],  # (contig, start, end), 1-based inclusive, already padded
      runner: Callable = subprocess.run,
  ) -> PolishModel | None:
      """Run exonerate --model protein2genome --refine region against the
      sliced window (via _extract_window, offset re-based exactly like the
      existing windowed second pass), parsing per-exon coordinates from the
      GFF `exon` feature lines (not just the `gene` line) so PolishModel's
      exons reflect the real intron/exon structure. Returns None when no
      usable model is produced for this gene in this window."""
  ```

The existing `search_genomic` function's gene-level GFF parsing (line-level
`\tgene\t` filtering) does not currently capture per-exon structure — this
task adds that capture as a **new** function rather than modifying
`search_genomic`'s existing behavior, since `search_genomic` remains used
elsewhere (Task 8 removes its callers from `pipeline.py`, but the function
itself, and its existing tests, stay intact until Task 10's cleanup pass
confirms nothing else needs it — do not delete `search_genomic` in this
task).

- [ ] **Step 1: Write the failing test**

```python
# tests/detect/test_search.py (add)
from MATPredict.detect.polish import ExonSpan
from MATPredict.detect.search import polish_with_exonerate

EXONERATE_REFINE_GFF = (
    "c1\texonerate\tgene\t1\t400\t.\t+\t.\t"
    "gene_id 1 ; sequence rec1|gene0|mfa1 ; gene_orientation . ; identity 95.00 ; similarity 96.00\n"
    "c1\texonerate\texon\t1\t150\t.\t+\t.\tinsertions 0 ; deletions 0\n"
    "c1\texonerate\texon\t200\t400\t.\t+\t.\tinsertions 0 ; deletions 0\n"
)


def fake_exonerate_refine_runner(cmd, **kwargs):
    assert "--refine" in cmd and cmd[cmd.index("--refine") + 1] == "region"
    class Result:
        returncode = 0
        stdout = EXONERATE_REFINE_GFF
        stderr = ""
    return Result()


def test_polish_with_exonerate_parses_exon_structure(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")
    model = polish_with_exonerate(
        genome_fasta=tmp_path / "genome.fa",
        family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa",
        record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500),
        runner=fake_exonerate_refine_runner,
    )
    assert model.exons == [ExonSpan(1, 150), ExonSpan(200, 400)]
    assert model.identity == 95.0
    assert model.method == "exonerate_refine"


def test_polish_with_exonerate_returns_none_when_no_model(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")

    def empty_runner(cmd, **kwargs):
        class Result:
            returncode = 0
            stdout = ""
            stderr = ""
        return Result()

    model = polish_with_exonerate(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa", record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500), runner=empty_runner,
    )
    assert model is None
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/detect/test_search.py -v -k polish_with_exonerate`
Expected: FAIL (`ImportError`).

- [ ] **Step 3: Implement**

Add to `search.py`, reusing `_extract_window`'s existing offset-rebasing
pattern:

```python
from MATPredict.detect.polish import ExonSpan, PolishModel


def polish_with_exonerate(
    genome_fasta: Path,
    family: Family,
    gene_name: str,
    reference_fasta: Path,
    record_families: dict[str, FamilyKey],
    window: tuple[str, int, int],
    runner: Callable = subprocess.run,
) -> "PolishModel | None":
    roles_by_family = _roles_by_family([family])
    contig, win_start, _win_end = window

    with tempfile.TemporaryDirectory() as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        target_fasta = _extract_window(genome_fasta, window, tmp_dir)
        offset = win_start - 1

        cmd = [
            "exonerate", "--model", "protein2genome",
            "--query", str(reference_fasta), "--target", str(target_fasta),
            "--refine", "region", "--showtargetgff", "yes", "--showalignment", "no",
        ]
        result = _run_checked(runner, cmd)

        gene_line = None
        exon_lines = []
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            if fields[2] == "gene" and gene_line is None:
                gene_line = fields
            elif fields[2] == "exon":
                exon_lines.append(fields)

        if gene_line is None:
            return None

        attrs = gene_line[8].split(" ; ")
        query_id = next(p.split(" ")[1] for p in attrs if p.startswith("sequence "))
        record_id, matched_gene = _parse_reference_header(query_id)
        if matched_gene != gene_name:
            return None
        attribution = _attribute(record_id, matched_gene, record_families, roles_by_family)
        if attribution is None:
            return None
        family_key, role = attribution
        identity = 0.0
        for part in attrs:
            if part.startswith("identity "):
                identity = float(part.split(" ")[1])
                break

        exons = [
            ExonSpan(int(e[3]) + offset, int(e[4]) + offset)
            for e in sorted(exon_lines, key=lambda e: int(e[3]))
        ]
        return PolishModel(
            gene_name=matched_gene, family_key=family_key, role=role, contig=contig,
            start=int(gene_line[3]) + offset, end=int(gene_line[4]) + offset,
            strand=gene_line[6], exons=exons, identity=identity,
            reference_record_id=record_id, method="exonerate_refine",
        )
```

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/detect/test_search.py -v -k polish_with_exonerate`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/detect/search.py tests/detect/test_search.py
git commit -m "$(cat <<'EOF'
feat: add exonerate --refine region polishing wrapper with exon parsing

Runs against a sliced, offset-rebased window (reusing _extract_window
unchanged) and captures per-exon GFF lines, not just the gene-level
span, so PolishModel carries real intron/exon structure for
cross-tool agreement comparison. Returns None rather than a
zero-identity placeholder when no usable model is produced.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: `miniprot` polishing wrapper

**Files:**
- Modify: `src/MATPredict/detect/search.py`
- Test: `tests/detect/test_search.py`

**Interfaces:**
- Produces:
  ```python
  def polish_with_miniprot(
      genome_fasta: Path,
      family: Family,
      gene_name: str,
      reference_fasta: Path,
      record_families: dict[str, FamilyKey],
      window: tuple[str, int, int],
      runner: Callable = subprocess.run,
  ) -> PolishModel | None:
      """Run miniprot against the same sliced, offset-rebased window as
      polish_with_exonerate (reusing _extract_window), parsing miniprot's
      GFF/PAF output into a PolishModel with real exon structure."""
  ```

- [ ] **Step 1: Write the failing test**

Follow the exact same test shape as Task 4's `test_polish_with_exonerate_*`
tests, adapted for miniprot's real output format. Before writing the fake
runner's output string, check miniprot's actual GFF output format (miniprot
supports `--gff` for GFF3-shaped output with `mRNA`/`CDS` features) --
confirm the real feature/attribute names via `miniprot --help` or by
running it once manually if available in this environment, rather than
guessing the format. Write the test's fake output to match the REAL format
you find, the same way Task 4a/Task 4 in the original SDD run corrected
the plan's illustrative header formats against the real tools.

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/detect/test_search.py -v -k polish_with_miniprot`
Expected: FAIL.

- [ ] **Step 3: Implement**

Mirror `polish_with_exonerate`'s structure (slice window via
`_extract_window`, invoke miniprot with `--gff` output, parse `mRNA`/`CDS`
features into `ExonSpan`s, re-base coordinates by the window's offset,
attribute via `_attribute`). Use your judgment on the exact miniprot
CLI invocation and GFF field parsing based on the real tool's actual
interface (verified in Step 1), following the same offset-rebasing and
attribution pattern as `polish_with_exonerate`.

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/detect/test_search.py -v -k polish_with_miniprot`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/detect/search.py tests/detect/test_search.py
git commit -m "$(cat <<'EOF'
feat: add miniprot polishing wrapper with exon parsing

Mirrors polish_with_exonerate's window-slicing and offset-rebasing
approach; parses miniprot's real GFF mRNA/CDS output into PolishModel
exon structure for cross-tool agreement comparison.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Per-gene status classification and canonical evidence selection

**Files:**
- Modify: `src/MATPredict/detect/polish.py`
- Test: `tests/detect/test_polish.py`

**Interfaces:**
- Consumes: `PolishModel`/`boundaries_agree` (Task 3), `GeneEvidence` (`pipeline.py`, existing).
- Produces:
  ```python
  STATUS_AGREE = "polished_agree"
  STATUS_DISAGREE = "polished_disagree"
  STATUS_SINGLE = "polished_single"
  STATUS_UNPOLISHED = "unpolished"

  @dataclass(frozen=True)
  class PolishOutcome:
      status: str  # one of the four STATUS_* constants
      canonical: PolishModel | None  # None only when status == STATUS_UNPOLISHED
      exonerate_model: PolishModel | None
      miniprot_model: PolishModel | None

  def classify(
      exonerate_model: PolishModel | None,
      miniprot_model: PolishModel | None,
      tolerance_bp: int = 10,
  ) -> PolishOutcome:
      """Classify a gene's two-tool polish outcome. exonerate_model is
      preferred as canonical when both models exist (more directly
      integrated, matches sub-project 1's curation conventions) --
      miniprot's model is retained for comparison either way, never
      discarded. Both None with no fallback is an error at the caller
      (an unpolished status needs the raw tblastn hit, which this
      function does not have -- the caller in pipeline.py builds
      STATUS_UNPOLISHED's canonical from the tblastn SearchHit directly,
      not through this function)."""
  ```

  Note: `classify` only ever returns `STATUS_UNPOLISHED` when BOTH inputs
  are `None`, with `canonical=None` — `pipeline.py` (Task 8) is
  responsible for substituting the raw `tblastn` hit's coordinates in
  that case, since this module has no access to the original `SearchHit`.

- [ ] **Step 1: Write the failing test**

```python
# tests/detect/test_polish.py (add)
from MATPredict.detect.polish import (
    STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE, STATUS_UNPOLISHED, classify,
)


def _model(start, end, method):
    return PolishModel(
        gene_name="mfa1", family_key=KEY, role="core_MAT", contig="c1",
        start=start, end=end, strand="+", exons=[ExonSpan(start, end)],
        identity=95.0, reference_record_id="rec1", method=method,
    )


def test_classify_agree_prefers_exonerate_as_canonical():
    ex = _model(100, 400, "exonerate_refine")
    mp = _model(102, 398, "miniprot")
    outcome = classify(ex, mp, tolerance_bp=10)
    assert outcome.status == STATUS_AGREE
    assert outcome.canonical is ex
    assert outcome.miniprot_model is mp


def test_classify_disagree_keeps_both_models():
    ex = _model(100, 400, "exonerate_refine")
    mp = _model(500, 900, "miniprot")
    outcome = classify(ex, mp, tolerance_bp=10)
    assert outcome.status == STATUS_DISAGREE
    assert outcome.canonical is ex
    assert outcome.exonerate_model is ex and outcome.miniprot_model is mp


def test_classify_single_tool_uses_whichever_succeeded():
    mp = _model(100, 400, "miniprot")
    outcome = classify(None, mp, tolerance_bp=10)
    assert outcome.status == STATUS_SINGLE
    assert outcome.canonical is mp


def test_classify_unpolished_when_neither_tool_succeeds():
    outcome = classify(None, None, tolerance_bp=10)
    assert outcome.status == STATUS_UNPOLISHED
    assert outcome.canonical is None
```

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/detect/test_polish.py -v -k classify`
Expected: FAIL.

- [ ] **Step 3: Implement**

```python
# add to src/MATPredict/detect/polish.py
STATUS_AGREE = "polished_agree"
STATUS_DISAGREE = "polished_disagree"
STATUS_SINGLE = "polished_single"
STATUS_UNPOLISHED = "unpolished"


@dataclass(frozen=True)
class PolishOutcome:
    status: str
    canonical: PolishModel | None
    exonerate_model: PolishModel | None
    miniprot_model: PolishModel | None


def classify(
    exonerate_model: PolishModel | None,
    miniprot_model: PolishModel | None,
    tolerance_bp: int = 10,
) -> PolishOutcome:
    if exonerate_model is not None and miniprot_model is not None:
        status = STATUS_AGREE if boundaries_agree(exonerate_model, miniprot_model, tolerance_bp) else STATUS_DISAGREE
        return PolishOutcome(status, exonerate_model, exonerate_model, miniprot_model)
    if exonerate_model is not None or miniprot_model is not None:
        canonical = exonerate_model or miniprot_model
        return PolishOutcome(STATUS_SINGLE, canonical, exonerate_model, miniprot_model)
    return PolishOutcome(STATUS_UNPOLISHED, None, None, None)
```

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/detect/test_polish.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/detect/polish.py tests/detect/test_polish.py
git commit -m "$(cat <<'EOF'
feat: classify per-gene polish outcomes into 4 explicit statuses

polished_agree/polished_disagree/polished_single/unpolished, with
exonerate preferred as canonical when both tools succeed and neither
model ever silently discarded. unpolished's canonical is left None --
pipeline.py substitutes the raw tblastn hit for that case, since this
module has no access to it.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Retire `second_pass_used` in favor of `unpolished`-driven tiering

**Files:**
- Modify: `src/MATPredict/detect/tiering.py`
- Test: `tests/detect/test_tiering.py`

**Interfaces:**
- Produces: `assign_tier`'s `second_pass_used: bool` parameter is renamed
  `any_gene_unpolished: bool`, with identical Medium-tier-cap semantics —
  this is a rename reflecting the new source of truth (Task 8 computes it
  from `PolishOutcome.status == STATUS_UNPOLISHED` across a family's
  genes, not from a relaxed-exonerate flag), not a behavior change to
  `assign_tier`'s own logic.

- [ ] **Step 1: Update the existing tests**

In `tests/detect/test_tiering.py`, rename every `second_pass_used=` keyword
argument to `any_gene_unpolished=` (the boolean values passed do not
change — a `True` in the old tests meant "the relaxed second pass was
needed", which is the same situation Task 8 will now flag as `True` when
any of a family's genes has `PolishOutcome.status == STATUS_UNPOLISHED`).

- [ ] **Step 2: Run to verify it fails**

Run: `pixi run pytest tests/detect/test_tiering.py -v`
Expected: FAIL (`TypeError: assign_tier() got an unexpected keyword argument`).

- [ ] **Step 3: Implement**

In `tiering.py`, rename the parameter (function body logic is otherwise
unchanged — this is purely a rename to reflect the new semantics):

```python
def assign_tier(
    score: FamilyScore,
    family: Family,
    cluster: GeneCluster,
    any_gene_unpolished: bool,
    fragmented: bool,
) -> str:
    ...
    elif any_gene_unpolished:
        tier = "medium"
    ...
```

Update the module docstring's reference to "second pass" language to
describe the new `unpolished` status instead.

- [ ] **Step 4: Run to verify it passes**

Run: `pixi run pytest tests/detect/test_tiering.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/detect/tiering.py tests/detect/test_tiering.py
git commit -m "$(cat <<'EOF'
refactor: rename assign_tier's second_pass_used to any_gene_unpolished

Same Medium-tier-cap semantics, new name reflecting the localize-
then-polish revision's source of truth (a gene left "unpolished" by
both miniprot and exonerate --refine), not a relaxed-exonerate flag
that no longer exists after this revision.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Rewire `pipeline.py`'s Stage 0/1/2/3 flow

This is the integration task — the highest-risk task in this plan, given
this project's history of subtle cross-family/cross-cluster attribution
bugs in exactly this file (Tasks 6, 7, and 9 of the original SDD run each
shipped one, all caught by careful review).

**Files:**
- Modify: `src/MATPredict/detect/pipeline.py`
- Test: `tests/detect/test_pipeline.py`

**Interfaces:**
- Consumes: `search_localize` (Task 2), `polish_with_exonerate`/
  `polish_with_miniprot` (Tasks 4-5), `polish.classify`/`PolishOutcome`
  (Task 6), `assign_tier`'s renamed parameter (Task 7).
- Produces: `run_pipeline`'s public signature and `DetectionOutcome`
  shape are unchanged at the type level; `GeneEvidence` gains new fields
  (see Task 9) that `run_pipeline` now populates from `PolishOutcome`
  instead of directly from a `SearchHit`.
- Removes: `_families_with_a_foothold`'s relaxed-exonerate windowed
  second-pass block (the `for family in _families_with_a_foothold(...)`
  loop calling `search_genomic(..., relaxed=True, window=...)`), the
  whole-genome `genome_wide_second_pass` block (the `if missing_families:`
  block calling `search_genomic(..., relaxed=True)` with no window), and
  `second_pass_used_for`/`genome_wide_second_pass`/`_second_pass_used` —
  all replaced by the logic below.

### New control flow

1. **Genome-only path** (`proteome_fasta is None`): call `search_localize`
   (Task 2) once, batched across every routed family — this replaces the
   `else: hits.extend(search_genomic(...))` branch entirely. Cluster the
   resulting `tblastn` hits with the existing `cluster_hits` (unchanged).
2. **Fast-path** (`proteome_fasta is not None`): call `search_fast_path`
   as today (unchanged), cluster, then for each cluster/family with a
   foothold missing a core gene (reusing `_families_with_a_foothold` and
   `_missing_core_genes` unchanged), **skip Stage 1 entirely** and go
   directly to Stage 2 (polish) with a window padded around the existing
   cluster's span, querying only the missing gene's curated proteins.
3. **Polish** (both paths, for every gene in every cluster that needs it —
   in the genome-only path, this is every gene the family expects but
   whose canonical location still needs precision refinement beyond a raw
   `tblastn` HSP; in the fast-path rescue, only the missing gene): for
   each such gene, call both `polish_with_exonerate` and
   `polish_with_miniprot` against the same padded window, then
   `polish.classify` the pair.
4. **Fold `PolishOutcome`s back into evidence**: `GeneEvidence` for a
   polished gene uses `PolishOutcome.canonical`'s coordinates/identity/
   method (Task 9 extends `GeneEvidence` to also carry the `status` and
   both tools' models). A gene with `status == STATUS_UNPOLISHED` uses
   the original raw `tblastn` (or, for the fast-path rescue, whichever
   hit localized it) `SearchHit`'s coordinates directly, flagged
   `unpolished` — this is the one case where `GeneEvidence` is built from
   a `SearchHit`, not a `PolishModel`.
5. **Tiering input**: for each family's `DetectionResult`, compute
   `any_gene_unpolished = any(status == STATUS_UNPOLISHED for status in
   this family's genes' PolishOutcome.status)` and pass it to
   `assign_tier` in place of the old `second_pass_used` argument.
   `polished_agree` vs `polished_disagree` must have **zero** effect on
   this computation or on `assign_tier`'s result — write a test proving
   this explicitly (see Step 1 below).

**Padding margin**: derive the window padding per gene from that gene's
own curated reference protein's length (per the spec: "a multiple of the
protein's nucleotide-equivalent length, plus a max-intron allowance") —
add a small helper function for this rather than a single hardcoded
constant; make the multiplier and max-intron-allowance both named,
overridable parameters with sensible defaults (document your chosen
defaults and rationale in the function's docstring).

- [ ] **Step 1: Write the failing tests**

In `tests/detect/test_pipeline.py`, add (at minimum) these cases, each
using injected stub functions for `search_localize`,
`polish_with_exonerate`, `polish_with_miniprot` (following this file's
existing pattern of injecting fake search functions):

```python
def test_genome_only_path_uses_search_localize_not_search_genomic():
    """The genome-only path must call the injected search_localize, and
    must never call search_genomic at all (search_genomic's old
    whole-genome/windowed-relaxed roles are fully retired)."""

def test_fast_path_missing_gene_rescue_skips_localization_and_polishes_directly():
    """A fast-path cluster missing one core gene must go straight to
    polish_with_exonerate/polish_with_miniprot against a window padded
    around the EXISTING cluster's span, without ever calling
    search_localize."""

def test_polished_agree_and_disagree_produce_identical_tier():
    """Two otherwise-identical scenarios, one where the two polish tools
    agree and one where they disagree on the SAME gene's boundaries,
    must produce the SAME confidence tier -- proving agreement is
    reported but not consulted by assign_tier."""

def test_unpolished_gene_caps_tier_at_medium():
    """A family whose genes are all found but one gene's polish outcome
    is STATUS_UNPOLISHED (neither tool produced a model) reaches at most
    Medium, never High -- same effect the old second_pass_used had."""

def test_gene_evidence_for_unpolished_gene_uses_raw_localization_hit():
    """A STATUS_UNPOLISHED gene's GeneEvidence carries the raw tblastn
    (or fast-path localization) hit's coordinates, not a fabricated or
    missing value."""
```

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/detect/test_pipeline.py -v`
Expected: multiple FAILures (old functions/parameters still referenced,
new behavior not implemented yet).

- [ ] **Step 3: Implement**

Rewrite `run_pipeline`'s body per the control flow above. Reuse
`_missing_core_genes`, `_families_with_a_foothold`, `_fragmented_family_segments`,
`_contig_lengths`, `_segments_for` unchanged. Replace `_gene_evidence`'s
current "best hit by identity across `SearchHit`s" logic with a version
that consumes `PolishOutcome`s (falling back to raw `SearchHit`s for
`STATUS_UNPOLISHED` genes) per point 4 above — this is a real rewrite of
`_gene_evidence`, not an extension, since its input shape changes from
"list of SearchHit" to "per-gene PolishOutcome, or a raw SearchHit for
unpolished genes".

Remove the `genome_wide_second_pass`/`second_pass_used_for`/
`_second_pass_used` machinery entirely (per "Removes" above), and delete
the module docstring's references to "relaxed genomic second pass" in
favor of describing the new localize-then-polish flow.

- [ ] **Step 4: Run to verify they pass**

Run: `pixi run pytest tests/detect/test_pipeline.py -v`
Expected: all PASS, zero regressions in tests not touched by this task
(the fragmentation, ambiguity, and idiomorph-assignment tests from the
original SDD run should be unaffected by this rewrite — verify this
explicitly, since this task touches the same file those tests exercise).

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/detect/pipeline.py tests/detect/test_pipeline.py
git commit -m "$(cat <<'EOF'
feat: rewire run_pipeline for localize-then-polish (Stage 0-3)

Genome-only path now uses batched tblastn localization instead of an
unrestricted whole-genome exonerate fallback. Fast-path per-family
rescue skips localization entirely (a rough location already exists
from the existing cluster) and polishes directly. Every gene's
two-tool polish outcome (agree/disagree/single/unpolished) feeds
GeneEvidence; only "unpolished" affects confidence tiering, proven by
an explicit test that agree vs disagree produce identical tiers.
Retires the relaxed-exonerate second-pass machinery entirely.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Extend `GeneEvidence`/`DetectionResult`/`report.py` for per-gene status

**Files:**
- Modify: `src/MATPredict/detect/pipeline.py` (dataclass definitions only — `run_pipeline`'s logic was already updated in Task 8 to populate these; this task adds the fields those changes needed)
- Modify: `src/MATPredict/detect/report.py`
- Test: `tests/detect/test_report.py`, `tests/detect/test_pipeline.py`

**Interfaces:**
- Produces: `GeneEvidence` gains:
  ```python
  @dataclass(frozen=True)
  class GeneEvidence:
      gene_name: str
      role: str
      contig: str
      start: int
      end: int
      strand: str
      identity: float
      coverage: float | None
      reference_record_id: str
      method: str
      status: str = "polished_agree"  # one of polish.STATUS_* -- default
                                        # only for backward-compatible
                                        # construction in existing tests
                                        # that predate this field
      alternate_model: dict | None = None  # the OTHER tool's model
                                             # (contig/start/end/exons/
                                             # identity/method) when
                                             # status == polished_disagree,
                                             # else None
  ```
  (This is a genuine shape change, per the spec's explicit note — treat
  it as such, not a same-shape extension. Existing tests constructing
  `GeneEvidence` without `status`/`alternate_model` continue to work via
  the defaults, but Task 8's real pipeline code must always pass a real
  `status` explicitly, never rely on the default.)

- [ ] **Step 1: Update/add tests**

In `tests/detect/test_report.py`, add a case where a `GeneEvidence` has
`status="polished_disagree"` and a populated `alternate_model`, asserting
both the canonical coordinates AND the alternate model's data round-trip
through `write_detection_report`'s YAML output. Add a GFF3 case asserting
the gene feature's attributes include the new `status` field (e.g.
`status=polished_disagree`).

- [ ] **Step 2: Run to verify they fail**

Run: `pixi run pytest tests/detect/test_report.py -v`
Expected: FAIL (`TypeError`/missing keys in the assertions).

- [ ] **Step 3: Implement**

Add the two fields to `GeneEvidence` in `pipeline.py` per the Interfaces
block above. In `report.py`:
- `write_detection_gff3`'s gene-feature attribute string gains
  `;status={evidence.status}` (and, when `alternate_model` is present, a
  compact encoding of it — e.g.
  `;alt_method={m};alt_start={s};alt_end={e}` — your judgment on the
  exact attribute names, consistent with the file's existing style).
- `_result_doc`'s `gene_evidence` dict gains `"status"` and
  `"alternate_model"` keys (the latter as a nested dict or `None`).

- [ ] **Step 4: Run to verify they pass**

Run: `pixi run pytest tests/detect/test_report.py -v tests/detect/test_pipeline.py -v`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/detect/pipeline.py src/MATPredict/detect/report.py tests/detect/test_report.py tests/detect/test_pipeline.py
git commit -m "$(cat <<'EOF'
feat: surface per-gene polish status and disagreeing model in output

GeneEvidence carries the polish status (agree/disagree/single/
unpolished) and, when the two tools disagreed, the non-canonical
tool's full model -- both the YAML report and GFF3 now expose this so
a human reviewer can inspect a disagreement rather than only ever
seeing the canonical (exonerate-preferred) coordinates.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Retire dead code and audit every `method`/reason-text reference

**Files:**
- Modify: `src/MATPredict/detect/search.py` (remove `search_genomic`'s
  `relaxed`/`window` parameters if Task 8 left them unused — confirm
  first; `search_genomic` itself may still be needed by `polish_with_exonerate`
  indirectly or may be fully superseded — check before deleting anything)
- Modify: `src/MATPredict/detect/pipeline.py` (the `NotDetectedFamily`
  reason text referencing "the relaxed whole-genome second pass")
- Modify: `tests/detect/test_search.py`, `tests/detect/test_pipeline.py`
  (any test still asserting the retired `exonerate_genome`/
  `exonerate_genome_relaxed` method strings or the old
  `second_pass_used_for`/`genome_wide_second_pass` mechanics)
- Modify: `pixi.toml` docstring/comment if `exonerate`'s stated role
  needs updating now that its genomic-fallback role moved to
  `polish_with_exonerate`

**Interfaces:** none new — this is a cleanup/consistency pass.

- [ ] **Step 1: Grep for every retired reference**

Run, from the repo root:
```bash
grep -rn "exonerate_genome\b\|exonerate_genome_relaxed\|second_pass_used_for\|genome_wide_second_pass\|relaxed whole-genome second pass" src/ tests/
```
Each hit is either: (a) genuinely retired code/text that must be
removed/updated, or (b) a legitimate remaining reference (e.g. if
`search_genomic`'s original whole-genome mode is still used somewhere
this plan didn't anticipate — investigate before deleting).

- [ ] **Step 2: Decide `search_genomic`'s fate**

Confirm whether anything still calls `search_genomic` with `window=None`
(the old unrestricted whole-genome mode) after Task 8's rewrite. If
nothing does, either delete `search_genomic` and its now-dead tests, or
keep it only if `polish_with_exonerate`/`polish_with_miniprot` internally
reuse parts of it (they were specified as new, standalone functions in
Tasks 4-5, so this is likely dead code to remove — verify, don't assume).

- [ ] **Step 3: Update the `NotDetectedFamily` reason text**

In `pipeline.py`, change the reason string that currently reads
`"...including after the relaxed whole-genome second pass"` to describe
the new flow, e.g. `"...including after genome-wide tblastn localization
and polishing"`.

- [ ] **Step 4: Run the full suite**

Run: `pixi run pytest -v`
Expected: 100% pass, and the grep from Step 1 returns no remaining
unintentional hits.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "$(cat <<'EOF'
chore: retire dead exonerate-relaxed code and update stale reason text

Removes references to the deleted relaxed-exonerate second-pass
mechanism (method strings, reason text, and search_genomic's
unrestricted-genome mode if nothing else calls it) now that
localize-then-polish has fully replaced it.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: End-to-end integration test for the localize-then-polish flow

**Files:**
- Create or extend: `tests/detect/test_integration_localize_polish.py`

**Interfaces:**
- Consumes: the real, already-accepted *Mycosarcoma maydis* a-locus
  record (`db/Basidiomycota/Ustilaginales/5270_521_aLocus_a1/`, confirmed
  in the original SDD run's Task 14 — verify the path still exists before
  writing this test, per that task's own lesson about the plan's path
  guesses needing verification) as real curated reference data.

- [ ] **Step 1: Write the test**

Build on the pattern from the original `test_integration_real_record.py`:
stub `search_localize`, `polish_with_exonerate`, and `polish_with_miniprot`
(no live binaries) to return realistic data for the real `aLocus`
family's genes (`mfa1`, `pra1`), run `run_pipeline` with `taxid=5270`
against the real `db/`, and assert:
- the result reaches "high" confidence with both genes in `genes_found`,
- `gene_evidence` entries carry `status="polished_agree"` (from stubbed
  matching models) rather than a default/placeholder value,
- switching one gene's stubbed models to disagree (different
  coordinates) still produces the same overall tier (proving the
  integration-level version of Task 8's agree/disagree-tier-independence
  test).

- [ ] **Step 2: Run to verify it passes**

Run: `pixi run pytest tests/detect/test_integration_localize_polish.py -v`
Expected: PASS (or a clear, graceful `pytest.mark.skipif` if the real
record's path has moved — confirm first, don't guess).

- [ ] **Step 3: Run the full suite one more time**

Run: `pixi run pytest -v`
Expected: 100% pass, no regressions anywhere in the suite.

- [ ] **Step 4: Commit**

```bash
git add tests/detect/test_integration_localize_polish.py
git commit -m "$(cat <<'EOF'
test: add end-to-end integration check for localize-then-polish

Exercises the real routing/family/order.yml data for the M. maydis
a-locus record through the full Stage 0-3 flow with stubbed (not
live) search_localize/polish_with_exonerate/polish_with_miniprot,
including an explicit check that a disagreeing-models scenario
produces the same tier as an agreeing one.

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Benchmark/timing acceptance evidence (manual, not unit-tested)

The spec's acceptance criteria require (a) running the existing
leave-one-out benchmark against `assembly`-type curated records to
compare gene recall/boundary accuracy against a pre-revision baseline,
and (b) a measured wall-clock comparison on a named real assembly. Both
require live binaries and real data, so this task is explicitly NOT
automated pytest coverage — it is a documented, one-time verification
step for whoever merges this plan's work.

**Files:**
- Create: `docs/superpowers/plans/2026-09-17-mat-detection-search-localization-benchmark-notes.md`
  (a short, plain-language record of what was run and what it showed —
  not a spec, not code, just an evidence log)

- [ ] **Step 1: Identify a real `assembly`-type curated record**

Query `db/matpredict.duckdb` (rebuild via `pixi run matpredict curate-db
build-duckdb` if stale) or grep `db/**/metadata.yaml` for a record whose
`locus.core.segments[].sequence_source.type == "assembly"` — this is the
kind of record the leave-one-out benchmark needs to test against a real,
full genome rather than a locus-only fragment.

- [ ] **Step 2: Run `matpredict detect` against that assembly, before and after this plan's changes**

If feasible in this environment (real `tblastn`/`miniprot`/`exonerate`
binaries installed, per Task 1), time an actual `matpredict detect`
invocation against the chosen assembly on the pre-revision code
(`git stash`/a separate checkout of the commit before Task 1) and again
on the post-revision code. Record both wall-clock times and both runs'
confidence/gene-recall outcomes for the same family.

- [ ] **Step 3: Record the results**

Write the findings (timings, outcome comparison, and — if either
measurement wasn't feasible in this environment, say so plainly rather
than fabricating a number) to the notes file created above.

- [ ] **Step 4: Commit**

```bash
git add docs/superpowers/plans/2026-09-17-mat-detection-search-localization-benchmark-notes.md
git commit -m "$(cat <<'EOF'
docs: record localize-then-polish benchmark/timing evidence

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>
EOF
)"
```

---

## Final Review Checklist (run once all 12 tasks are complete)

- [ ] `pixi run pytest -v` passes with zero failures across the whole suite.
- [ ] No test in `tests/detect/` shells out to a live `tblastn`/`miniprot`/
  `exonerate`/`makeblastdb` binary — grep for bare `subprocess.run(` calls
  without an injected `runner=` parameter.
- [ ] `matpredict detect --genome <fasta> --taxid <n> --out-dir <dir>` runs
  end-to-end against a real genome (even a small held-out contig) using
  the real installed binaries at least once, confirming the wrapper
  invocations are syntactically valid against the real tools (not just
  passing against fakes).
- [ ] Every Fable-review decision this plan encodes is verifiably present
  in the shipped code: `unpolished` (not a relaxed-exonerate flag) drives
  the Medium-tier cap; `polished_agree`/`polished_disagree` produce
  identical tiers; the fast-path rescue skips Stage 1; the genome-only
  path batches one `tblastn` call across all families; every retired
  `method` string and reason-text reference is gone.
- [ ] The Task 12 benchmark notes file exists and honestly reports
  whatever could and could not be measured in this environment.
