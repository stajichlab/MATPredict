# Curation and Calibration Follow-ups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close out the follow-up work from the detection-rollout-fixes plan, in the user's specified order: (1-3) resolve the 3 curated records with no protein sequence reaching the reference FASTA, (4) wire the Stage1→Stage2 evidence-floor gate into real entry points and fix its `min_hits` semantics, (5) investigate why genome 5501's own real MAT locus was only weakly detected, and (6) re-run the full 13-genome pilot rollout with all of this session's fixes in place.

**Architecture:** Six tasks. Tasks 1-3 are a code task (extend the DB build tooling to derive a protein sequence from already-curated genomic coordinates when no `protein_accession` exists) followed by two data/curation tasks that use it for real records. Task 4 is a code task wiring an already-built, already-tested mechanism into two real entry points and fixing one real semantic bug. Tasks 5-6 are analysis/execution steps using the now-complete, now-corrected pipeline.

**Tech Stack:** Python, this project's existing `NcbiClient`/`_independent_translation`/`gff_export`/`family_registry`/`pipeline` modules, real `tblastn` (already a pixi dependency), `pixi run pytest`.

**Spec:** No separate spec doc — this plan implements the follow-up items recorded in the previous plan's SDD ledger (`docs/superpowers/plans/2026-09-18-detection-rollout-fixes.md`'s own final-review findings I-1/I-2, and the 3-null-accession-record curation gap first surfaced in that plan's Task 2).

## Global Constraints

- Never build a second NCBI client — reuse `MATPredict.db.ncbi_client.NcbiClient` and its existing `_independent_translation`/`fetch_nucleotide_sequence` methods everywhere in this plan.
- Never fabricate a coordinate, accession, or protein sequence — every curated-data change in this plan must be backed by a real, reproducible fetch or search, with the evidence documented in the record's own `definition_note`/`evidence` blocks, matching this project's existing curation discipline (see the 3 records' own current notes for the tone and rigor expected).
- `db/candidates/` records are never authoritative and are out of scope for this plan.
- No change to `run_pipeline`'s default behavior for callers that don't opt into `EvidenceFloor`/`evidence_diagnostics_path` — Task 4 only wires existing, already-default-no-op parameters into new call sites.
- Live-verify any new accession/coordinate claim against real NCBI data before writing it into `order.yml` or a curated record, per this project's standing rule against unverified claims.

---

### Task 1: Support deriving a protein sequence from curated genomic coordinates when no `protein_accession` exists

**Files:**
- Modify: `src/MATPredict/db/cli.py` (`build_gff_for_record`)
- Test: `tests/db/test_backfill_gff.py`

**Interfaces:**
- Produces: `build_gff_for_record`'s existing per-gene loop gains a fallback path — when a present gene has no `protein_accession` but has real fetchable genomic coordinates (a `segment_index` pointing at a segment whose `sequence_source.type == "insdc_nucleotide"`), it derives the sequence via `MATPredict.db.validate._independent_translation(record, gene, ncbi)` instead of skipping the gene.

**Background:** Two accepted records (`db/Ascomycota/Teloschistales/2903220_liq80xsp_MAT_combined`, `db/Ascomycota/Teloschistales/2903222_liq146xsp_MAT_combined`) have real, curator-verified genomic coordinates for all 4 of their genes (each confirmed to begin with ATG and end with a stop codon under this project's 1-based fully-closed convention — see each gene's segment/coordinates and the record's own `definition_note`), but these are unannotated MAG assemblies with no NCBI protein records at all — there is no accession to fetch, ever, for these genes. `_independent_translation` (already implemented in `src/MATPredict/db/validate.py:54-103`, currently used only to cross-check a claimed `protein_accession`'s sequence against the genome) already does exactly the fetch-and-translate work needed here — it takes a `record` dict and a `gene` dict, resolves the gene's own segment's `sequence_source.accession`, fetches the sequence (single span, or spliced across `exons` if present), applies `codon_start` (default 1) and `transl_table` (default 1, standard code — correct for these nuclear MAT genes), and returns a real translated protein string, or `None` if the coordinate data isn't fetchable. This task wires that existing function into the build path as a fallback, rather than writing new translation logic.

- [ ] **Step 1: Write the failing test**

```python
# append to tests/db/test_backfill_gff.py
from unittest.mock import MagicMock

from MATPredict.db.cli import build_gff_for_record


def _write_record_with_coordinates_no_accession(db_root, phylum, order_or_family, record_id):
    """A gene with real segment coordinates but protein_accession: null -- the
    Xanthoria MAG case: a real genomic span exists, no NCBI protein record does."""
    record_dir = db_root / phylum / order_or_family / record_id
    record_dir.mkdir(parents=True)
    metadata = {
        "record_id": record_id,
        "locus": {"core": {"segments": [
            {"sequence_source": {"type": "insdc_nucleotide", "accession": "ACC1.1"}},
        ]}},
        "genes": [
            {
                "gene_index": 0, "name": "G1", "role": "core_MAT", "present": True,
                "protein_accession": None, "segment_index": 0,
                "start": 100, "end": 108, "strand": "+",
            },
        ],
    }
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(metadata))
    return record_dir


def test_build_gff_for_record_derives_sequence_from_coordinates_when_no_accession(tmp_path, monkeypatch):
    db_root = tmp_path / "db"
    record_dir = _write_record_with_coordinates_no_accession(
        db_root, "Ascomycota", "Teloschistales", "111_a_MAT_combined"
    )

    monkeypatch.setattr(
        "MATPredict.db.cli._independent_translation",
        lambda record, gene, ncbi: "MSEQ",
    )

    build_gff_for_record(db_root, "Ascomycota", "Teloschistales", "111_a_MAT_combined", ncbi=MagicMock(), uniprot=MagicMock())

    faa_text = (record_dir / "proteins.faa").read_text()
    assert ">111_a_MAT_combined|gene_index=0|name=G1|role=core_MAT" in faa_text
    assert "MSEQ" in faa_text


def test_build_gff_for_record_skips_gene_when_translation_returns_none(tmp_path, monkeypatch):
    db_root = tmp_path / "db"
    _write_record_with_coordinates_no_accession(db_root, "Ascomycota", "Teloschistales", "222_b_MAT_combined")

    monkeypatch.setattr("MATPredict.db.cli._independent_translation", lambda record, gene, ncbi: None)

    build_gff_for_record(db_root, "Ascomycota", "Teloschistales", "222_b_MAT_combined", ncbi=MagicMock(), uniprot=MagicMock())

    faa_text = (db_root / "Ascomycota" / "Teloschistales" / "222_b_MAT_combined" / "proteins.faa").read_text()
    assert faa_text.strip() == ""
```

- [ ] **Step 2: Run to verify it fails**

```bash
pixi run pytest tests/db/test_backfill_gff.py -v
```

Expected: FAIL (`ImportError: cannot import name '_independent_translation' from 'MATPredict.db.cli'` or an `AttributeError` on the monkeypatch target, since it isn't imported into `cli.py` yet).

- [ ] **Step 3: Implement the fallback in `build_gff_for_record`**

In `src/MATPredict/db/cli.py`, add `from MATPredict.db.validate import _client_for, _independent_translation` (extending the existing `from MATPredict.db.validate import _client_for, validate_record` import line — check its exact current form first) and change `build_gff_for_record`'s gene loop:

```python
def build_gff_for_record(
    db_root: Path, phylum: str, order_or_family: str, record_id: str,
    ncbi: NcbiClient, uniprot: UniprotClient,
) -> None:
    """Fetch every present gene's protein sequence and write
    locus.gff3/locus.gbk/proteins.faa for one accepted record. A gene with a
    real `protein_accession` is fetched directly; a gene with none but real
    curated genomic coordinates (an unannotated MAG assembly with no NCBI
    protein record to cite) has its sequence independently derived via
    `_independent_translation`, which fetches and translates from the
    record's own recorded coordinates -- the same function `validate.py`
    already uses to cross-check a claimed accession's sequence, reused here
    as the sequence SOURCE rather than a cross-check target. A gene with
    neither a protein_accession nor derivable coordinates is skipped, same
    as before. The single-record CLI command and the batch backfill command
    both call this so there is exactly one place this logic lives."""
    record_dir = db_root / phylum / order_or_family / record_id
    record = yaml.safe_load((record_dir / "metadata.yaml").read_text())

    sequences: dict[int, str] = {}
    for gene in record.get("genes", []):
        if not gene.get("present", True):
            continue
        if gene.get("protein_accession"):
            client, bare_accession = _client_for(gene["protein_accession"], ncbi, uniprot)
            sequences[gene["gene_index"]] = client.fetch_protein_sequence(bare_accession)
            continue
        derived = _independent_translation(record, gene, ncbi)
        if derived:
            sequences[gene["gene_index"]] = derived

    gff_export.write_gff3(record, out_path=record_dir / "locus.gff3")
    gff_export.write_genbank(record, sequences, out_path=record_dir / "locus.gbk")
    gff_export.write_proteins_fasta(record, sequences, out_path=record_dir / "proteins.faa")
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run pytest tests/db/test_backfill_gff.py -v
```

Expected: PASS, all tests including the 2 new ones.

- [ ] **Step 5: Run the full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/db/cli.py tests/db/test_backfill_gff.py
git commit -m "feat: derive a gene's protein sequence from curated genomic coordinates when no protein_accession exists"
```

---

### Task 2: Apply Task 1's capability to the 2 Xanthoria records (data step)

**Files:** none new — this task runs Task 1's tooling for real against 2 specific records and verifies the output.

- [ ] **Step 1: Confirm each of the 8 genes' coordinate data is real and already-verified** (it is — read both records' `definition_note`s again yourself, don't re-derive from scratch) before running the backfill: `2903220_liq80xsp_MAT_combined`'s APN2/MAT1-2-1/MAT1-1-1/SLA2 and `2903222_liq146xsp_MAT_combined`'s same 4 genes, all confirmed ATG-to-stop under 1-based fully-closed coordinates per the curator's own notes.

- [ ] **Step 2: Run the real backfill**

```bash
pixi run matpredict curate-db backfill-gff
```

Since both records already have `proteins.faa` present (empty, from the previous plan's backfill), they will NOT be picked up by `find_records_missing_proteins_faa` (which — after this plan's Task 1 — still uses the "empty file counts as missing" check from the prior plan). Confirm this is the case, then run each one directly instead:

```bash
pixi run python -c "
from pathlib import Path
from MATPredict.db.cli import build_gff_for_record, _make_clients
from MATPredict.config import MatpredictConfig
config = MatpredictConfig.from_env(repo_root=Path.cwd())
ncbi, uniprot = _make_clients(config)
for record_id in ['2903220_liq80xsp_MAT_combined', '2903222_liq146xsp_MAT_combined']:
    build_gff_for_record(config.db_root, 'Ascomycota', 'Teloschistales', record_id, ncbi, uniprot)
    print('rebuilt', record_id)
"
```

- [ ] **Step 3: Verify the real output**

```bash
cat db/Ascomycota/Teloschistales/2903220_liq80xsp_MAT_combined/proteins.faa
cat db/Ascomycota/Teloschistales/2903222_liq146xsp_MAT_combined/proteins.faa
```

Confirm: 4 real FASTA entries in each file (APN2, MAT1-2-1, MAT1-1-1, SLA2), each a plausible amino-acid sequence (standard 20-letter alphabet, no error text), each ending without a trailing stop codon (`_translate_cds` already stops at and drops the first stop codon per its own docstring). Cross-check at least one gene's length against its recorded genomic span: e.g. `2903222`'s `APN2` is `32002..33995` (1994 nt) on the `+` strand with no `exons` recorded (single span) — expect a translated protein of `(1994 // 3) - 1` ≈ 663 residues (minus 1 for the dropped stop codon), give or take a few residues depending on where the actual in-frame stop lands. **Do not accept a suspiciously short (a few residues) or empty sequence without investigating** — a genuinely wrong reading frame or an off-by-one in the stored coordinates would show up exactly this way.

- [ ] **Step 4: Update each record's own `validation.sequence_match` block to reflect the new, real independently-derived sequence**

Both records' current `sequence_match.notes` say "No protein-level check was possible... this assembly carries no NCBI protein annotation, so protein_accession is null for every gene" — this is now only half true (there's still no protein_accession, but there IS now a real derived sequence). Update each record's `validation.sequence_match.notes` to state plainly that the protein sequence was independently derived via translation of the record's own already-verified genomic coordinates (not fetched from an external protein accession, since none exists), and name the commit/date. Do not change `sequence_match.status` (leave as `pass`, since the underlying coordinate verification this reflects is unchanged) — only the notes' explanation of *why* a protein sequence now exists needs updating.

- [ ] **Step 5: Rebuild the reference FASTA and confirm the new sequences reach it**

```bash
pixi run python -c "
from pathlib import Path
from MATPredict.detect.reference_fasta import build_reference_fasta
build_reference_fasta(Path('db'), Path('/tmp/ref_check.faa'))
"
grep -c "^>2903220_liq80xsp_MAT_combined\|^>2903222_liq146xsp_MAT_combined" /tmp/ref_check.faa
```

Expected: 8 (4 genes × 2 records).

- [ ] **Step 6: Run the full suite, commit**

```bash
pixi run pytest -v
git add db/Ascomycota/Teloschistales/2903220_liq80xsp_MAT_combined db/Ascomycota/Teloschistales/2903222_liq146xsp_MAT_combined
git commit -m "feat: derive real protein sequences for the 2 Xanthoria MAG records from their own curated coordinates"
```

---

### Task 3: Re-derive Taphrina deformans' matMc/matPi coordinates via a real tblastn search, then derive their protein sequences (data/analysis step)

**Files:** none new — this task performs a real bioinformatic search and updates one curated record.

**Background:** `db/Ascomycota/Taphrinales/5011_unknown-1_PM_combined` (*Taphrina deformans*, `PM` family) is `validation.status: needs_review` — a previous curation error (citing coordinates from a chronologically-impossible source) was already caught and corrected to `coordinate_provenance: not_available`, with the record's own note recommending "a tblastn/exonerate search of the matMc/matPi reference proteins against RHGF01000006.1" as the real fix. Real, usable reference proteins already exist in this DB, curated from the SAME source paper (Almeida et al. 2015) for the related genus *Pneumocystis*: `db/Ascomycota/Pneumocystidales/42068_ru7_PM_combined` (`matPi: ncbi_protein:XP_018229374.1`, `matMc: ncbi_protein:XP_018229377.1`) and `db/Ascomycota/Pneumocystidales/4754_b80_PM_combined` (`matPi: ncbi_protein:XP_018224792.1`, `matMc: ncbi_protein:XP_018224795.1`).

- [ ] **Step 1: Fetch the 4 real reference proteins and the real Taphrina scaffold**

```bash
pixi run python -c "
from pathlib import Path
from MATPredict.config import MatpredictConfig
from MATPredict.db.cli import _make_clients
config = MatpredictConfig.from_env(repo_root=Path.cwd())
ncbi, _uniprot = _make_clients(config)

refs = {
    'matPi_ru7': 'XP_018229374.1', 'matMc_ru7': 'XP_018229377.1',
    'matPi_b80': 'XP_018224792.1', 'matMc_b80': 'XP_018224795.1',
}
with open('/tmp/pm_refs.faa', 'w') as f:
    for name, acc in refs.items():
        seq = ncbi.fetch_protein_sequence(acc)
        f.write(f'>{name}|{acc}\n{seq}\n')
        print(name, acc, len(seq), 'aa')
"
```

- [ ] **Step 2: Determine RHGF01000006.1's real length, fetch the full scaffold**

```bash
pixi run python -c "
import urllib.request
resp = urllib.request.urlopen(
    'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=RHGF01000006.1&retmode=json',
    timeout=15,
).read().decode()
print(resp)
"
```

Read the real `slen` (sequence length) field from the response, then fetch the whole scaffold:

```bash
pixi run python -c "
from pathlib import Path
from MATPredict.config import MatpredictConfig
from MATPredict.db.cli import _make_clients
config = MatpredictConfig.from_env(repo_root=Path.cwd())
ncbi, _uniprot = _make_clients(config)
seq = ncbi.fetch_nucleotide_sequence('RHGF01000006.1', 1, <REAL_LENGTH_FROM_STEP_2>, strand=None)
with open('/tmp/rhgf6.fa', 'w') as f:
    f.write('>RHGF01000006.1\n' + seq + '\n')
print(len(seq))
"
```

(Substitute the real length you found; do not guess it.)

- [ ] **Step 3: Run the real tblastn search**

```bash
pixi run makeblastdb -in /tmp/rhgf6.fa -dbtype nucl -out /tmp/rhgf6_db
pixi run tblastn -query /tmp/pm_refs.faa -db /tmp/rhgf6_db -seg no \
  -outfmt "6 qseqid sseqid pident length sstart send evalue bitscore" \
  -out /tmp/tblastn_results.tsv
cat /tmp/tblastn_results.tsv
```

- [ ] **Step 4: Interpret the results and determine real coordinates**

For each of `matPi`/`matMc`, identify the best-scoring hit (by bitscore, requiring a sane e-value — this project's own `search.py` uses `-seg no` with no other pre-filter, so apply your own judgment on what counts as a credible hit given these are cross-genus, likely fast-evolving MAT genes — a hit in the tens-of-percent identity range with a clear best-scoring peak over other candidates is expected, not a near-exact match). Convert the tblastn subject start/end (which may be reported start>end for a minus-strand hit) into this project's 1-based fully-closed, `start<end`-with-a-separate-`strand`-field convention (matching every other record's `genes[].start/end/strand` — read a couple of real existing records for the exact convention if unsure).

**If no credible hit is found for one or both genes**, do not fabricate a result — update the record to state plainly that the tblastn re-derivation was attempted and did not produce a credible hit, keep `coordinate_provenance: not_available` and `validation.status: needs_review`, and stop this task here (do not proceed to Steps 5-7) — report this as a real, honest outcome to the human reviewer rather than forcing a low-confidence coordinate into the record.

- [ ] **Step 5: If a credible hit was found for both genes, update the record**

Edit `db/Ascomycota/Taphrinales/5011_unknown-1_PM_combined/metadata.yaml`:
- Update `genes[].start`/`end`/`strand` for `matPi`/`matMc` to the real, re-derived coordinates.
- Update `locus.core.coordinate_provenance` from `not_available` to a value reflecting real independent re-derivation (check `db/_schema/metadata.schema.yaml`'s enum for `coordinate_provenance` — use whichever existing value means "independently re-derived by this project," e.g. something like `independently_derived` — read the schema's actual enum list yourself rather than guessing a value that doesn't exist in it).
- Update `locus.core.excluded_from_coordinate_benchmark` to `false` (it was `true` specifically because the old coordinates were unverified; that condition no longer holds once real coordinates are established).
- Append a new paragraph to `locus.core.definition_note` documenting exactly what was done: the 4 real Pneumocystis reference accessions used, the real tblastn command and parameters, the real hit coordinates/identity/e-value for each gene, and today's date — matching the existing note's own level of rigor and citation density.
- Update `evidence.boundaries` from `tier: 2, experimental_method: UNVERIFIED...` to a tier/method reflecting real bioinformatic re-derivation (this project's own tier table treats genome-annotation/homology evidence as tier 2 — keep it tier 2, but replace the "UNVERIFIED" method text with a real description of the tblastn re-derivation performed).
- Update `validation.status` from `needs_review` to `accepted` — but this is a curation acceptance decision; per this project's standing rule ("Never self-accept a proposed record — always a human decision"), do NOT change `validation.status` yourself. Leave it as `needs_review` and flag this explicitly to the human reviewer as ready for their accept/reject call, with the real tblastn evidence summarized for them.

- [ ] **Step 6: Once the coordinates are updated (regardless of the human's later accept/reject decision on validation.status), derive the protein sequences using Task 1's new capability**

```bash
pixi run python -c "
from pathlib import Path
from MATPredict.db.cli import build_gff_for_record, _make_clients
from MATPredict.config import MatpredictConfig
config = MatpredictConfig.from_env(repo_root=Path.cwd())
ncbi, uniprot = _make_clients(config)
build_gff_for_record(config.db_root, 'Ascomycota', 'Taphrinales', '5011_unknown-1_PM_combined', ncbi, uniprot)
"
cat db/Ascomycota/Taphrinales/5011_unknown-1_PM_combined/proteins.faa
```

Confirm 2 real, plausible protein sequences (matPi, matMc), and sanity-check each against the reference proteins fetched in Step 1 (a same-genus-family-level BLAST identity is expected — not exact identity, since these are different species/genera, but the two sequences should be broadly comparable in length/composition to their Pneumocystis counterparts, not wildly different).

- [ ] **Step 7: Run the full suite, commit**

```bash
pixi run pytest -v
git add db/Ascomycota/Taphrinales/5011_unknown-1_PM_combined
git commit -m "fix: re-derive Taphrina deformans matMc/matPi coordinates via real tblastn search against Pneumocystis references"
```

---

### Task 4: Wire the evidence-floor gate into real entry points; fix `min_hits` semantics

**Files:**
- Modify: `src/MATPredict/detect/cli.py` (`_cmd_detect`)
- Modify: `src/MATPredict/detect/batch_runner.py` (`run_batch`)
- Modify: `src/MATPredict/detect/pipeline.py` (`_families_meeting_evidence_floor`, `_write_evidence_diagnostics`)
- Test: `tests/detect/test_pipeline_evidence_floor.py`, `tests/test_cli_smoke.py` (or wherever `_cmd_detect`'s argparse wiring is already tested — check first)

**Interfaces:**
- Produces: `matpredict detect --evidence-diagnostics <path> [--min-hits N] [--min-identity F] [--require-core-role]` new optional CLI flags; `run_batch`'s existing per-genome `run_pipeline` call gains `evidence_floor`/`evidence_diagnostics_path` (written to `genome_out_dir / "evidence_diagnostics.jsonl"` per genome, resolving the final review's M-1 finding about missing genome identity in diagnostics rows for free — one file per genome directory needs no genome id in the row).
- Fixes: `_families_meeting_evidence_floor`'s `min_hits` now counts **distinct genes** (`len({h.gene_name for h in own_hits})`), matching the field's name and its own existing test's name (`test_min_hits_floor_admits_a_family_with_enough_distinct_genes`) rather than raw HSP count. `_write_evidence_diagnostics` gains a new `hit_count` field (the raw HSP count, kept separately) so a future calibration pass can see both numbers rather than losing the raw-HSP-count information the old (buggy) `min_hits` semantics happened to expose.

**Background:** `search.py`'s tblastn loop appends one `SearchHit` per HSP line with no per-gene deduplication, so a single real gene matched by N curated reference records (routine after Task 2's earlier `Ascomycota:MAT` backfill, which took that family from 12 to 80 reference sequences) produces N hits. The prior plan's final review found `min_hits` was gating on this raw, curation-density-dependent count while the diagnostics `gene_count` field (meant to help calibrate it) counted distinct genes — a latent mismatch that would make a threshold "calibrated" from `gene_count` wrong when plugged into `min_hits`.

- [ ] **Step 1: Write the failing test for the `min_hits` fix**

```python
# append to tests/detect/test_pipeline_evidence_floor.py
def test_min_hits_now_counts_distinct_genes_not_raw_hsps():
    family = _family("Ascomycota", "MAT")
    # Two HSPs for the SAME gene (e.g. two curated reference records both hit
    # by tblastn) -- this is exactly the case search.py's own lack of
    # per-gene dedup produces, and it must NOT count as 2 toward min_hits.
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1", identity=40.0),
        _hit(family.key, gene_name="G1", identity=42.0),
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_hits=2))

    assert result == []  # only 1 DISTINCT gene, even though there are 2 hit objects


def test_min_hits_admits_when_distinct_gene_count_clears_the_floor():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1", identity=40.0),
        _hit(family.key, gene_name="G1", identity=42.0),  # 2nd HSP, same gene
        _hit(family.key, gene_name="G2", identity=38.0),  # a genuinely different gene
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_hits=2))

    assert [f.key for f in result] == [family.key]  # 2 distinct genes (G1, G2)
```

- [ ] **Step 2: Run to verify it fails**

```bash
pixi run pytest tests/detect/test_pipeline_evidence_floor.py -v
```

Expected: the first new test FAILS (old code counts 2 raw hits, admits the family when it shouldn't).

- [ ] **Step 3: Fix `_families_meeting_evidence_floor`**

Change the `min_hits` check in `pipeline.py`:

```python
    admitted = []
    for family in families:
        own_hits = [h for h in cluster.hits if h.family_key == family.key]
        distinct_genes = {h.gene_name for h in own_hits}
        if len(distinct_genes) < floor.min_hits:
            continue
```

(the rest of the function is unchanged). Update `EvidenceFloor.min_hits`'s docstring/field-level comment to say explicitly that it counts distinct genes, not raw hit objects, and why (cross-reference `search.py`'s lack of per-HSP-per-reference-record dedup).

- [ ] **Step 4: Run tests to verify they pass**

```bash
pixi run pytest tests/detect/test_pipeline_evidence_floor.py -v
```

Expected: PASS, all tests (including the pre-existing `test_min_hits_floor_admits_a_family_with_enough_distinct_genes`/`test_min_hits_floor_rejects_a_single_hit_family`, which used distinct gene names anyway so are unaffected by this fix).

- [ ] **Step 5: Add the `hit_count` field to diagnostics**

In `_write_evidence_diagnostics`, add `"hit_count": len(own_hits)` alongside the existing `"gene_count": len({h.gene_name for h in own_hits})` field. Update its test (`test_evidence_diagnostics_written_when_path_given`) to assert on the new field too — construct a cluster with 2 hits for the same gene and confirm `hit_count == 2` while `gene_count == 1`.

- [ ] **Step 6: Run tests to verify they pass**

```bash
pixi run pytest tests/detect/test_pipeline_evidence_floor.py -v
```

- [ ] **Step 7: Wire the CLI flags into `_cmd_detect`**

Read `src/MATPredict/detect/cli.py`'s current `register_subcommands`/argument-registration for the `detect` subparser (find the `detect.add_argument(...)` block) and add:

```python
    detect.add_argument("--evidence-diagnostics", required=False)
    detect.add_argument("--min-hits", type=int, default=1)
    detect.add_argument("--min-identity", type=float, default=None)
    detect.add_argument("--require-core-role", action="store_true")
```

In `_cmd_detect`, construct the floor and pass both new parameters through to `run_pipeline`:

```python
    from MATPredict.detect.pipeline import EvidenceFloor  # add to the top-of-file imports instead if cleaner -- check existing import grouping style first

    evidence_floor = EvidenceFloor(
        min_hits=args.min_hits, min_identity=args.min_identity,
        require_core_role=args.require_core_role,
    )
    outcome = run_pipeline(
        genome_fasta=Path(args.genome),
        proteome_fasta=Path(args.proteins) if args.proteins else None,
        taxid=args.taxid,
        db_root=config.db_root,
        reference_fasta=reference_fasta,
        evidence_floor=evidence_floor,
        evidence_diagnostics_path=Path(args.evidence_diagnostics) if args.evidence_diagnostics else None,
    )
```

- [ ] **Step 8: Wire diagnostics into `batch_runner.run_batch`**

In `src/MATPredict/detect/batch_runner.py`'s `run_batch`, add `evidence_diagnostics_path=genome_out_dir / "evidence_diagnostics.jsonl"` to its existing `run_pipeline(...)` call (the block already shown in this plan's Background section) — do NOT add an `evidence_floor` override here; leave it at `run_pipeline`'s own default (`EvidenceFloor()`, still a no-op) so a real batch run collects diagnostics without changing detection behavior. `genome_out_dir` already exists at this point in the function (it's what `write_detection_gff3`/`write_detection_report` already write into) — reuse it, don't compute a new path.

- [ ] **Step 9: Write/extend a CLI smoke test**

Find how `_cmd_detect`'s existing argparse wiring is tested (check `tests/test_cli_smoke.py` first, per the file list above) and add a test confirming the new flags parse correctly and are threaded through — following whatever mocking convention the existing `test_detect_subcommand_registered`-style tests already use in that file (read 1-2 of them first) rather than inventing a new test style.

- [ ] **Step 10: Run the full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/detect/cli.py src/MATPredict/detect/batch_runner.py src/MATPredict/detect/pipeline.py tests/detect/test_pipeline_evidence_floor.py tests/test_cli_smoke.py
git commit -m "feat: wire evidence-floor diagnostics into detect CLI and batch runner; fix min_hits to count distinct genes"
```

---

### Task 5: Investigate genome 5501's weak real-family detection (data/analysis step)

**Files:** none new — this task re-runs detection against the real, already-acquired genome 5501 (*C. immitis*, `GCA_004115165.2`) with all of this session's fixes in place, and investigates why only 1 of ~10 core genes was found.

**Background:** In the original rollout, genome 5501's own real family (`Ascomycota:MAT`) detected only 1 gene (`APN2`, a flanking gene, at low confidence, matched against the WRONG reference record — a Diaporthales record, not this genome's own curated Onygenales records `5501_h538-4_MAT_MAT1-1`/`5501_rs_MAT_MAT1-2`). None of the core idiomorph genes (`MAT1-1-1/2/3/5`, `MAT1-2-1/10/4`, `SLA2`, `COX13`, `CIMG_00407`) were found. This was flagged as needing investigation "independent of the routing bug" — this task is that investigation.

- [ ] **Step 1: Re-run `matpredict detect` against genome 5501 with all current fixes**

```bash
pixi run matpredict detect \
  --genome <path to GCA_004115165.2's decompressed FASTA -- check Task 2 of the genome-scale-detection-rollout plan's real acquisition output, or re-run acquire_genomes for this one taxid; it resolves via the local BFD library> \
  --taxid 5501 \
  --out-dir /tmp/5501_recheck \
  --evidence-diagnostics /tmp/5501_recheck/evidence_diagnostics.jsonl
```

Confirm via the printed `families_attempted` count that routing now correctly narrows to `Ascomycota:MAT` alone (or very few families), not all 19 — this alone should already dramatically improve the signal-to-noise ratio versus the original run.

- [ ] **Step 2: Compare the new result against the curated ground truth**

Read `db/Ascomycota/Onygenales/5501_h538-4_MAT_MAT1-1/metadata.yaml` and `db/Ascomycota/Onygenales/5501_rs_MAT_MAT1-2/metadata.yaml`'s own gene list/coordinates/protein_accessions. For each core gene NOT found in the new run, investigate why using the actual raw evidence:
- Was the gene searched at all? Check the new run's `evidence_diagnostics.jsonl` for any row mentioning `Ascomycota:MAT` and whatever cluster(s) exist near the curated genes' real genomic locations.
- If a gene WAS searched but not admitted/detected, what was its real tblastn identity/coverage against the curated reference proteins? Genome 5501's real accession is `GCA_004115165.2` (strain `WA_211`, per this session's earlier ground-truth work) — a genuinely different isolate from either curated Onygenales record (`H538.4`/`RS`), already established and accepted as a real strain difference, not a bug (see Task 5 of the earlier genome-scale-detection-rollout plan's `match_ground_truth` work). A real, biologically-expected level of sequence divergence between strains could itself explain weak-but-real hits that don't clear detection thresholds — distinguish this from a genuine technical bug (wrong window, wrong reference record matched, a translation/frame error, or a real assembly gap/misassembly at this genome's own MAT locus) using the actual evidence, not assumption.
- Specifically investigate why the one detected hit (`APN2`) matched a Diaporthales reference record instead of one of this genome's own curated Onygenales records — read `search.py`'s `_attribute`/reference-selection logic to understand how a gene's best-matching reference record is chosen, and check whether the Onygenales records' own `APN2` sequences (from `5501_h538-4_MAT_MAT1-1`/`5501_rs_MAT_MAT1-2`, or check whether those records even declare an `APN2` gene at all — some MAT records may only declare core genes, not this particular flanking gene) were even in the search's candidate set at the time, or whether a real higher-identity match to a distant reference record legitimately outscored the "right" one.

- [ ] **Step 3: Write up the findings**

Document the real cause (or causes, if more than one factor contributes) in a short findings note — append to `docs/superpowers/plans/2026-09-18-genome-scale-detection-rollout-findings.md` as a new dated section (do not create a new file for this) with: whether this is (a) a real, expected consequence of strain-level sequence divergence, (b) a real reference-record-selection gap in the pipeline worth its own fix, (c) a real assembly-quality issue specific to this genome, or (d) something else — cite the actual evidence for whichever conclusion(s) you reach. Do not fix any pipeline bug found here inline as part of this task — list it as a new, separate triage item for a future plan, matching this project's established practice of separating investigation from fixing.

- [ ] **Step 4: Commit the findings update**

```bash
git add docs/superpowers/plans/2026-09-18-genome-scale-detection-rollout-findings.md
git commit -m "docs: investigate genome 5501's weak real-family detection"
```

---

### Task 6: Re-run the full 13-genome pilot rollout with all fixes in place (data/execution step)

**Files:** none new — this task executes Tasks 1-3 of the prior `2026-09-18-detection-rollout-fixes` plan's already-built tooling, now with this plan's Tasks 1-4 also landed, against the same confirmed 13-taxid pilot list.

- [ ] **Step 1: Re-acquire the 13 pilot genomes** (cheap — all resolve locally, no network download; see the original rollout plan's Task 2 for the exact taxid list: 5501, 199306, 162425, 746128, 5061, 5059, 5076, 27334, 36651, 5141, 5518, 5507, 510951).

- [ ] **Step 2: Plan and run the batch, with diagnostics enabled**

Use `scripts/run_detection_batch.py` (or call `batch_runner.run_batch` directly, per the original plan's own guidance on avoiding redundant per-genome setup) against all 13 genomes. Since Task 4 of THIS plan wired `evidence_diagnostics_path` into `run_batch` by default, no extra flag is needed — every genome's run will now produce its own `evidence_diagnostics.jsonl` for free, which is real, valuable calibration data this rollout previously could not collect (per the prior plan's own stated rationale for why the diagnostics mode was built in the first place).

Budget real time for this based on what Task 1 of the prior plan's routing fix should now deliver: most genomes should route to 1-2 families (not the 19-family exhaustive fallback that caused the original ~2h20m/genome cost), so the batch should now plausibly complete close to the original ~19-minute estimate — but confirm this empirically rather than assuming it; if any genome still takes dramatically longer than the ~96s/39Mb reference figure, investigate why before assuming it's fine, using the same live-diagnosis approach this session already used successfully once.

- [ ] **Step 3: Aggregate and score**

Run `aggregate_reports`/`matpredict detect rollout-summary` over the batch output, and `benchmark.match_ground_truth`/`score_self_consistency` for the 4 ground-truth-tier genomes, exactly as the prior plan's Task 6 did.

- [ ] **Step 4: Write the real findings**

Update `docs/superpowers/plans/2026-09-18-genome-scale-detection-rollout-findings.md` with a new dated section covering the real, complete 13-genome result: tier distribution, per-family detection rates, every anomaly, every ground-truth score obtained (now potentially including the 2 Xanthoria and Taphrina-adjacent families this plan's Tasks 1-3 made searchable for the first time, if any pilot genome happens to route to them), and an honest comparison against the original single-genome result (was the routing fix's real-world speedup as large as expected? did detection sensitivity improve for genome 5501 specifically, given Task 5's investigation?). Present a new triage list for whatever this full run surfaces — the whole point of running at this scale, same as the first time.

- [ ] **Step 5: Present the findings to the human reviewer** rather than closing this plan silently, per this project's established practice for every rollout run.

## Self-review notes (controller, at plan-writing time)

- **Spec coverage:** the user's 4-item sequence (source accessions → wire EvidenceFloor/resolve min_hits → investigate 5501 → kick off the pilot run) maps to Tasks 1-3, Task 4, Task 5, and Task 6 respectively, in the same order.
- **No placeholders in code steps:** Tasks 1 and 4's code is real, complete, and grounded in this session's actual current source (`validate.py`'s real `_independent_translation`, `pipeline.py`'s real current `_families_meeting_evidence_floor`, `cli.py`'s real current `_cmd_detect`/`build_gff_for_record`, `batch_runner.py`'s real current `run_pipeline` call site) — no invented signatures. Tasks 2, 3, 5, 6 are data/analysis steps (matching the established style of this session's earlier plans' own data-step tasks) with concrete, executable commands rather than vague direction, and explicit "stop and report honestly" instructions for the genuinely open-ended parts (Task 3's tblastn interpretation, Task 5's root-cause investigation) rather than a placeholder telling the implementer to "figure it out."
- **Type/signature consistency:** Task 1's `build_gff_for_record` change is consumed unchanged by Task 2/3 (same function, same call signature). Task 4's `EvidenceFloor`/`evidence_diagnostics_path` parameters are the exact ones already defined by the prior plan's Task 3 — this plan only adds call sites and one real bug fix, it does not redefine the interface.
- **Task dependencies:** Task 1 must land before Tasks 2/3 (they use its new fallback path). Task 4 should land before Task 6 (so the pilot run collects real diagnostics and uses the corrected `min_hits` semantics) but is independent of Tasks 1-3. Task 5 can run any time after the `taxonomic_scope` fix already merged in the prior plan; running it before Task 6 (per the user's stated order) lets its findings inform what to watch for in the full pilot run.
