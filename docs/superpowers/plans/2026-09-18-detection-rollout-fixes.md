# Detection Rollout Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the three real, root-caused gaps the 2026-09-18 genome-scale detection rollout surfaced: (1) a stale `taxonomic_scope` value that puts 11 of 13 pilot genomes into `family_registry.route()`'s exhaustive-fallback path, (2) a `proteins.faa` materialization gap that leaves 2 families completely unsearchable and starves the `Ascomycota:MAT` family's own reference set, and (3) an unfiltered Stage 1 → Stage 2 handoff that lets a single spurious tblastn hit trigger dozens of expensive polish subprocess calls.

**Architecture:** Three independent fixes, each landing as its own task with its own tests: (1) a new, reusable `scope_audit` module + CLI subcommand that measures and recommends `taxonomic_scope` fixes, applied to `db/Ascomycota/order.yml`; (2) a `backfill-gff` CLI subcommand that regenerates `proteins.faa`/`locus.gff3`/`locus.gbk` for every accepted record missing them, run for real against the live `db/` tree; (3) a configurable, default-no-op evidence-floor gate on the Stage 1 → Stage 2 handoff in `pipeline.py`, shipped with a diagnostics mode so its real threshold can be calibrated from data once (1) has landed — this plan does NOT invent an untested threshold number, because real data from the rollout shows a naive identity or hit-count cutoff would have filtered the one genuinely correct signal `matpredict detect` found while keeping several of the highest-confidence spurious hits (see Task 3's own notes).

**Tech Stack:** Python, this project's existing `family_registry`/`taxonomy`/`ncbi_client`/`gff_export` modules, `pixi run pytest`.

**Spec:** No separate spec doc — this plan implements the triage list from `docs/superpowers/plans/2026-09-18-genome-scale-detection-rollout-findings.md` (read that file for full root-cause detail before Task 1) plus the reorganized recommendation from this session's screening-strategy evaluation (not a saved doc; summarized in each task's own notes below).

## Global Constraints

- Never build a second NCBI/UniProt client — reuse `src/MATPredict/db/cli.py`'s `_make_clients`/`_client_for` and `src/MATPredict/db/ncbi_client.py`'s `NcbiClient` for any live fetch.
- Never let one bad record/genome/cluster abort a batch — every loop introduced in this plan isolates failures per-item (log + continue), matching this project's established convention (see the genome-acquisition and batch-runner fixes from the immediately preceding rollout plan).
- No changes to `run_pipeline`'s existing behavior for callers that don't opt into the new Task 3 parameter — its default must reproduce today's exact behavior.
- Live-verify any new taxid/taxonomy claim against real NCBI data before writing it into `order.yml` or a commit message — this project's standing rule against unverified claims applies to curated data, not just code.
- `db/candidates/` records are never authoritative and must stay excluded from every audit/backfill this plan adds, matching the existing exclusion convention in `family_registry.load_record_families`, `benchmark._load_records`, and `reference_fasta.build_reference_fasta`.

---

### Task 1: Audit and fix `taxonomic_scope` data across all three `order.yml` files

**Files:**
- Create: `src/MATPredict/detect/scope_audit.py`
- Modify: `src/MATPredict/detect/cli.py` (register `matpredict detect audit-scope`)
- Modify: `db/Ascomycota/order.yml` (fix the `MAT` locus's `taxonomic_scope`)
- Test: `tests/detect/test_scope_audit.py`

**Interfaces:**
- Produces: `record_taxids_by_family(db_root: Path) -> dict[FamilyKey, list[int]]`, `deepest_common_ancestor(taxids: list[int], lineage_taxids_resolver: Callable[[int], list[int]]) -> int | None`, `ScopeAuditResult` (frozen dataclass: `family_key: FamilyKey`, `total_records: int`, `uncovered_taxids: list[int]`, `recommended_scope_taxid: int | None`), `audit_scope(families: list[Family], record_taxids: dict[FamilyKey, list[int]], lineage_taxids_resolver: Callable[[int], list[int]] = default_lineage_taxids) -> list[ScopeAuditResult]`. This is the same live-lineage-overlap check `route()` itself uses, run offline against every curated record instead of one query taxid at a time — no genome data or detection run required.

**Background this task fixes:** `db/Ascomycota/order.yml`'s `MAT` locus declares `taxonomic_scope: [222544, 5180, 28548, 40559]` (four Helotiales-adjacent taxa), even though its own code comment says these genes are "shared across Pezizomycotina" and its own curated records span at least 7 orders (Helotiales, Onygenales, Eurotiales, Hypocreales, Sordariales, Teloschistales, Lecanorales). `family_registry.route()` therefore cannot match most of this family's own genomes by direct membership or lineage overlap, and falls through to its documented exhaustive fallback — searching all 19 curated families against a genome that should search 1. This was independently found and quantified twice this session: a live audit found 17 of 22 curated `Ascomycota:MAT` records (77%) fall outside their own family's declared scope, and 11 of the rollout's 13 pilot genomes (85%) hit the exhaustive fallback as a direct result.

- [ ] **Step 1: Write the failing test for `deepest_common_ancestor`**

```python
# tests/detect/test_scope_audit.py
from __future__ import annotations

from MATPredict.detect.scope_audit import (
    ScopeAuditResult,
    audit_scope,
    deepest_common_ancestor,
    record_taxids_by_family,
)
from MATPredict.detect.family_registry import Family, FamilyKey


def _family(phylum, name, scope, genes=None):
    return Family(
        key=FamilyKey(phylum, name),
        vocabulary_type="enum",
        idiomorph_values=["a", "alpha"],
        idiomorph_pattern=None,
        genes=genes or [{"name": "STE3", "role": "core_MAT"}],
        taxonomic_scope=scope,
    )


# NCBI-shaped fake lineages: each list is root-to-parent order, most-specific last,
# matching the real order `default_lineage_taxids` returns (verified live this
# session: taxid 5501's lineage ends ..., 147545 (Eurotiomycetes), 451871, 33183,
# 33184 (Onygenales), 5500 (Onygenaceae) -- broad to narrow).
_FAKE_LINEAGES = {
    111: [1, 10, 100],          # e.g. species 111 under genus 100 under order 10
    112: [1, 10, 100],          # same genus/order as 111
    113: [1, 10, 200],          # same order (10) but a different genus (200)
    999: [1, 20, 300],          # unrelated order entirely
}


def _fake_resolver(taxid: int) -> list[int]:
    return _FAKE_LINEAGES[taxid]


def test_deepest_common_ancestor_of_two_taxids_sharing_a_genus():
    assert deepest_common_ancestor([111, 112], _fake_resolver) == 100


def test_deepest_common_ancestor_falls_back_to_shared_order():
    assert deepest_common_ancestor([111, 113], _fake_resolver) == 10


def test_deepest_common_ancestor_of_a_single_taxid_is_itself():
    assert deepest_common_ancestor([111], _fake_resolver) == 111


def test_deepest_common_ancestor_returns_none_for_empty_list():
    assert deepest_common_ancestor([], _fake_resolver) is None


def test_deepest_common_ancestor_returns_none_when_genuinely_disjoint():
    # 111's lineage is [1, 10, 100]; 999's is [1, 20, 300] -- only root taxid 1
    # is shared, so the deepest common ancestor is 1, not None: two real fungi
    # always share at least a root/kingdom-level ancestor. None is reserved for
    # an empty input, not for "very distantly related."
    assert deepest_common_ancestor([111, 999], _fake_resolver) == 1
```

- [ ] **Step 2: Run it to verify it fails**

```bash
pixi run pytest tests/detect/test_scope_audit.py -v
```

Expected: FAIL with `ModuleNotFoundError: No module named 'MATPredict.detect.scope_audit'`.

- [ ] **Step 3: Write the failing test for `audit_scope`**

Append to the same file:

```python
def test_audit_scope_flags_uncovered_records_and_recommends_a_fix():
    # Family's own scope [999] covers neither of its two real records (111, 112),
    # whose real common ancestor is 100.
    family = _family("Ascomycota", "MAT", scope=[999])
    record_taxids = {family.key: [111, 112]}

    results = audit_scope([family], record_taxids, lineage_taxids_resolver=_fake_resolver)

    assert results == [
        ScopeAuditResult(
            family_key=family.key, total_records=2,
            uncovered_taxids=[111, 112], recommended_scope_taxid=100,
        )
    ]


def test_audit_scope_reports_no_uncovered_taxids_when_scope_already_correct():
    family = _family("Ascomycota", "MATsc", scope=[100])  # 100 covers both via lineage
    record_taxids = {family.key: [111, 112]}

    results = audit_scope([family], record_taxids, lineage_taxids_resolver=_fake_resolver)

    assert results == [
        ScopeAuditResult(
            family_key=family.key, total_records=2,
            uncovered_taxids=[], recommended_scope_taxid=None,
        )
    ]


def test_audit_scope_direct_membership_counts_as_covered():
    # 111 is listed literally in scope -- must not require lineage resolution
    # (mirrors route()'s own direct-membership short-circuit).
    def unused_resolver(taxid):
        raise AssertionError("direct membership match must short-circuit before any lineage lookup")

    family = _family("Ascomycota", "MAT", scope=[111])
    record_taxids = {family.key: [111]}

    results = audit_scope([family], record_taxids, lineage_taxids_resolver=unused_resolver)

    assert results[0].uncovered_taxids == []
```

- [ ] **Step 4: Run tests to verify they fail**

```bash
pixi run pytest tests/detect/test_scope_audit.py -v
```

Expected: FAIL (same import error, or `NameError`/`AttributeError` once the module exists but the functions don't).

- [ ] **Step 5: Implement `scope_audit.py`**

```python
"""Offline audit of every curated record's taxid against its own family's
declared taxonomic_scope -- the same direct-membership/lineage-overlap check
family_registry.route() uses at detection time, run against the whole curated
DB instead of one query taxid, so a stale scope value can be found and fixed
without running the detection pipeline at all."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import yaml

from MATPredict.db.taxonomy import default_lineage_taxids
from MATPredict.detect.family_registry import Family, FamilyKey


def record_taxids_by_family(db_root: Path) -> dict[FamilyKey, list[int]]:
    """Every accepted (non-candidate) record's own taxid, grouped by the
    FamilyKey it belongs to. Mirrors family_registry.load_record_families's
    glob and candidates/ exclusion, but collects taxids instead of an index."""
    by_family: dict[FamilyKey, list[int]] = {}
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        if meta_path.relative_to(db_root).parts[0] == "candidates":
            continue
        doc = yaml.safe_load(meta_path.read_text())
        locus_name = (doc.get("mating_type") or {}).get("locus_name")
        taxid = (doc.get("taxonomy") or {}).get("taxid")
        if not locus_name or taxid is None:
            continue
        key = FamilyKey(meta_path.parents[2].name, locus_name)
        by_family.setdefault(key, []).append(taxid)
    return by_family


def deepest_common_ancestor(
    taxids: list[int], lineage_taxids_resolver: Callable[[int], list[int]]
) -> int | None:
    """The most specific taxid shared by every lineage in `taxids`, or the
    single taxid itself when only one is given, or None for an empty list.

    `lineage_taxids_resolver` (normally default_lineage_taxids) returns each
    taxid's ancestor chain broad-to-narrow (root first), never including the
    taxid itself -- so each taxid's own full lineage path, for this
    comparison, is `[*ancestors, taxid]`. Two real taxa always share at least
    a root-level ancestor, so this returns None only for an empty input, not
    for "very distantly related" -- callers that want to reject a too-broad
    recommendation should check the returned taxid's own rank/scientific name
    (e.g. via a live NCBI esummary lookup) before writing it into order.yml,
    the same live-verification discipline this project already requires
    everywhere else.
    """
    if not taxids:
        return None
    paths = [[*lineage_taxids_resolver(t), t] for t in taxids]
    # A taxonomic lineage is a genuine tree path (each taxid has exactly one
    # parent, never rejoining a sibling branch), so the deepest ancestor
    # shared by every path is exactly its longest common PREFIX -- position
    # by position, stop at the first position where the paths disagree.
    # zip(*paths) already stops at the shortest path's length.
    common: list[int] = []
    for level in zip(*paths):
        if len(set(level)) != 1:
            break
        common.append(level[0])
    return common[-1] if common else None


@dataclass(frozen=True)
class ScopeAuditResult:
    family_key: FamilyKey
    total_records: int
    uncovered_taxids: list[int]
    recommended_scope_taxid: int | None


def _covered(taxid: int, scope: list[int], ancestors: list[int]) -> bool:
    return taxid in scope or bool(set(ancestors) & set(scope))


def audit_scope(
    families: list[Family],
    record_taxids: dict[FamilyKey, list[int]],
    lineage_taxids_resolver: Callable[[int], list[int]] = default_lineage_taxids,
) -> list[ScopeAuditResult]:
    """For every family with at least one curated record, check whether
    route()'s own direct-membership-or-lineage-overlap test would cover each
    record's taxid, and recommend a single replacement scope taxid (the
    deepest common ancestor of every UNCOVERED record's taxid) when it would
    not. Families with zero uncovered records get `recommended_scope_taxid =
    None` -- there is nothing to fix."""
    results = []
    for family in families:
        taxids = record_taxids.get(family.key, [])
        uncovered = []
        for taxid in taxids:
            ancestors = [] if taxid in family.taxonomic_scope else lineage_taxids_resolver(taxid)
            if not _covered(taxid, family.taxonomic_scope, ancestors):
                uncovered.append(taxid)
        recommended = (
            deepest_common_ancestor(uncovered, lineage_taxids_resolver) if uncovered else None
        )
        results.append(ScopeAuditResult(
            family_key=family.key, total_records=len(taxids),
            uncovered_taxids=uncovered, recommended_scope_taxid=recommended,
        ))
    return results
```

- [ ] **Step 6: Run tests to verify they pass**

```bash
pixi run pytest tests/detect/test_scope_audit.py -v
```

Expected: PASS, all 7 tests.

- [ ] **Step 7: Wire a `matpredict detect audit-scope` CLI subcommand**

Read `src/MATPredict/detect/cli.py`'s existing subcommand registration (the `rollout-summary` subcommand added in the previous rollout plan is the closest precedent — follow its exact style) and add:

```python
def _cmd_audit_scope(args: argparse.Namespace) -> int:
    db_root = Path(args.db_root)
    families = load_all_families(db_root)
    record_taxids = record_taxids_by_family(db_root)
    results = audit_scope(families, record_taxids)
    any_uncovered = False
    for result in sorted(results, key=lambda r: (-len(r.uncovered_taxids), r.family_key.phylum, r.family_key.locus_name)):
        if not result.uncovered_taxids:
            continue
        any_uncovered = True
        pct = 100.0 * len(result.uncovered_taxids) / result.total_records
        print(
            f"{result.family_key.phylum}:{result.family_key.locus_name}: "
            f"{len(result.uncovered_taxids)}/{result.total_records} records ({pct:.0f}%) "
            f"uncovered by taxonomic_scope -- recommend {result.recommended_scope_taxid}"
        )
    if not any_uncovered:
        print("every family's own curated records are covered by its taxonomic_scope")
    return 0
```

Add the corresponding `action.add_parser("audit-scope")` block with a `--db-root` argument (`default="db"`), following the exact `add_argument`/`set_defaults(func=...)` pattern the `rollout-summary` subcommand already uses in this same file. Import `record_taxids_by_family`, `audit_scope`, and `load_all_families` at the top of `cli.py` alongside its existing `MATPredict.detect.*` imports.

- [ ] **Step 8: Run the audit for real against the live `db/` tree**

```bash
pixi run matpredict detect audit-scope --db-root db
```

Confirm the real output reproduces (in substance, not necessarily exact counts, since this is a fresh independent implementation) the same finding this session already made twice: `Ascomycota:MAT` is the only family with a non-trivial uncovered fraction, and every other family in all three `order.yml` files reports 0 uncovered.

- [ ] **Step 9: Live-verify the tool's recommended taxid before writing it into `order.yml`**

```bash
pixi run python -c "
import urllib.request
resp = urllib.request.urlopen(
    f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=<RECOMMENDED_TAXID>&retmode=json',
    timeout=15,
).read().decode()
print(resp)
"
```

Confirm the recommended taxid resolves to a real, sensible fungal taxon (this session already live-confirmed taxid `147538` = *Pezizomycotina*, "filamentous ascomycetes", rank `subphylum` — matching `order.yml`'s own pre-existing code comment that these genes are "shared across Pezizomycotina"; if Step 8's tool recommends a different but still-sensible taxid due to a slightly different common-ancestor computation, prefer the tool's real output over this pre-computed value, but confirm it the same way before using it).

- [ ] **Step 10: Apply the fix to `db/Ascomycota/order.yml`**

Change the `MAT` locus's `taxonomic_scope` line from:

```yaml
taxonomic_scope: [222544, 5180, 28548, 40559]
```

to:

```yaml
taxonomic_scope: [147538]  # Pezizomycotina -- live-confirmed via NCBI esummary
    # 2026-XX-XX; covers every order this family's own curated records span
    # (Helotiales, Onygenales, Eurotiales, Hypocreales, Sordariales,
    # Teloschistales, Lecanorales). KNOWN, ACCEPTED TRADEOFF: this also
    # covers Pezizales, so a Tuber/Pezizales genome will now route to BOTH
    # this family AND Ascomycota:MATtub (Pezizales' own dedicated, narrowly-
    # scoped family) -- 2 families searched instead of 1, still vastly
    # cheaper than the 19-family exhaustive fallback this replaces. route()'s
    # OR-based multi-family matching already tolerates this; see
    # docs/superpowers/plans/2026-09-18-genome-scale-detection-rollout-findings.md
    # for the audit that found this value stale.
```

(Use the actual date and the actual recommended taxid from Step 8/9's real output if it differs from 147538.)

- [ ] **Step 11: Re-run the audit to confirm the fix**

```bash
pixi run matpredict detect audit-scope --db-root db
```

Expected: `every family's own curated records are covered by its taxonomic_scope`.

- [ ] **Step 12: Run the full suite, commit**

```bash
pixi run pytest -v
git add src/MATPredict/detect/scope_audit.py src/MATPredict/detect/cli.py tests/detect/test_scope_audit.py db/Ascomycota/order.yml
git commit -m "fix: correct stale Ascomycota:MAT taxonomic_scope; add reusable scope audit"
```

---

### Task 2: Backfill missing `proteins.faa` for every accepted record

**Files:**
- Modify: `src/MATPredict/db/cli.py` (extract a reusable `build_gff_for_record` helper from `_cmd_build_gff`'s body; add `matpredict curate-db backfill-gff`)
- Test: `tests/db/test_backfill_gff.py`

**Interfaces:**
- Produces: `find_records_missing_proteins_faa(db_root: Path) -> list[tuple[str, str, str]]` (each tuple is `(phylum, order_or_family, record_id)`), `build_gff_for_record(db_root: Path, phylum: str, order_or_family: str, record_id: str, ncbi: NcbiClient, uniprot: UniprotClient) -> None` (the exact logic `_cmd_build_gff` already runs per-record, extracted so both the single-record command and the new batch command call the same code), `backfill_missing_proteins_faa(db_root: Path, ncbi: NcbiClient, uniprot: UniprotClient) -> tuple[list[tuple[str, str, str]], list[tuple[tuple[str, str, str], str]]]` (returns `(succeeded, failed)`, where `failed` pairs each record identifier with the exception's string message — one record's live-fetch failure must never abort the batch).

**Background this task fixes:** `db/**/*/proteins.faa` is generated by `matpredict curate-db build-gff`, which fetches every gene's curated protein sequence live from NCBI/UniProt and writes `locus.gff3`, `locus.gbk`, and `proteins.faa`. 31 of 59 accepted records (all of them under `db/Ascomycota/`, matching this session's newly-added Eurotiales/Onygenales/Hypocreales/Sordariales/Helotiales/Teloschistales/Lecanorales/Pezizales records plus one older Saccharomycetales record) have never had this command run, so they contribute zero sequences to `reference_fasta.build_reference_fasta`'s output -- the concatenated FASTA every search actually uses. Two consequences, both measured this session: `Ascomycota:MATtub` and `Ascomycota:MTL` are completely unsearchable (0 of their curated sequences reach the reference FASTA), and `Ascomycota:MAT` -- the family Task 1 just fixed the routing for -- searches with only 12 of a possible ~80 curated sequences.

- [ ] **Step 1: Write the failing test for `find_records_missing_proteins_faa`**

```python
# tests/db/test_backfill_gff.py
from __future__ import annotations

import yaml

from MATPredict.db.cli import build_gff_for_record, find_records_missing_proteins_faa


def _write_record(db_root, phylum, order_or_family, record_id, genes=None):
    record_dir = db_root / phylum / order_or_family / record_id
    record_dir.mkdir(parents=True)
    metadata = {
        "record_id": record_id,
        "genes": genes if genes is not None else [
            {"gene_index": 0, "name": "G1", "present": True, "protein_accession": "ncbi_protein:ABC1.1"},
        ],
    }
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(metadata))
    return record_dir


def test_find_records_missing_proteins_faa_skips_records_that_already_have_one(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, "Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    complete_dir = _write_record(db_root, "Ascomycota", "Eurotiales", "222_b_MAT_MAT1-2")
    (complete_dir / "proteins.faa").write_text(">already|gene_index=0|name=G1|role=core_MAT\nMSEQ\n")

    missing = find_records_missing_proteins_faa(db_root)

    assert missing == [("Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")]


def test_find_records_missing_proteins_faa_excludes_candidates(tmp_path):
    db_root = tmp_path / "db"
    _write_record(db_root, "candidates", "Ascomycota", "999_c_MAT_MAT1-1")

    assert find_records_missing_proteins_faa(db_root) == []
```

- [ ] **Step 2: Run to verify it fails**

```bash
pixi run pytest tests/db/test_backfill_gff.py -v
```

Expected: FAIL with `ImportError: cannot import name 'find_records_missing_proteins_faa'`.

- [ ] **Step 3: Write the failing test for the batch backfill's failure isolation**

Append:

```python
from unittest.mock import MagicMock


def test_backfill_missing_proteins_faa_isolates_one_record_failure(tmp_path, monkeypatch):
    db_root = tmp_path / "db"
    _write_record(db_root, "Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    _write_record(db_root, "Ascomycota", "Onygenales", "222_b_MAT_MAT1-2")

    calls = []

    def fake_build(db_root, phylum, order_or_family, record_id, ncbi, uniprot):
        calls.append(record_id)
        if record_id == "111_a_MAT_MAT1-1":
            raise ConnectionError("simulated NCBI outage")

    monkeypatch.setattr("MATPredict.db.cli.build_gff_for_record", fake_build)

    from MATPredict.db.cli import backfill_missing_proteins_faa

    succeeded, failed = backfill_missing_proteins_faa(db_root, ncbi=MagicMock(), uniprot=MagicMock())

    assert calls == ["111_a_MAT_MAT1-1", "222_b_MAT_MAT1-2"]  # second record still attempted
    assert succeeded == [("Ascomycota", "Onygenales", "222_b_MAT_MAT1-2")]
    assert len(failed) == 1
    assert failed[0][0] == ("Ascomycota", "Eurotiales", "111_a_MAT_MAT1-1")
    assert "simulated NCBI outage" in failed[0][1]
```

- [ ] **Step 4: Run to verify it fails**

```bash
pixi run pytest tests/db/test_backfill_gff.py -v
```

Expected: FAIL (`backfill_missing_proteins_faa` doesn't exist yet).

- [ ] **Step 5: Implement the extraction and the new functions in `src/MATPredict/db/cli.py`**

Replace `_cmd_build_gff`'s body with a call to a new extracted helper, and add the two new functions near it:

```python
def find_records_missing_proteins_faa(db_root: Path) -> list[tuple[str, str, str]]:
    """(phylum, order_or_family, record_id) for every accepted (non-candidate)
    record whose proteins.faa does not exist on disk yet."""
    missing = []
    for meta_path in sorted(db_root.glob("*/*/*/metadata.yaml")):
        parts = meta_path.relative_to(db_root).parts
        if parts[0] == "candidates":
            continue
        record_dir = meta_path.parent
        if not (record_dir / "proteins.faa").exists():
            missing.append((parts[0], parts[1], parts[2]))
    return missing


def build_gff_for_record(
    db_root: Path, phylum: str, order_or_family: str, record_id: str,
    ncbi: NcbiClient, uniprot: UniprotClient,
) -> None:
    """Fetch every present gene's curated protein sequence and write
    locus.gff3/locus.gbk/proteins.faa for one accepted record. The single-
    record CLI command and the batch backfill command both call this so
    there is exactly one place this logic lives."""
    record_dir = db_root / phylum / order_or_family / record_id
    record = yaml.safe_load((record_dir / "metadata.yaml").read_text())

    sequences: dict[int, str] = {}
    for gene in record.get("genes", []):
        if not gene.get("present", True) or not gene.get("protein_accession"):
            continue
        client, bare_accession = _client_for(gene["protein_accession"], ncbi, uniprot)
        sequences[gene["gene_index"]] = client.fetch_protein_sequence(bare_accession)

    gff_export.write_gff3(record, out_path=record_dir / "locus.gff3")
    gff_export.write_genbank(record, sequences, out_path=record_dir / "locus.gbk")
    gff_export.write_proteins_fasta(record, sequences, out_path=record_dir / "proteins.faa")


def backfill_missing_proteins_faa(
    db_root: Path, ncbi: NcbiClient, uniprot: UniprotClient,
) -> tuple[list[tuple[str, str, str]], list[tuple[tuple[str, str, str], str]]]:
    """Run build_gff_for_record for every record find_records_missing_proteins_faa
    reports, isolating each record's failure so one live-fetch error (a
    suppressed accession, a transient NCBI outage) never aborts the rest of
    the batch -- the same discipline this project's genome-acquisition and
    batch-runner fixes already established."""
    succeeded: list[tuple[str, str, str]] = []
    failed: list[tuple[tuple[str, str, str], str]] = []
    for phylum, order_or_family, record_id in find_records_missing_proteins_faa(db_root):
        identifier = (phylum, order_or_family, record_id)
        try:
            build_gff_for_record(db_root, phylum, order_or_family, record_id, ncbi, uniprot)
        except Exception as exc:  # noqa: BLE001 -- recorded, never swallowed silently
            failed.append((identifier, str(exc)))
            continue
        succeeded.append(identifier)
    return succeeded, failed


def _cmd_build_gff(args: argparse.Namespace) -> int:
    config = _config(args)
    ncbi, uniprot = _make_clients(config)
    build_gff_for_record(config.db_root, args.phylum, args.order_or_family, args.record_id, ncbi, uniprot)
    record_dir = config.db_root / args.phylum / args.order_or_family / args.record_id
    print(f"wrote {record_dir / 'locus.gff3'}, {record_dir / 'locus.gbk'}, {record_dir / 'proteins.faa'}")
    return 0


def _cmd_backfill_gff(args: argparse.Namespace) -> int:
    config = _config(args)
    ncbi, uniprot = _make_clients(config)
    succeeded, failed = backfill_missing_proteins_faa(config.db_root, ncbi, uniprot)
    for phylum, order_or_family, record_id in succeeded:
        print(f"backfilled {phylum}/{order_or_family}/{record_id}")
    for (phylum, order_or_family, record_id), message in failed:
        print(f"FAILED {phylum}/{order_or_family}/{record_id}: {message}")
    print(f"{len(succeeded)} succeeded, {len(failed)} failed")
    return 0
```

Add the `backfill-gff` subparser next to `build_gff`'s existing registration:

```python
    backfill_gff = action.add_parser("backfill-gff")
    backfill_gff.set_defaults(func=_cmd_backfill_gff)
```

- [ ] **Step 6: Run tests to verify they pass**

```bash
pixi run pytest tests/db/test_backfill_gff.py -v
```

Expected: PASS, all 4 tests.

- [ ] **Step 7: Run the full suite**

```bash
pixi run pytest -v
```

Expected: all green, no regressions from the `_cmd_build_gff` extraction.

- [ ] **Step 8: Run the backfill for real against the live `db/` tree (a data step, not a code change)**

```bash
pixi run matpredict curate-db backfill-gff
```

If the failure count is more than a couple of the 31 records, stop and investigate before committing anything -- a high failure rate likely means a real bug in the extraction (e.g. a suppressed/withdrawn accession genuinely needs `validate.py`'s own handling, not a silent skip) rather than something to paper over.

- [ ] **Step 9: Spot-check a few real generated files**

```bash
head -4 db/Ascomycota/Onygenales/5501_h538-4_MAT_MAT1-1/proteins.faa
head -4 db/Ascomycota/Onygenales/199306_rmscc1040_MAT_MAT1-1/proteins.faa
```

Confirm real headers (`>{record_id}|gene_index=N|name=...|role=...`) and real amino-acid sequences, not empty files.

- [ ] **Step 10: Commit the code and the newly-generated data together**

```bash
git add src/MATPredict/db/cli.py tests/db/test_backfill_gff.py \
  db/Ascomycota/*/*/locus.gff3 db/Ascomycota/*/*/locus.gbk db/Ascomycota/*/*/proteins.faa
git commit -m "feat: backfill proteins.faa/locus.gff3/locus.gbk for 31 accepted records missing them"
```

---

### Task 3: A configurable, default-no-op Stage 1 → Stage 2 evidence gate, with calibration diagnostics

**Files:**
- Modify: `src/MATPredict/detect/pipeline.py`
- Test: `tests/detect/test_pipeline_evidence_floor.py`

**Interfaces:**
- Produces: `EvidenceFloor` (frozen dataclass: `min_hits: int = 1`, `min_identity: float | None = None`, `require_core_role: bool = False` -- these defaults exactly reproduce today's `_families_with_a_foothold` behavior, so shipping this changes nothing until a caller opts in), `_families_meeting_evidence_floor(cluster: GeneCluster, families: list[Family], floor: EvidenceFloor) -> list[Family]` (replaces `_families_with_a_foothold`, which becomes `_families_meeting_evidence_floor(cluster, families, EvidenceFloor())`), a new `run_pipeline(..., evidence_floor: EvidenceFloor = EvidenceFloor(), evidence_diagnostics_path: Path | None = None)` parameter pair.

**Why this task does NOT ship a tuned threshold:** the natural first instinct -- require a minimum tblastn identity before letting a family into the expensive polish loop -- was checked against this rollout's own real data before writing this plan, and it does not safely separate signal from noise yet. In the one genome that completed (5501, *C. immitis*), the family's own single real hit (`APN2`, `flanking_conserved`, matched against the wrong-order reference record) scored **34.9%** identity, while several of the highest-confidence *spurious* cross-phylum hits scored 50-73% (`Basidiomycota:bLocus bE` at 72.7%, `Ascomycota:PM matPi` at 64.3%, `Basidiomycota:HD HD1` at 63.6%). A naive identity floor set high enough to reject that noise would also have rejected the genome's one genuine signal. Gene-count is no better here: 55 of 67 detected results (including the real `Ascomycota:MAT` result) are single-gene hits, so "require >=2 genes" would have discarded the real signal too, on this example. **This task ships the mechanism, tested and off by default, plus a diagnostics mode to collect the real calibration data Task 1's fix makes possible** (once routing is fixed, most genomes will only ever see their 1-2 truly relevant families, so the diagnostics collected from THOSE runs -- not this noisy exhaustive-fallback run -- are the real calibration set). Setting non-default `EvidenceFloor` values is explicitly out of scope for this task; it is the natural fast-follow once that data exists.

- [ ] **Step 1: Write the failing test for `EvidenceFloor`'s default reproducing current behavior**

```python
# tests/detect/test_pipeline_evidence_floor.py
from __future__ import annotations

from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.pipeline import EvidenceFloor, _families_meeting_evidence_floor
from MATPredict.detect.search import SearchHit


def _family(phylum, name, genes=None):
    return Family(
        key=FamilyKey(phylum, name), vocabulary_type="enum",
        idiomorph_values=["a", "alpha"], idiomorph_pattern=None,
        genes=genes or [{"name": "G1", "role": "core_MAT"}],
        taxonomic_scope=[1],
    )


def _hit(family_key, gene_name="G1", role="core_MAT", identity=30.0):
    return SearchHit(
        family_key=family_key, gene_name=gene_name, role=role, contig="c1",
        start=1, end=100, strand="+", identity=identity,
        reference_record_id="rec1", method="tblastn_genome",
    )


def test_default_evidence_floor_admits_any_single_hit_regardless_of_identity():
    real = _family("Ascomycota", "MAT")
    spurious = _family("Basidiomycota", "bLocus")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(real.key, identity=34.9), _hit(spurious.key, identity=72.7),
    ])

    result = _families_meeting_evidence_floor(cluster, [real, spurious], EvidenceFloor())

    assert {f.key for f in result} == {real.key, spurious.key}


def test_evidence_floor_with_no_hits_admits_nothing():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[])

    assert _families_meeting_evidence_floor(cluster, [family], EvidenceFloor()) == []
```

- [ ] **Step 2: Run to verify it fails**

```bash
pixi run pytest tests/detect/test_pipeline_evidence_floor.py -v
```

Expected: FAIL with `ImportError: cannot import name 'EvidenceFloor'`.

- [ ] **Step 3: Write the failing tests for the non-default floor options**

Append:

```python
def test_min_hits_floor_rejects_a_single_hit_family():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[_hit(family.key)])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_hits=2))

    assert result == []


def test_min_hits_floor_admits_a_family_with_enough_distinct_genes():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1"), _hit(family.key, gene_name="G2"),
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_hits=2))

    assert [f.key for f in result] == [family.key]


def test_min_identity_floor_rejects_a_family_whose_best_hit_is_below_it():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[_hit(family.key, identity=34.9)])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(min_identity=50.0))

    assert result == []


def test_require_core_role_rejects_a_family_with_only_flanking_hits():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, role="flanking_conserved", identity=90.0),
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(require_core_role=True))

    assert result == []


def test_require_core_role_admits_a_family_with_at_least_one_core_hit():
    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, role="flanking_conserved", identity=90.0),
        _hit(family.key, gene_name="G1", role="core_MAT", identity=40.0),
    ])

    result = _families_meeting_evidence_floor(cluster, [family], EvidenceFloor(require_core_role=True))

    assert [f.key for f in result] == [family.key]
```

- [ ] **Step 4: Run to verify they fail, then implement**

```bash
pixi run pytest tests/detect/test_pipeline_evidence_floor.py -v
```

In `pipeline.py`, add the dataclass and replace `_families_with_a_foothold`:

```python
@dataclass(frozen=True)
class EvidenceFloor:
    """Minimum evidence a family must clear in a cluster before it is admitted
    to the (expensive, per-gene, two-subprocess-per-tool) Stage 2 polish loop.

    Defaults exactly reproduce the pre-existing `_families_with_a_foothold`
    behavior (any single hit of any role, any identity) -- this gate changes
    nothing until a caller sets a stricter floor. See run_pipeline's own
    `evidence_floor`/`evidence_diagnostics_path` parameters and this task's
    plan notes for why no non-default value is set here yet: this rollout's
    own real data shows a naive identity or hit-count cutoff is not yet safely
    separable from real signal, so tuning these is left to a calibration pass
    once db/Ascomycota/order.yml's taxonomic_scope fix (a separate task) lets
    most genomes run with correctly-narrowed routing instead of the
    exhaustive fallback that produced this rollout's noisy example.
    """

    min_hits: int = 1
    min_identity: float | None = None
    require_core_role: bool = False


def _families_meeting_evidence_floor(
    cluster: GeneCluster, families: list[Family], floor: EvidenceFloor
) -> list[Family]:
    """Families whose OWN hits in this cluster clear `floor` -- generalizes
    the old `_families_with_a_foothold` (which is exactly
    `_families_meeting_evidence_floor(cluster, families, EvidenceFloor())`)."""
    admitted = []
    for family in families:
        own_hits = [h for h in cluster.hits if h.family_key == family.key]
        if len(own_hits) < floor.min_hits:
            continue
        if floor.require_core_role:
            # SearchHit.role already carries "core_MAT | flanking_conserved |
            # flanking_variable" directly (search.py) -- filter on the HIT's
            # own role, not on whether the gene NAME happens to be one the
            # family defines as core_MAT, which would not actually test
            # anything (a hit's gene_name is only ever one the family
            # declares in the first place).
            own_hits = [h for h in own_hits if h.role == "core_MAT"]
            if not own_hits:
                continue
        if floor.min_identity is not None and max(h.identity for h in own_hits) < floor.min_identity:
            continue
        admitted.append(family)
    return admitted
```

Update the polish loop's call site (the `for family in _families_with_a_foothold(cluster, families):` line found earlier in this file) to `for family in _families_meeting_evidence_floor(cluster, families, evidence_floor):`, and add `evidence_floor: EvidenceFloor = EvidenceFloor()` to `run_pipeline`'s parameter list, threading it through to that call site.

- [ ] **Step 5: Run tests to verify they pass**

```bash
pixi run pytest tests/detect/test_pipeline_evidence_floor.py -v
```

Expected: PASS, all 7 tests.

- [ ] **Step 6: Write the failing test for the diagnostics dump**

```python
def test_evidence_diagnostics_written_when_path_given(tmp_path):
    from MATPredict.detect.pipeline import _write_evidence_diagnostics

    family = _family("Ascomycota", "MAT")
    cluster = GeneCluster(contig="c1", start=1, end=100, hits=[
        _hit(family.key, gene_name="G1", role="core_MAT", identity=34.9),
    ])
    out_path = tmp_path / "diagnostics.jsonl"

    _write_evidence_diagnostics(out_path, cluster, family, admitted=True)
    _write_evidence_diagnostics(out_path, cluster, family, admitted=False)

    import json
    lines = out_path.read_text().splitlines()
    assert len(lines) == 2
    row = json.loads(lines[0])
    assert row["family"] == "Ascomycota:MAT"
    assert row["contig"] == "c1"
    assert row["gene_count"] == 1
    assert row["best_identity"] == 34.9
    assert row["admitted"] is True
```

- [ ] **Step 7: Run to verify it fails, then implement**

```python
def _write_evidence_diagnostics(
    out_path: Path, cluster: GeneCluster, family: Family, admitted: bool
) -> None:
    """Append one JSON line describing this (cluster, family) admission
    decision -- the real calibration dataset `EvidenceFloor`'s docstring
    refers to. Never raises on a write failure; diagnostics are best-effort
    and must never abort a real detection run."""
    own_hits = [h for h in cluster.hits if h.family_key == family.key]
    row = {
        "family": f"{family.key.phylum}:{family.key.locus_name}",
        "contig": cluster.contig, "cluster_start": cluster.start, "cluster_end": cluster.end,
        "gene_count": len({h.gene_name for h in own_hits}),
        "roles": sorted({h.role for h in own_hits}),
        "best_identity": max((h.identity for h in own_hits), default=None),
        "admitted": admitted,
    }
    try:
        with out_path.open("a") as f:
            f.write(json.dumps(row) + "\n")
    except OSError:
        pass
```

Add `import json` to `pipeline.py`'s imports if not already present. Wire it into the polish loop: when `run_pipeline`'s new `evidence_diagnostics_path: Path | None = None` parameter is set, call `_write_evidence_diagnostics` once per (cluster, family) for every family `_families_with_a_foothold`'s OLD unconditional-admit set would have considered (i.e. every family with >=1 own hit in the cluster, not just the ones `_families_meeting_evidence_floor` actually admits) -- so the diagnostics capture both admitted and would-have-been-rejected cases for later threshold analysis, tagged with the real `evidence_floor` decision (`admitted=family in _families_meeting_evidence_floor(...)`).

- [ ] **Step 8: Run tests to verify they pass**

```bash
pixi run pytest tests/detect/test_pipeline_evidence_floor.py -v
```

Expected: PASS, all 8 tests.

- [ ] **Step 9: Run the full suite to confirm no regressions**

```bash
pixi run pytest -v
```

Expected: all green -- `run_pipeline`'s default `EvidenceFloor()` and `evidence_diagnostics_path=None` must not change any existing test's outcome, since both new parameters default to today's exact behavior.

- [ ] **Step 10: Commit**

```bash
git add src/MATPredict/detect/pipeline.py tests/detect/test_pipeline_evidence_floor.py
git commit -m "feat: add a configurable, default-no-op Stage1->Stage2 evidence floor with calibration diagnostics"
```

---

## Self-review notes (controller, at plan-writing time)

- **Spec coverage:** all three of the rollout findings report's triage items (stale `taxonomic_scope`, `proteins.faa` gap, no Stage1->Stage2 gate) map 1:1 to Tasks 1-3.
- **No placeholders:** every code block is real and complete against this session's actual, freshly-read source (`family_registry.py`, `reference_fasta.py`, `pipeline.py` lines ~192-960, `db/cli.py`, `db/gff_export.py`) -- no invented function names or signatures. The one deliberately-undetermined value (Task 3's real threshold numbers) is not a placeholder in the "TBD" sense the writing-plans skill forbids: it is a documented design decision (ship the mechanism, not a guessed number) backed by real data from this session's own rollout run (genome 5501's real per-hit identity table), not an unwritten detail.
- **Type/signature consistency:** `ScopeAuditResult`/`audit_scope` (Task 1) are self-contained, consumed only by Task 1's own CLI wiring. `build_gff_for_record`/`backfill_missing_proteins_faa` (Task 2) reuse `_client_for`/`_make_clients`/`NcbiClient`/`UniprotClient` exactly as `_cmd_build_gff` already does -- no new client type introduced. `EvidenceFloor`/`_families_meeting_evidence_floor` (Task 3) is a strict generalization of the existing `_families_with_a_foothold`, verified by a test asserting the default floor reproduces its exact behavior.
- **Task independence:** Tasks 1-3 touch disjoint files (`db/Ascomycota/order.yml` + new `scope_audit.py`; `db/cli.py` + generated `db/Ascomycota/*/*/proteins.faa` etc.; `pipeline.py`) and can be executed and reviewed in any order, though Task 1 should land first in practice since it's what makes Task 3's eventual real calibration data meaningful.
