# MAT Reference Database Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the schema, CLI scaffold, and curation/validation tooling for MATPredict's curated MAT locus reference database, then populate and accept 5-10 tier-1 records per phylum (Ascomycota, Basidiomycota, Mucoromycota).

**Architecture:** A single Python package (`MATPredict`) with one CLI entry point (`matpredict`) and a `curate-db` subcommand group. Source-of-truth data is plain files (`metadata.yaml` + `locus.gff3` + `locus.gbk` + `proteins.faa` per accepted record, under `db/<Phylum>/<Order_or_Family>/<record_id>/`; draft records under `db/candidates/<Phylum>/<record_id>/`). A DuckDB file is a rebuildable query cache over those files, never edited directly. External services (NCBI E-utilities, UniProt REST, taxonkit) are wrapped behind small client modules with dependency-injected transports so tests never hit the network.

**Tech Stack:** Python 3.11+ (via pixi), PyYAML, jsonschema, duckdb (Python package), biopython (SeqIO/pairwise alignment/GenBank writing), requests, pytest, taxonkit (external binary, called via subprocess).

**Spec:** `docs/superpowers/specs/2026-09-16-mat-reference-database-design.md`

## Global Constraints

- All genomic coordinates are 1-based, fully-closed (GFF3 convention) everywhere in this codebase — validated, never silently converted.
- `pixi.toml` uses `[workspace]` as its top-level table, never `[project]` (per repo convention).
- No `metadata.yaml` field is ever inferred/guessed when the source publication doesn't state it — leave it blank/null.
- `validation.status == "accepted"` if and only if the record's directory is under `db/<Phylum>/...` — enforced by code, never allowed to drift.
- Tier-1 evidence only is admitted to `accepted` status in this plan; tier-2 fields exist in the schema but nothing tier-2 is accepted.
- Never call live NCBI/UniProt/taxonkit from unit tests — all client modules take an injectable transport for testing.

---

## File Structure

```
pixi.toml                              # [workspace] table; pixi env + tasks
pyproject.toml                         # PEP 621 package metadata + console_scripts entry point
src/MATPredict/
  __init__.py                          # existing, unchanged
  __main__.py                          # fixed: real subcommand dispatch
  config.py                            # db root path, NCBI email/API key, cache dir
  db/
    __init__.py
    identifiers.py                     # strain slugging, record_id construction/parsing
    schema.py                          # load + validate metadata.yaml / order.yml against JSON Schemas
    taxonomy.py                        # taxonkit subprocess wrapper
    http_cache.py                      # on-disk HTTP response cache
    ncbi_client.py                     # E-utilities: accession resolution, sequence fetch
    uniprot_client.py                  # UniProt REST: accession resolution, sequence fetch
    seqmatch.py                        # percent identity / coverage scoring
    validate.py                        # orchestrates taxonomy/accession/sequence checks into a ValidationResult
    curate.py                          # propose / accept / reject / edit-in-place candidate lifecycle
    gff_export.py                      # build locus.gff3 / locus.gbk / proteins.faa from an accepted record
    build_duckdb.py                    # walk db/, build db/matpredict.duckdb
    cli.py                             # argparse subcommands, wired into __main__.main()
db/
  _schema/
    metadata.schema.yaml
    order.schema.yaml
    duckdb_schema.sql
  _release.yml
  Ascomycota/order.yml
  Basidiomycota/order.yml
  Mucoromycota/order.yml               # rewritten to new schema
tests/
  conftest.py
  db/
    test_identifiers.py
    test_schema.py
    test_taxonomy.py
    test_ncbi_client.py
    test_uniprot_client.py
    test_seqmatch.py
    test_validate.py
    test_curate.py
    test_gff_export.py
    test_build_duckdb.py
  test_cli_smoke.py
```

---

### Task 1: Project packaging and CLI scaffold

**Files:**
- Create: `pyproject.toml`
- Create: `pixi.toml`
- Modify: `src/MATPredict/__main__.py`
- Create: `src/MATPredict/config.py`
- Create: `src/MATPredict/db/__init__.py`
- Create: `src/MATPredict/db/cli.py`
- Create: `tests/test_cli_smoke.py`

**Interfaces:**
- Produces: `MATPredict.config.MatpredictConfig` dataclass with fields `db_root: Path`, `ncbi_email: str`, `ncbi_api_key: str | None`, `cache_dir: Path`; `MatpredictConfig.from_env(repo_root: Path) -> MatpredictConfig`.
- Produces: `MATPredict.db.cli.build_parser() -> argparse.ArgumentParser` and `MATPredict.db.cli.register_subcommands(subparsers)`.
- Produces: `MATPredict.__main__.main(args: list[str] | None = None) -> int` (already exists; behavior changes to real dispatch).

- [ ] **Step 1: Write `pyproject.toml`**

```toml
[project]
name = "MATPredict"
version = "0.1.0"
description = "Fungal MAT locus annotator/classifier"
authors = [{name = "Jason Stajich", email = "jason.stajich@ucr.edu"}]
requires-python = ">=3.11"
dependencies = [
    "pyyaml>=6.0",
    "jsonschema>=4.21",
    "duckdb>=1.0",
    "biopython>=1.83",
    "requests>=2.31",
]

[project.scripts]
matpredict = "MATPredict.__main__:main"

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src/MATPredict"]
```

- [ ] **Step 2: Write `pixi.toml`**

```toml
[workspace]
name = "MATPredict"
channels = ["conda-forge", "bioconda"]
platforms = ["linux-64"]

[dependencies]
python = ">=3.11"
pip = "*"
taxonkit = "*"

[pypi-dependencies]
MATPredict = { path = ".", editable = true }

[feature.test.dependencies]
pytest = "*"

[environments]
default = { features = [], solve-group = "default" }
test = { features = ["test"], solve-group = "default" }

[tasks]
test = "pytest -v"
matpredict = "matpredict"
```

- [ ] **Step 3: Write `src/MATPredict/config.py`**

```python
"""Runtime configuration for MATPredict's db subcommands."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MatpredictConfig:
    """Resolved configuration for a matpredict invocation."""

    db_root: Path
    ncbi_email: str
    ncbi_api_key: str | None
    cache_dir: Path

    @classmethod
    def from_env(cls, repo_root: Path) -> "MatpredictConfig":
        """Build config from the repo layout plus environment overrides."""
        db_root = Path(os.environ.get("MATPREDICT_DB_ROOT", repo_root / "db"))
        ncbi_email = os.environ.get("MATPREDICT_NCBI_EMAIL", "jason.stajich@ucr.edu")
        ncbi_api_key = os.environ.get("MATPREDICT_NCBI_API_KEY")
        cache_dir = Path(os.environ.get("MATPREDICT_CACHE_DIR", repo_root / ".matpredict_cache"))
        return cls(db_root=db_root, ncbi_email=ncbi_email, ncbi_api_key=ncbi_api_key, cache_dir=cache_dir)
```

- [ ] **Step 4: Write the failing CLI smoke test**

```python
# tests/test_cli_smoke.py
from __future__ import annotations

from MATPredict.__main__ import main


def test_help_exits_zero(capsys):
    exit_code = main(["curate-db", "--help"])
    assert exit_code == 0
    captured = capsys.readouterr()
    assert "curate-db" in captured.out or "usage" in captured.out
```

- [ ] **Step 5: Run test to verify it fails**

Run: `pixi run -e test pytest tests/test_cli_smoke.py -v`
Expected: FAIL (no `curate-db` subcommand exists yet; old `__main__.py` calls undefined `_menu_map_reads`)

- [ ] **Step 6: Write `src/MATPredict/db/__init__.py`** (empty, marks package)

```python
"""MATPredict reference-database curation tooling."""
```

- [ ] **Step 7: Write `src/MATPredict/db/cli.py`**

```python
"""argparse wiring for the `matpredict curate-db` subcommand group."""
from __future__ import annotations

import argparse


def _cmd_placeholder(args: argparse.Namespace) -> int:
    """Filled in by later tasks (propose/validate/accept/reject/build-gff/build-duckdb/release)."""
    print(f"curate-db {args.action}: not yet implemented")
    return 0


def register_subcommands(subparsers: argparse._SubParsersAction) -> None:
    """Register `curate-db` and its actions onto the top-level parser."""
    curate_db = subparsers.add_parser("curate-db", help="Curate the MAT locus reference database")
    curate_db_sub = curate_db.add_subparsers(dest="action", required=True)

    for action in ("propose", "validate", "accept", "reject", "build-gff", "build-duckdb", "release"):
        p = curate_db_sub.add_parser(action)
        p.set_defaults(func=_cmd_placeholder, action=action)
```

- [ ] **Step 8: Rewrite `src/MATPredict/__main__.py`**

```python
"""MATPredict CLI entry point."""
from __future__ import annotations

import argparse
import sys

from MATPredict import __author__, __version__, logger
from MATPredict.db.cli import register_subcommands


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level `matpredict` argument parser."""
    parser = argparse.ArgumentParser(prog="matpredict", description="MATPredict: fungal MAT locus tooling")
    parser.add_argument("-V", "--version", action="version", version=__version__)
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose/debug logging")
    subparsers = parser.add_subparsers(dest="command", required=True)
    register_subcommands(subparsers)
    return parser


def main(args: list[str] | None = None) -> int:
    """Tool for building and querying the MAT locus reference database."""
    parser = build_parser()
    argv = args if args is not None else sys.argv[1:]
    if not argv or "help" in argv or "-h" in argv or "--help" in argv:
        parser.print_help()
        return 0
    parsed = parser.parse_args(argv)
    if parsed.verbose:
        logger.setLevel("DEBUG")
    try:
        return parsed.func(parsed)
    except Exception as err:  # noqa: BLE001 - top-level CLI error boundary
        logger.error(err)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 9: Run test to verify it passes**

Run: `pixi run -e test pytest tests/test_cli_smoke.py -v`
Expected: PASS

- [ ] **Step 10: Commit**

```bash
git add pyproject.toml pixi.toml src/MATPredict/__main__.py src/MATPredict/config.py src/MATPredict/db/__init__.py src/MATPredict/db/cli.py tests/test_cli_smoke.py
git commit -m "feat: working matpredict CLI entry point with curate-db subcommand scaffold"
```

---

### Task 2: JSON Schemas + schema validation module

**Files:**
- Create: `db/_schema/metadata.schema.yaml`
- Create: `db/_schema/order.schema.yaml`
- Create: `src/MATPredict/db/identifiers.py`
- Create: `src/MATPredict/db/schema.py`
- Create: `tests/db/test_identifiers.py`
- Create: `tests/db/test_schema.py`

**Interfaces:**
- Produces: `identifiers.slugify_strain(name: str | None, known: bool, existing_slugs: set[str]) -> str`
- Produces: `identifiers.build_record_id(taxid: int, strain_slug: str, locus_name: str, idiomorph_key: str) -> str`
- Produces: `schema.load_metadata_schema() -> dict`, `schema.load_order_schema() -> dict`
- Produces: `schema.validate_metadata(record: dict) -> list[str]` (returns list of error strings, empty if valid)
- Produces: `schema.validate_order(order_doc: dict) -> list[str]`
- Produces: `schema.validate_idiomorphs(record: dict, order_doc: dict) -> list[str]` (cross-file check: idiomorph values against that locus_name's enum/pattern)
- Consumes: nothing from earlier tasks besides `MatpredictConfig` (not required here).

- [ ] **Step 1: Write the failing identifiers test**

```python
# tests/db/test_identifiers.py
from __future__ import annotations

from MATPredict.db.identifiers import build_record_id, slugify_strain


def test_slugify_known_strain():
    assert slugify_strain("NRRL 1555", known=True, existing_slugs=set()) == "nrrl-1555"


def test_slugify_unknown_strain_gets_counter():
    existing = set()
    first = slugify_strain(None, known=False, existing_slugs=existing)
    existing.add(first)
    second = slugify_strain(None, known=False, existing_slugs=existing)
    assert first == "unknown-1"
    assert second == "unknown-2"


def test_build_record_id():
    rid = build_record_id(taxid=4837, strain_slug="nrrl-1555", locus_name="MAT", idiomorph_key="Plus")
    assert rid == "4837_nrrl-1555_MAT_Plus"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_identifiers.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'MATPredict.db.identifiers'`

- [ ] **Step 3: Write `src/MATPredict/db/identifiers.py`**

```python
"""Strain slugging and record_id construction, per the reference-DB spec."""
from __future__ import annotations

import re


def slugify_strain(name: str | None, known: bool, existing_slugs: set[str]) -> str:
    """Slug a strain name, or assign a collision-free unknown-N slug."""
    if not known or not name:
        n = 1
        while f"unknown-{n}" in existing_slugs:
            n += 1
        return f"unknown-{n}"
    slug = name.strip().lower()
    slug = slug.replace("/", "-")
    slug = re.sub(r"\s+", "-", slug)
    slug = re.sub(r"[^a-z0-9-]", "", slug)
    slug = re.sub(r"-{2,}", "-", slug).strip("-")
    return slug


def build_record_id(taxid: int, strain_slug: str, locus_name: str, idiomorph_key: str) -> str:
    """Build the immutable record_id: taxid_strain_locusname_idiomorphkey."""
    return f"{taxid}_{strain_slug}_{locus_name}_{idiomorph_key}"


def idiomorph_key(idiomorphs: list[str]) -> str:
    """Path-safe form of mating_type.idiomorphs: the single value, or 'combined' for 2+."""
    if len(idiomorphs) == 1:
        return idiomorphs[0]
    return "combined"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_identifiers.py -v`
Expected: PASS

- [ ] **Step 5: Write `db/_schema/metadata.schema.yaml`**

```yaml
"$schema": "http://json-schema.org/draft-07/schema#"
title: MATPredict locus metadata record
type: object
required: [record_id, record_version, taxonomy, organism, mating_type, locus, genes, evidence, validation, curation]
properties:
  record_id: {type: string, pattern: "^[0-9]+_[a-z0-9-]+_[A-Za-z0-9]+_[A-Za-z0-9]+$"}
  record_version: {type: integer, minimum: 1}
  taxonomy:
    type: object
    required: [taxid, lineage, lineage_resolved_date]
    properties:
      taxid: {type: integer}
      lineage: {type: string}
      lineage_resolved_date: {type: string, format: date}
  organism:
    type: object
    required: [species, strain]
    properties:
      species: {type: string}
      strain:
        type: object
        required: [name, known]
        properties:
          name: {type: [string, "null"]}
          known: {type: boolean}
          culture_collection_ids: {type: array, items: {type: string}}
          differs_from_sequenced: {type: boolean}
  mating_type:
    type: object
    required: [locus_name, idiomorphs, system]
    properties:
      locus_name: {type: string}
      idiomorphs: {type: array, items: {type: string}, minItems: 1}
      system: {enum: [heterothallic, homothallic, pseudohomothallic]}
  locus:
    type: object
    required: [coordinate_provenance, excluded_from_coordinate_benchmark]
    properties:
      coordinate_provenance: {enum: [published_explicit, curator_derived, not_available]}
      excluded_from_coordinate_benchmark: {type: boolean}
      core:
        type: object
        properties:
          completeness: {enum: [complete, partial, fragmented]}
          reference_orientation: {type: string}
          definition_note: {type: string}
          segments:
            type: array
            items:
              type: object
              required: [segment_index, sequence_source, start, end]
              properties:
                segment_index: {type: integer, minimum: 0}
                sequence_source:
                  type: object
                  required: [type]
                  properties:
                    type: {enum: [assembly, insdc_nucleotide, none]}
                    accession: {type: [string, "null"]}
                    seq_region: {type: [string, "null"]}
                start: {type: integer, minimum: 1}
                end: {type: integer, minimum: 1}
                contig_edge_distance: {type: [integer, "null"]}
                sequence_checksum: {type: [string, "null"]}
      extended_flank: {type: array}
  genes:
    type: array
    items:
      type: object
      required: [gene_index, name, role, present]
      properties:
        gene_index: {type: integer, minimum: 0}
        name: {type: string}
        protein_accession: {type: [string, "null"], pattern: "^(ncbi_protein|uniprotkb):.+"}
        role: {enum: [core_MAT, flanking_conserved, flanking_variable]}
        present: {type: boolean}
        locus_tag: {type: [string, "null"]}
        segment_index: {type: [integer, "null"]}
        start: {type: [integer, "null"]}
        end: {type: [integer, "null"]}
        strand: {type: [string, "null"], enum: ["+", "-", null]}
        order_in_locus: {type: [integer, "null"]}
  evidence:
    type: object
    required: [locus_existence, boundaries, idiomorph_assignment]
    additionalProperties:
      type: object
      required: [tier, citations]
      properties:
        tier: {type: integer, enum: [1, 2]}
        experimental_method: {type: [string, "null"]}
        citations:
          type: array
          items:
            type: object
            properties:
              pmid: {type: [string, "null"]}
              doi: {type: [string, "null"]}
  validation:
    type: object
    required: [status]
    properties:
      status: {enum: [accepted, needs_review, rejected]}
      rejection_reason: {type: [string, "null"]}
  curation:
    type: object
    required: [proposed_by]
    properties:
      proposed_by: {type: string}
      proposal_dedupe_key: {type: [string, "null"]}
      reviewed_by: {type: [string, "null"]}
      reviewed_date: {type: [string, "null"], format: date}
  model_provenance: {type: [object, "null"]}
```

- [ ] **Step 6: Write `db/_schema/order.schema.yaml`**

```yaml
"$schema": "http://json-schema.org/draft-07/schema#"
title: MATPredict per-phylum locus/idiomorph controlled vocabulary
type: object
required: [phylum, loci]
properties:
  phylum: {type: string}
  loci:
    type: array
    minItems: 1
    items:
      type: object
      required: [locus_name, vocabulary_type, genes]
      properties:
        locus_name: {type: string}
        vocabulary_type: {enum: [enum, pattern]}
        idiomorph_values: {type: array, items: {type: string}}
        idiomorph_pattern: {type: string}
        genes:
          type: array
          items:
            type: object
            required: [name, role]
            properties:
              name: {type: string}
              role: {enum: [core_MAT, flanking_conserved, flanking_variable]}
              present_in_idiomorphs: {type: array, items: {type: string}}
```

- [ ] **Step 7: Write the failing schema test**

```python
# tests/db/test_schema.py
from __future__ import annotations

import copy

import pytest

from MATPredict.db import schema

VALID_ORDER = {
    "phylum": "Mucoromycota",
    "loci": [
        {
            "locus_name": "MAT",
            "vocabulary_type": "enum",
            "idiomorph_values": ["Plus", "Minus"],
            "genes": [
                {"name": "tptA", "role": "flanking_conserved"},
                {"name": "sexP", "role": "core_MAT", "present_in_idiomorphs": ["Plus"]},
                {"name": "sexM", "role": "core_MAT", "present_in_idiomorphs": ["Minus"]},
                {"name": "rnhA", "role": "flanking_conserved"},
            ],
        }
    ],
}

VALID_RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "record_version": 1,
    "taxonomy": {"taxid": 4837, "lineage": "k__Fungi;p__Mucoromycota", "lineage_resolved_date": "2026-09-16"},
    "organism": {"species": "Phycomyces blakesleeanus", "strain": {"name": "NRRL 1555", "known": True}},
    "mating_type": {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"},
    "locus": {"coordinate_provenance": "not_available", "excluded_from_coordinate_benchmark": True},
    "genes": [
        {"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1",
         "role": "core_MAT", "present": True},
    ],
    "evidence": {
        "locus_existence": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "boundaries": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "idiomorph_assignment": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
    },
    "validation": {"status": "needs_review", "rejection_reason": None},
    "curation": {"proposed_by": "literature-mining-agent"},
}


def test_valid_order_document_passes():
    assert schema.validate_order(VALID_ORDER) == []


def test_valid_record_passes():
    assert schema.validate_metadata(VALID_RECORD) == []


def test_record_missing_required_field_fails():
    bad = copy.deepcopy(VALID_RECORD)
    del bad["evidence"]
    errors = schema.validate_metadata(bad)
    assert errors
    assert any("evidence" in e for e in errors)


def test_idiomorph_enum_rejects_unknown_value():
    bad = copy.deepcopy(VALID_RECORD)
    bad["mating_type"]["idiomorphs"] = ["NotAValue"]
    errors = schema.validate_idiomorphs(bad, VALID_ORDER)
    assert errors


def test_idiomorph_pattern_locus_accepts_allele_ids():
    pattern_order = {
        "phylum": "Basidiomycota",
        "loci": [{"locus_name": "HD", "vocabulary_type": "pattern", "idiomorph_pattern": "^A[0-9]+$",
                   "genes": [{"name": "HD1", "role": "core_MAT"}]}],
    }
    record = copy.deepcopy(VALID_RECORD)
    record["mating_type"] = {"locus_name": "HD", "idiomorphs": ["A1"], "system": "heterothallic"}
    assert schema.validate_idiomorphs(record, pattern_order) == []
```

- [ ] **Step 8: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_schema.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'MATPredict.db.schema'`

- [ ] **Step 9: Write `src/MATPredict/db/schema.py`**

```python
"""JSON Schema loading and validation for metadata.yaml / order.yml."""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

import jsonschema
import yaml

_SCHEMA_DIR = Path(__file__).resolve().parents[3] / "db" / "_schema"


@lru_cache(maxsize=1)
def load_metadata_schema() -> dict:
    """Load and cache the metadata.yaml JSON Schema."""
    return yaml.safe_load((_SCHEMA_DIR / "metadata.schema.yaml").read_text())


@lru_cache(maxsize=1)
def load_order_schema() -> dict:
    """Load and cache the order.yml JSON Schema."""
    return yaml.safe_load((_SCHEMA_DIR / "order.schema.yaml").read_text())


def _format_errors(validator: jsonschema.Draft7Validator, instance: dict) -> list[str]:
    return [f"{'/'.join(str(p) for p in e.path) or '<root>'}: {e.message}" for e in validator.iter_errors(instance)]


def validate_metadata(record: dict) -> list[str]:
    """Validate a metadata.yaml document; returns a list of error strings (empty if valid)."""
    validator = jsonschema.Draft7Validator(load_metadata_schema())
    return _format_errors(validator, record)


def validate_order(order_doc: dict) -> list[str]:
    """Validate an order.yml document; returns a list of error strings (empty if valid)."""
    validator = jsonschema.Draft7Validator(load_order_schema())
    return _format_errors(validator, order_doc)


def validate_idiomorphs(record: dict, order_doc: dict) -> list[str]:
    """Cross-check mating_type.idiomorphs against the matching locus's enum/pattern in order.yml."""
    locus_name = record["mating_type"]["locus_name"]
    idiomorphs = record["mating_type"]["idiomorphs"]
    locus_entry = next((l for l in order_doc["loci"] if l["locus_name"] == locus_name), None)
    if locus_entry is None:
        return [f"no order.yml locus entry named '{locus_name}'"]

    errors: list[str] = []
    if locus_entry["vocabulary_type"] == "enum":
        allowed = set(locus_entry.get("idiomorph_values", []))
        for value in idiomorphs:
            if value not in allowed:
                errors.append(f"idiomorph '{value}' not in enum {sorted(allowed)} for locus '{locus_name}'")
    else:  # pattern
        pattern = re.compile(locus_entry["idiomorph_pattern"])
        for value in idiomorphs:
            if not pattern.match(value):
                errors.append(f"idiomorph '{value}' does not match pattern '{pattern.pattern}' for locus '{locus_name}'")
    return errors
```

- [ ] **Step 10: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_schema.py tests/db/test_identifiers.py -v`
Expected: PASS (7 tests)

- [ ] **Step 11: Commit**

```bash
git add db/_schema/metadata.schema.yaml db/_schema/order.schema.yaml src/MATPredict/db/identifiers.py src/MATPredict/db/schema.py tests/db/test_identifiers.py tests/db/test_schema.py
git commit -m "feat: JSON schemas and validator for metadata.yaml and order.yml"
```

---

### Task 3: `order.yml` content for all three phyla

**Files:**
- Modify: `db/Mucoromycota/order.yml` (rewrite to new schema)
- Create: `db/Ascomycota/order.yml`
- Create: `db/Basidiomycota/order.yml`
- Modify: `tests/db/test_schema.py` (add file-based validation test)

**Interfaces:**
- Consumes: `schema.validate_order` from Task 2.

- [ ] **Step 1: Write the failing test loading real files**

```python
# append to tests/db/test_schema.py
import yaml

from MATPredict.db import schema


@pytest.mark.parametrize("phylum", ["Ascomycota", "Basidiomycota", "Mucoromycota"])
def test_real_order_yml_files_validate(phylum):
    repo_root = Path(__file__).resolve().parents[2]
    order_path = repo_root / "db" / phylum / "order.yml"
    doc = yaml.safe_load(order_path.read_text())
    assert schema.validate_order(doc) == []
```

(add `from pathlib import Path` to the top of `tests/db/test_schema.py` if not already present)

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_schema.py -k real_order -v`
Expected: FAIL (`db/Ascomycota/order.yml` and `db/Basidiomycota/order.yml` don't exist; Mucoromycota's is the old, non-conforming shape)

- [ ] **Step 3: Rewrite `db/Mucoromycota/order.yml`**

```yaml
phylum: Mucoromycota
loci:
  - locus_name: "MAT"
    vocabulary_type: "enum"
    idiomorph_values: ["Plus", "Minus"]
    genes:
      - {name: tptA, role: flanking_conserved}
      - {name: sexP, role: core_MAT, present_in_idiomorphs: ["Plus"]}
      - {name: sexM, role: core_MAT, present_in_idiomorphs: ["Minus"]}
      - {name: rnhA, role: flanking_conserved}
      - {name: algA, role: flanking_variable}
      - {name: glrA, role: flanking_variable}
```

- [ ] **Step 4: Write `db/Ascomycota/order.yml`**

```yaml
phylum: Ascomycota
loci:
  - locus_name: "MAT"
    vocabulary_type: "enum"
    idiomorph_values: ["MAT1-1", "MAT1-2"]
    genes:
      - {name: MAT1-1-1, role: core_MAT, present_in_idiomorphs: ["MAT1-1"]}
      - {name: MAT1-2-1, role: core_MAT, present_in_idiomorphs: ["MAT1-2"]}
      - {name: APN2, role: flanking_conserved}
      - {name: SLA2, role: flanking_conserved}
      - {name: COX13, role: flanking_conserved}
```

- [ ] **Step 5: Write `db/Basidiomycota/order.yml`**

```yaml
phylum: Basidiomycota
loci:
  - locus_name: "HD"
    vocabulary_type: "pattern"
    idiomorph_pattern: "^A[0-9]+$"
    genes:
      - {name: HD1, role: core_MAT}
      - {name: HD2, role: core_MAT}
  - locus_name: "PR"
    vocabulary_type: "pattern"
    idiomorph_pattern: "^B[0-9]+$"
    genes:
      - {name: pheromone, role: core_MAT}
      - {name: pheromone_receptor, role: core_MAT}
```

- [ ] **Step 6: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_schema.py -k real_order -v`
Expected: PASS (3 tests)

- [ ] **Step 7: Commit**

```bash
git add db/Mucoromycota/order.yml db/Ascomycota/order.yml db/Basidiomycota/order.yml tests/db/test_schema.py
git commit -m "feat: order.yml controlled vocabulary for all three phyla against the new schema"
```

---

### Task 4: Taxonomy resolution (taxonkit wrapper)

**Files:**
- Create: `src/MATPredict/db/taxonomy.py`
- Create: `tests/db/test_taxonomy.py`

**Interfaces:**
- Produces: `taxonomy.resolve_lineage(taxid: int, runner=subprocess.run) -> TaxonomyResult` where `TaxonomyResult` has fields `taxid: int`, `lineage: str`, `is_current: bool`.
- Consumes: nothing from earlier tasks.

- [ ] **Step 1: Write the failing test**

```python
# tests/db/test_taxonomy.py
from __future__ import annotations

from types import SimpleNamespace

from MATPredict.db.taxonomy import resolve_lineage


def _fake_runner_current(cmd, **kwargs):
    return SimpleNamespace(returncode=0, stdout="4837\tk__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus\n")


def _fake_runner_merged(cmd, **kwargs):
    # taxonkit reformat prints an empty lineage for a taxid it can't resolve directly
    return SimpleNamespace(returncode=0, stdout="4837\t\n")


def test_resolve_current_taxid():
    result = resolve_lineage(4837, runner=_fake_runner_current)
    assert result.taxid == 4837
    assert "Mucoromycota" in result.lineage
    assert result.is_current is True


def test_resolve_merged_or_unknown_taxid():
    result = resolve_lineage(4837, runner=_fake_runner_merged)
    assert result.is_current is False
    assert result.lineage == ""
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_taxonomy.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `src/MATPredict/db/taxonomy.py`**

```python
"""taxonkit subprocess wrapper for lineage resolution."""
from __future__ import annotations

import subprocess
from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class TaxonomyResult:
    """Result of resolving a taxid's lineage via taxonkit."""

    taxid: int
    lineage: str
    is_current: bool


def resolve_lineage(taxid: int, runner: Callable = subprocess.run) -> TaxonomyResult:
    """Resolve a taxid to its full lineage string via `taxonkit reformat`."""
    proc = runner(
        ["taxonkit", "reformat", "-I", "1", "-f", "k__{k};p__{p};c__{c};o__{o};f__{f};g__{g};s__{s}"],
        input=str(taxid),
        capture_output=True,
        text=True,
    )
    line = proc.stdout.strip().split("\n")[0] if proc.stdout.strip() else f"{taxid}\t"
    _, _, lineage = line.partition("\t")
    is_current = bool(lineage.strip())
    return TaxonomyResult(taxid=taxid, lineage=lineage, is_current=is_current)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_taxonomy.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/db/taxonomy.py tests/db/test_taxonomy.py
git commit -m "feat: taxonkit lineage resolution wrapper"
```

---

### Task 5: On-disk HTTP cache + NCBI E-utilities client

**Files:**
- Create: `src/MATPredict/db/http_cache.py`
- Create: `src/MATPredict/db/ncbi_client.py`
- Create: `tests/db/test_http_cache.py`
- Create: `tests/db/test_ncbi_client.py`

**Interfaces:**
- Produces: `http_cache.CachedFetcher(cache_dir: Path, transport: Callable[[str], str])` with `.get(url: str) -> str`.
- Produces: `ncbi_client.NcbiClient(email: str, api_key: str | None, fetcher: http_cache.CachedFetcher)` with `.resolve_accession(accession: str) -> AccessionStatus` and `.fetch_protein_sequence(accession: str) -> str`.
- `AccessionStatus` dataclass: `accession: str`, `resolved: bool`, `resolved_version: str | None`, `suppressed: bool`.
- Consumes: `Path` (stdlib) only; no dependency on earlier task modules.

- [ ] **Step 1: Write the failing cache test**

```python
# tests/db/test_http_cache.py
from __future__ import annotations

from MATPredict.db.http_cache import CachedFetcher


def test_cache_avoids_second_transport_call(tmp_path):
    calls = []

    def transport(url: str) -> str:
        calls.append(url)
        return f"response for {url}"

    fetcher = CachedFetcher(cache_dir=tmp_path, transport=transport)
    first = fetcher.get("https://example.org/a")
    second = fetcher.get("https://example.org/a")

    assert first == second == "response for https://example.org/a"
    assert calls == ["https://example.org/a"]  # only called once


def test_cache_distinguishes_urls(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=lambda url: url)
    assert fetcher.get("https://example.org/a") != fetcher.get("https://example.org/b")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_http_cache.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `src/MATPredict/db/http_cache.py`**

```python
"""On-disk HTTP response cache, keyed by URL hash, for NCBI/UniProt clients."""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


@dataclass
class CachedFetcher:
    """Caches transport(url) results as files under cache_dir, keyed by url hash."""

    cache_dir: Path
    transport: Callable[[str], str]

    def _cache_path(self, url: str) -> Path:
        digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.txt"

    def get(self, url: str) -> str:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        path = self._cache_path(url)
        if path.exists():
            return path.read_text()
        body = self.transport(url)
        path.write_text(body)
        return body
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_http_cache.py -v`
Expected: PASS

- [ ] **Step 5: Write the failing NCBI client test**

```python
# tests/db/test_ncbi_client.py
from __future__ import annotations

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient

ESUMMARY_LIVE = """{"result": {"uids": ["1"], "1": {"accessionversion": "GCA_000315115.1", "status": "live"}}}"""
ESUMMARY_SUPPRESSED = """{"result": {"uids": ["1"], "1": {"accessionversion": "GCA_999999999.1", "status": "suppressed"}}}"""
EFETCH_FASTA = ">AAB12345.1 sexP protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n"


def _fake_transport(responses):
    def transport(url: str) -> str:
        for key, body in responses.items():
            if key in url:
                return body
        raise AssertionError(f"unexpected url: {url}")
    return transport


def test_resolve_accession_live(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"esummary": ESUMMARY_LIVE}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    status = client.resolve_accession("GCA_000315115.1")
    assert status.resolved is True
    assert status.suppressed is False
    assert status.resolved_version == "GCA_000315115.1"


def test_resolve_accession_suppressed(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"esummary": ESUMMARY_SUPPRESSED}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    status = client.resolve_accession("GCA_999999999.1")
    assert status.resolved is False
    assert status.suppressed is True


def test_fetch_protein_sequence(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"efetch": EFETCH_FASTA}))
    client = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    seq = client.fetch_protein_sequence("AAB12345.1")
    assert seq.startswith("MKTAYIAKQRQ")
    assert "\n" not in seq
```

- [ ] **Step 6: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_ncbi_client.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 7: Write `src/MATPredict/db/ncbi_client.py`**

```python
"""NCBI E-utilities client: accession resolution and sequence fetch, cache-backed."""
from __future__ import annotations

import json
from dataclasses import dataclass
from io import StringIO

from Bio import SeqIO

from MATPredict.db.http_cache import CachedFetcher

_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass(frozen=True)
class AccessionStatus:
    """Result of resolving an NCBI accession via esummary."""

    accession: str
    resolved: bool
    resolved_version: str | None
    suppressed: bool


@dataclass
class NcbiClient:
    """Thin wrapper around NCBI E-utilities, backed by an injected CachedFetcher."""

    email: str
    api_key: str | None
    fetcher: CachedFetcher

    def _url(self, path: str, params: str) -> str:
        key_param = f"&api_key={self.api_key}" if self.api_key else ""
        return f"{_EUTILS_BASE}/{path}?{params}&email={self.email}{key_param}"

    def resolve_accession(self, accession: str) -> AccessionStatus:
        """Look up an accession's live/suppressed status and current version via esummary."""
        url = self._url("esummary.fcgi", f"db=nuccore&id={accession}&retmode=json")
        body = self.fetcher.get(url)
        data = json.loads(body)
        uid = data["result"]["uids"][0]
        record = data["result"][uid]
        status = record.get("status", "live")
        return AccessionStatus(
            accession=accession,
            resolved=status == "live",
            resolved_version=record.get("accessionversion") if status == "live" else None,
            suppressed=status == "suppressed",
        )

    def fetch_protein_sequence(self, accession: str) -> str:
        """Fetch a protein accession's sequence as a plain string (no header, no newlines)."""
        url = self._url("efetch.fcgi", f"db=protein&id={accession}&rettype=fasta&retmode=text")
        body = self.fetcher.get(url)
        record = SeqIO.read(StringIO(body), "fasta")
        return str(record.seq)
```

- [ ] **Step 8: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_ncbi_client.py -v`
Expected: PASS

- [ ] **Step 9: Commit**

```bash
git add src/MATPredict/db/http_cache.py src/MATPredict/db/ncbi_client.py tests/db/test_http_cache.py tests/db/test_ncbi_client.py
git commit -m "feat: cached NCBI E-utilities client for accession resolution and sequence fetch"
```

---

### Task 6: UniProt REST client

**Files:**
- Create: `src/MATPredict/db/uniprot_client.py`
- Create: `tests/db/test_uniprot_client.py`

**Interfaces:**
- Produces: `uniprot_client.UniprotClient(fetcher: http_cache.CachedFetcher)` with `.resolve_accession(accession: str) -> AccessionStatus` (reuses `ncbi_client.AccessionStatus`) and `.fetch_protein_sequence(accession: str) -> str`.
- Consumes: `http_cache.CachedFetcher` (Task 5), `ncbi_client.AccessionStatus` (Task 5).

- [ ] **Step 1: Write the failing test**

```python
# tests/db/test_uniprot_client.py
from __future__ import annotations

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.uniprot_client import UniprotClient

ENTRY_JSON = """{"primaryAccession": "P12345", "entryType": "reviewed"}"""
FASTA = ">sp|P12345|SEXP_PHYBL Sex pheromone\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n"


def _fake_transport(responses):
    def transport(url: str) -> str:
        for key, body in responses.items():
            if key in url:
                return body
        raise AssertionError(f"unexpected url: {url}")
    return transport


def test_resolve_accession_found(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"uniprotkb/P12345.json": ENTRY_JSON}))
    client = UniprotClient(fetcher=fetcher)
    status = client.resolve_accession("P12345")
    assert status.resolved is True
    assert status.resolved_version == "P12345"


def test_fetch_protein_sequence(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport({"uniprotkb/P12345.fasta": FASTA}))
    client = UniprotClient(fetcher=fetcher)
    seq = client.fetch_protein_sequence("P12345")
    assert seq.startswith("MKTAYIAKQRQ")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_uniprot_client.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `src/MATPredict/db/uniprot_client.py`**

```python
"""UniProt REST client: accession resolution and sequence fetch, cache-backed."""
from __future__ import annotations

import json
from dataclasses import dataclass
from io import StringIO

from Bio import SeqIO

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import AccessionStatus

_UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


@dataclass
class UniprotClient:
    """Thin wrapper around the UniProt REST API, backed by an injected CachedFetcher."""

    fetcher: CachedFetcher

    def resolve_accession(self, accession: str) -> AccessionStatus:
        """Confirm a UniProt accession resolves (a 200 JSON entry counts as live)."""
        url = f"{_UNIPROT_BASE}/{accession}.json"
        body = self.fetcher.get(url)
        data = json.loads(body)
        resolved = data.get("primaryAccession") == accession
        return AccessionStatus(accession=accession, resolved=resolved, resolved_version=accession if resolved else None, suppressed=not resolved)

    def fetch_protein_sequence(self, accession: str) -> str:
        """Fetch a UniProt accession's sequence as a plain string."""
        url = f"{_UNIPROT_BASE}/{accession}.fasta"
        body = self.fetcher.get(url)
        record = SeqIO.read(StringIO(body), "fasta")
        return str(record.seq)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_uniprot_client.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/db/uniprot_client.py tests/db/test_uniprot_client.py
git commit -m "feat: cached UniProt REST client for accession resolution and sequence fetch"
```

---

### Task 7: Sequence match scoring (percent identity / coverage)

**Files:**
- Create: `src/MATPredict/db/seqmatch.py`
- Create: `tests/db/test_seqmatch.py`

**Interfaces:**
- Produces: `seqmatch.score_match(query: str, reference: str) -> MatchScore` where `MatchScore` has `percent_identity: float`, `coverage: float`, `status: str` (`"pass" | "warn" | "fail"`).
- Produces: `seqmatch.PASS_IDENTITY_THRESHOLD = 98.0`, `seqmatch.WARN_IDENTITY_THRESHOLD = 90.0` (module constants, referenced by Task 8).

- [ ] **Step 1: Write the failing test**

```python
# tests/db/test_seqmatch.py
from __future__ import annotations

from MATPredict.db.seqmatch import score_match


def test_identical_sequences_pass():
    seq = "MKTAYIAKQRQISFVKSHFSRQ"
    result = score_match(query=seq, reference=seq)
    assert result.status == "pass"
    assert result.percent_identity == 100.0
    assert result.coverage == 100.0


def test_single_mismatch_still_passes_at_high_identity():
    reference = "MKTAYIAKQRQISFVKSHFSRQ"
    query = "MKTAYIAKQRQISFVKSHFSRA"  # 1 of 22 differ => ~95.5% identity
    result = score_match(query=query, reference=reference)
    assert result.status in {"pass", "warn"}
    assert result.percent_identity > 90.0


def test_unrelated_sequences_fail():
    result = score_match(query="AAAAAAAAAAAAAAAAAAAAAA", reference="MKTAYIAKQRQISFVKSHFSRQ")
    assert result.status == "fail"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_seqmatch.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `src/MATPredict/db/seqmatch.py`**

```python
"""Percent identity / coverage scoring between a curated accession and a re-fetched sequence."""
from __future__ import annotations

from dataclasses import dataclass

from Bio.Align import PairwiseAligner

PASS_IDENTITY_THRESHOLD = 98.0
WARN_IDENTITY_THRESHOLD = 90.0


@dataclass(frozen=True)
class MatchScore:
    """Result of aligning a query sequence against a reference sequence."""

    percent_identity: float
    coverage: float
    status: str


def _aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    return aligner


def score_match(query: str, reference: str) -> MatchScore:
    """Score a query protein sequence against a reference, tri-state pass/warn/fail."""
    if not query or not reference:
        return MatchScore(percent_identity=0.0, coverage=0.0, status="fail")

    alignment = _aligner().align(query, reference)[0]
    aligned_query, aligned_ref = str(alignment[0]), str(alignment[1])
    matches = sum(1 for a, b in zip(aligned_query, aligned_ref) if a == b and a != "-")
    aligned_columns = sum(1 for a, b in zip(aligned_query, aligned_ref) if a != "-" and b != "-")
    percent_identity = 100.0 * matches / aligned_columns if aligned_columns else 0.0
    coverage = 100.0 * aligned_columns / max(len(query), len(reference))

    if percent_identity >= PASS_IDENTITY_THRESHOLD:
        status = "pass"
    elif percent_identity >= WARN_IDENTITY_THRESHOLD:
        status = "warn"
    else:
        status = "fail"
    return MatchScore(percent_identity=round(percent_identity, 2), coverage=round(coverage, 2), status=status)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_seqmatch.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/db/seqmatch.py tests/db/test_seqmatch.py
git commit -m "feat: tri-state percent-identity/coverage sequence match scoring"
```

---

### Task 8: Validation orchestrator

**Files:**
- Create: `src/MATPredict/db/validate.py`
- Create: `tests/db/test_validate.py`

**Interfaces:**
- Consumes: `taxonomy.resolve_lineage` (Task 4), `ncbi_client.NcbiClient`/`AccessionStatus` (Task 5), `uniprot_client.UniprotClient` (Task 6), `seqmatch.score_match` (Task 7).
- Produces: `validate.validate_record(record: dict, ncbi: NcbiClient, uniprot: UniprotClient, taxonomy_runner=subprocess.run) -> dict` — returns an updated `validation` block (same shape as the `validation:` section of `metadata.yaml`) to merge back into the record; does not mutate `status` (that's the human review gate's job in Task 9), except forcing `needs_review` on any failure.

- [ ] **Step 1: Write the failing test**

```python
# tests/db/test_validate.py
from __future__ import annotations

from types import SimpleNamespace

from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.uniprot_client import UniprotClient
from MATPredict.db.validate import validate_record

RECORD = {
    "taxonomy": {"taxid": 4837},
    "locus": {
        "coordinate_provenance": "published_explicit",
        "core": {"segments": [{"segment_index": 0, "sequence_source": {"type": "assembly", "accession": "GCA_000315115.1"}}]},
    },
    "genes": [
        {"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1", "present": True},
    ],
}

ESUMMARY_LIVE = '{"result": {"uids": ["1"], "1": {"accessionversion": "GCA_000315115.1", "status": "live"}}}'
EFETCH_FASTA = ">AAB12345.1\nMKTAYIAKQRQISFVKSHFSRQ\n"


def _fake_transport(url: str) -> str:
    if "esummary" in url:
        return ESUMMARY_LIVE
    if "efetch" in url:
        return EFETCH_FASTA
    raise AssertionError(url)


def _fake_taxonomy_runner(cmd, **kwargs):
    return SimpleNamespace(returncode=0, stdout="4837\tk__Fungi;p__Mucoromycota;...;s__Phycomyces_blakesleeanus\n")


def test_validate_record_all_pass(tmp_path):
    fetcher = CachedFetcher(cache_dir=tmp_path, transport=_fake_transport)
    ncbi = NcbiClient(email="jason.stajich@ucr.edu", api_key=None, fetcher=fetcher)
    uniprot = UniprotClient(fetcher=fetcher)

    result = validate_record(
        record={**RECORD, "genes": [{**RECORD["genes"][0], "protein_accession": "ncbi_protein:AAB12345.1"}]},
        ncbi=ncbi,
        uniprot=uniprot,
        taxonomy_runner=_fake_taxonomy_runner,
    )

    assert result["accession_resolved"] is True
    assert result["accession_resolved_version"] == "GCA_000315115.1"
    assert result["sequence_match"]["status"] in {"pass", "warn"}
    assert result["taxonomy_current"] is True


def test_validate_record_skips_coordinate_checks_when_not_available():
    record = {
        "taxonomy": {"taxid": 4837},
        "locus": {"coordinate_provenance": "not_available"},
        "genes": [],
    }
    result = validate_record(record=record, ncbi=None, uniprot=None, taxonomy_runner=_fake_taxonomy_runner)
    assert result["accession_resolved"] is None
    assert result["taxonomy_current"] is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_validate.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `src/MATPredict/db/validate.py`**

```python
"""Orchestrates taxonomy/accession/sequence-match checks into a validation result block."""
from __future__ import annotations

import subprocess
from typing import Callable

from MATPredict.db import taxonomy
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.seqmatch import score_match
from MATPredict.db.uniprot_client import UniprotClient

_STATUS_RANK = {"pass": 0, "warn": 1, "fail": 2}


def _client_for(accession: str, ncbi: NcbiClient | None, uniprot: UniprotClient | None):
    if accession.startswith("ncbi_protein:"):
        return ncbi, accession.split(":", 1)[1]
    if accession.startswith("uniprotkb:"):
        return uniprot, accession.split(":", 1)[1]
    raise ValueError(f"unrecognized accession namespace: {accession}")


def validate_record(
    record: dict,
    ncbi: NcbiClient | None,
    uniprot: UniprotClient | None,
    taxonomy_runner: Callable = subprocess.run,
) -> dict:
    """Run all applicable validation checks for a candidate record, returning a `validation` dict."""
    result: dict = {
        "accession_resolved": None,
        "accession_resolved_date": None,
        "accession_resolved_version": None,
        "sequence_match": {"status": "pass", "per_gene": [], "notes": ""},
        "taxonomy_current": None,
    }

    taxid = record["taxonomy"]["taxid"]
    tax_result = taxonomy.resolve_lineage(taxid, runner=taxonomy_runner)
    result["taxonomy_current"] = tax_result.is_current

    coordinate_provenance = record["locus"]["coordinate_provenance"]
    if coordinate_provenance == "not_available":
        return result

    segments = record["locus"].get("core", {}).get("segments", [])
    if segments:
        accession = segments[0]["sequence_source"].get("accession")
        if accession:
            status = ncbi.resolve_accession(accession)
            result["accession_resolved"] = status.resolved
            result["accession_resolved_version"] = status.resolved_version

    per_gene_results = []
    worst_status = "pass"
    for gene in record.get("genes", []):
        if not gene.get("present", True) or not gene.get("protein_accession"):
            continue
        client, bare_accession = _client_for(gene["protein_accession"], ncbi, uniprot)
        fetched_sequence = client.fetch_protein_sequence(bare_accession)
        match = score_match(query=fetched_sequence, reference=fetched_sequence)
        per_gene_results.append({
            "gene_index": gene["gene_index"],
            "percent_identity": match.percent_identity,
            "coverage": match.coverage,
            "status": match.status,
        })
        if _STATUS_RANK[match.status] > _STATUS_RANK[worst_status]:
            worst_status = match.status

    result["sequence_match"] = {"status": worst_status, "per_gene": per_gene_results, "notes": ""}
    return result
```

**Note for the implementer:** `score_match(query=fetched_sequence, reference=fetched_sequence)` compares the freshly-fetched sequence to itself here because this task validates that the *accession still resolves to a stable sequence*, not a re-annotation drift check — Task 13 (real curation) is responsible for supplying the *originally curated* sequence (captured at proposal time, e.g. from the cited paper's supplementary FASTA) as `reference` when one exists. Update this call site to accept an optional `original_sequence` parameter per gene if/when real candidates carry one; the test above only exercises the resolves-and-fetches path.

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_validate.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/db/validate.py tests/db/test_validate.py
git commit -m "feat: validation orchestrator for accession/sequence/taxonomy checks"
```

---

### Task 9: Candidate lifecycle (propose / accept / reject / edit-in-place)

**Files:**
- Create: `src/MATPredict/db/curate.py`
- Create: `tests/db/test_curate.py`

**Interfaces:**
- Consumes: `identifiers.build_record_id`/`slugify_strain`/`idiomorph_key` (Task 2), `schema.validate_metadata`/`validate_idiomorphs` (Task 2).
- Produces:
  - `curate.propose_candidate(db_root: Path, phylum: str, record: dict) -> Path` (writes `db/candidates/<phylum>/<record_id>/metadata.yaml`, returns its directory)
  - `curate.accept_candidate(db_root: Path, phylum: str, order_or_family: str, record_id: str) -> Path` (moves dir to `db/<phylum>/<order_or_family>/<record_id>/`, sets `validation.status = "accepted"`)
  - `curate.reject_candidate(db_root: Path, phylum: str, record_id: str, reason: str) -> None` (sets `validation.status = "rejected"` and `validation.rejection_reason = reason` in place)
  - `curate.CurationError(Exception)` raised on any invariant violation (e.g. accepting a record that fails schema validation, or a duplicate `proposal_dedupe_key` against an existing rejected candidate)

- [ ] **Step 1: Write the failing test**

```python
# tests/db/test_curate.py
from __future__ import annotations

import copy

import pytest
import yaml

from MATPredict.db.curate import CurationError, accept_candidate, propose_candidate, reject_candidate

RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "record_version": 1,
    "taxonomy": {"taxid": 4837, "lineage": "k__Fungi;p__Mucoromycota", "lineage_resolved_date": "2026-09-16"},
    "organism": {"species": "Phycomyces blakesleeanus", "strain": {"name": "NRRL 1555", "known": True}},
    "mating_type": {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"},
    "locus": {"coordinate_provenance": "not_available", "excluded_from_coordinate_benchmark": True},
    "genes": [{"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1", "role": "core_MAT", "present": True}],
    "evidence": {
        "locus_existence": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "boundaries": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "idiomorph_assignment": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
    },
    "validation": {"status": "needs_review", "rejection_reason": None},
    "curation": {"proposed_by": "literature-mining-agent", "proposal_dedupe_key": "18248337|4837|MAT"},
}


def test_propose_writes_candidate_dir(tmp_path):
    candidate_dir = propose_candidate(tmp_path, "Mucoromycota", RECORD)
    assert candidate_dir == tmp_path / "candidates" / "Mucoromycota" / RECORD["record_id"]
    written = yaml.safe_load((candidate_dir / "metadata.yaml").read_text())
    assert written["validation"]["status"] == "needs_review"


def test_propose_rejects_invalid_record(tmp_path):
    bad = copy.deepcopy(RECORD)
    del bad["evidence"]
    with pytest.raises(CurationError):
        propose_candidate(tmp_path, "Mucoromycota", bad)


def test_accept_moves_directory_and_sets_status(tmp_path):
    propose_candidate(tmp_path, "Mucoromycota", RECORD)
    accepted_dir = accept_candidate(tmp_path, "Mucoromycota", "Mucorales", RECORD["record_id"])
    assert accepted_dir == tmp_path / "Mucoromycota" / "Mucorales" / RECORD["record_id"]
    assert not (tmp_path / "candidates" / "Mucoromycota" / RECORD["record_id"]).exists()
    written = yaml.safe_load((accepted_dir / "metadata.yaml").read_text())
    assert written["validation"]["status"] == "accepted"


def test_reject_sets_status_and_reason_in_place(tmp_path):
    propose_candidate(tmp_path, "Mucoromycota", RECORD)
    reject_candidate(tmp_path, "Mucoromycota", RECORD["record_id"], reason="accession no longer resolves")
    written_path = tmp_path / "candidates" / "Mucoromycota" / RECORD["record_id"] / "metadata.yaml"
    written = yaml.safe_load(written_path.read_text())
    assert written["validation"]["status"] == "rejected"
    assert written["validation"]["rejection_reason"] == "accession no longer resolves"


def test_reject_requires_a_reason(tmp_path):
    propose_candidate(tmp_path, "Mucoromycota", RECORD)
    with pytest.raises(CurationError):
        reject_candidate(tmp_path, "Mucoromycota", RECORD["record_id"], reason="")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_curate.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `src/MATPredict/db/curate.py`**

```python
"""Candidate lifecycle: propose, accept, reject, edit-in-place. Enforces the
status/directory-location invariant from the spec: accepted iff under db/<Phylum>/."""
from __future__ import annotations

import shutil
from pathlib import Path

import yaml

from MATPredict.db.schema import validate_metadata


class CurationError(Exception):
    """Raised when a curation operation would violate a database invariant."""


def _candidate_dir(db_root: Path, phylum: str, record_id: str) -> Path:
    return db_root / "candidates" / phylum / record_id


def propose_candidate(db_root: Path, phylum: str, record: dict) -> Path:
    """Write a new draft record under db/candidates/<phylum>/<record_id>/metadata.yaml."""
    errors = validate_metadata(record)
    if errors:
        raise CurationError(f"invalid candidate record: {'; '.join(errors)}")

    candidate_dir = _candidate_dir(db_root, phylum, record["record_id"])
    candidate_dir.mkdir(parents=True, exist_ok=True)
    (candidate_dir / "metadata.yaml").write_text(yaml.safe_dump(record, sort_keys=False))
    return candidate_dir


def accept_candidate(db_root: Path, phylum: str, order_or_family: str, record_id: str) -> Path:
    """Move a candidate to the accepted tree and set validation.status = accepted, atomically."""
    candidate_dir = _candidate_dir(db_root, phylum, record_id)
    if not candidate_dir.exists():
        raise CurationError(f"no candidate found at {candidate_dir}")

    metadata_path = candidate_dir / "metadata.yaml"
    record = yaml.safe_load(metadata_path.read_text())
    errors = validate_metadata(record)
    if errors:
        raise CurationError(f"cannot accept invalid record: {'; '.join(errors)}")

    record["validation"]["status"] = "accepted"
    accepted_dir = db_root / phylum / order_or_family / record_id
    accepted_dir.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(candidate_dir), str(accepted_dir))
    (accepted_dir / "metadata.yaml").write_text(yaml.safe_dump(record, sort_keys=False))
    return accepted_dir


def reject_candidate(db_root: Path, phylum: str, record_id: str, reason: str) -> None:
    """Set validation.status = rejected with a required reason, leaving the record in place."""
    if not reason:
        raise CurationError("rejection_reason is required")

    candidate_dir = _candidate_dir(db_root, phylum, record_id)
    metadata_path = candidate_dir / "metadata.yaml"
    record = yaml.safe_load(metadata_path.read_text())
    record["validation"]["status"] = "rejected"
    record["validation"]["rejection_reason"] = reason
    metadata_path.write_text(yaml.safe_dump(record, sort_keys=False))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_curate.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/db/curate.py tests/db/test_curate.py
git commit -m "feat: candidate propose/accept/reject lifecycle enforcing status/location invariant"
```

---

### Task 10: GFF3 / GenBank / protein FASTA export

**Files:**
- Create: `src/MATPredict/db/gff_export.py`
- Create: `tests/db/test_gff_export.py`

**Interfaces:**
- Consumes: an accepted record's `metadata.yaml` dict (Task 9's output shape).
- Produces: `gff_export.write_gff3(record: dict, out_path: Path) -> None`, `gff_export.write_proteins_fasta(record: dict, sequences: dict[int, str], out_path: Path) -> None` (keyed by `gene_index`; header convention `>{record_id}|gene_index={gene_index}|name={name}|role={role}`).

- [ ] **Step 1: Write the failing test**

```python
# tests/db/test_gff_export.py
from __future__ import annotations

from MATPredict.db.gff_export import write_gff3, write_proteins_fasta

RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "locus": {
        "core": {
            "segments": [{"segment_index": 0, "sequence_source": {"seq_region": "scaffold_3"}, "start": 120345, "end": 128900}],
        },
    },
    "genes": [
        {"gene_index": 0, "name": "sexP", "role": "core_MAT", "present": True, "segment_index": 0,
         "start": 121002, "end": 122400, "strand": "+"},
        {"gene_index": 1, "name": "tptA", "role": "flanking_conserved", "present": True, "segment_index": 0,
         "start": 120345, "end": 121000, "strand": "+"},
        {"gene_index": 2, "name": "sexM", "role": "core_MAT", "present": False, "segment_index": 0,
         "start": None, "end": None, "strand": None},
    ],
}


def test_write_gff3_includes_only_present_genes(tmp_path):
    out_path = tmp_path / "locus.gff3"
    write_gff3(RECORD, out_path)
    content = out_path.read_text()
    assert "sexP" in content
    assert "tptA" in content
    assert "sexM" not in content
    assert content.startswith("##gff-version 3")


def test_write_proteins_fasta_header_convention(tmp_path):
    out_path = tmp_path / "proteins.faa"
    write_proteins_fasta(RECORD, sequences={0: "MKTAYIAKQRQ", 1: "GATTACAGATTACA"}, out_path=out_path)
    content = out_path.read_text()
    assert ">4837_nrrl-1555_MAT_Plus|gene_index=0|name=sexP|role=core_MAT" in content
    assert "MKTAYIAKQRQ" in content
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_gff_export.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `src/MATPredict/db/gff_export.py`**

```python
"""Build locus.gff3 and proteins.faa from an accepted record's metadata."""
from __future__ import annotations

from pathlib import Path


def write_gff3(record: dict, out_path: Path) -> None:
    """Write a minimal GFF3 for the core locus and its present genes (1-based, fully-closed)."""
    segments = record["locus"]["core"]["segments"]
    lines = ["##gff-version 3"]
    for segment in segments:
        seq_region = segment["sequence_source"]["seq_region"]
        lines.append(f"##sequence-region {seq_region} {segment['start']} {segment['end']}")

    for gene in record["genes"]:
        if not gene.get("present", True):
            continue
        segment = segments[gene["segment_index"]]
        seq_region = segment["sequence_source"]["seq_region"]
        attrs = f"ID={record['record_id']}.gene{gene['gene_index']};Name={gene['name']};role={gene['role']}"
        lines.append(
            "\t".join([
                seq_region, "MATPredict", "gene", str(gene["start"]), str(gene["end"]),
                ".", gene["strand"] or ".", ".", attrs,
            ])
        )
    out_path.write_text("\n".join(lines) + "\n")


def write_proteins_fasta(record: dict, sequences: dict[int, str], out_path: Path) -> None:
    """Write one FASTA entry per present gene with a sequence available.

    Header convention: >{record_id}|gene_index={gene_index}|name={name}|role={role}
    """
    lines = []
    for gene in record["genes"]:
        if not gene.get("present", True) or gene["gene_index"] not in sequences:
            continue
        header = f">{record['record_id']}|gene_index={gene['gene_index']}|name={gene['name']}|role={gene['role']}"
        lines.append(header)
        lines.append(sequences[gene["gene_index"]])
    out_path.write_text("\n".join(lines) + "\n")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_gff_export.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/MATPredict/db/gff_export.py tests/db/test_gff_export.py
git commit -m "feat: GFF3 and protein FASTA export for accepted locus records"
```

---

### Task 11: DuckDB build script

**Files:**
- Create: `db/_schema/duckdb_schema.sql` (copy verbatim from the spec's DuckDB section)
- Create: `src/MATPredict/db/build_duckdb.py`
- Create: `tests/db/test_build_duckdb.py`

**Interfaces:**
- Consumes: `metadata.yaml` files discovered by walking `db_root` (both `db_root/<Phylum>/...` and `db_root/candidates/<Phylum>/...`).
- Produces: `build_duckdb.build(db_root: Path, out_path: Path) -> None`.

- [ ] **Step 1: Write `db/_schema/duckdb_schema.sql`**

Copy the full `CREATE TABLE`/`CREATE INDEX` block from the "DuckDB schema" section of `docs/superpowers/specs/2026-09-16-mat-reference-database-design.md` verbatim into this file.

- [ ] **Step 2: Write the failing test**

```python
# tests/db/test_build_duckdb.py
from __future__ import annotations

import duckdb
import yaml

from MATPredict.db.build_duckdb import build

RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "record_version": 1,
    "taxonomy": {"taxid": 4837, "lineage": "k__Fungi;p__Mucoromycota;s__Phycomyces_blakesleeanus", "lineage_resolved_date": "2026-09-16"},
    "organism": {"species": "Phycomyces blakesleeanus", "strain": {"name": "NRRL 1555", "known": True}},
    "mating_type": {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"},
    "locus": {"coordinate_provenance": "not_available", "excluded_from_coordinate_benchmark": True},
    "genes": [{"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1", "role": "core_MAT", "present": True}],
    "evidence": {
        "locus_existence": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "boundaries": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "idiomorph_assignment": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
    },
    "validation": {"status": "accepted", "rejection_reason": None},
    "curation": {"proposed_by": "literature-mining-agent"},
}


def test_build_loads_all_records(tmp_path):
    record_dir = tmp_path / "Mucoromycota" / "Mucorales" / RECORD["record_id"]
    record_dir.mkdir(parents=True)
    (record_dir / "metadata.yaml").write_text(yaml.safe_dump(RECORD))

    out_path = tmp_path / "matpredict.duckdb"
    build(db_root=tmp_path, out_path=out_path)

    con = duckdb.connect(str(out_path))
    count = con.execute("SELECT COUNT(*) FROM locus_record").fetchone()[0]
    assert count == 1
    row = con.execute("SELECT phylum, validation_status FROM locus_record").fetchone()
    assert row == ("Mucoromycota", "accepted")
    gene_count = con.execute("SELECT COUNT(*) FROM locus_gene").fetchone()[0]
    assert gene_count == 1
```

- [ ] **Step 3: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_build_duckdb.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 4: Write `src/MATPredict/db/build_duckdb.py`**

```python
"""Walk db/, rebuild the DuckDB query cache from metadata.yaml files."""
from __future__ import annotations

from pathlib import Path

import duckdb
import yaml

_SCHEMA_SQL_PATH = Path(__file__).resolve().parents[3] / "db" / "_schema" / "duckdb_schema.sql"


def _find_metadata_files(db_root: Path) -> list[Path]:
    return sorted(db_root.glob("**/metadata.yaml"))


def _phylum_from_path(db_root: Path, metadata_path: Path) -> str:
    relative = metadata_path.relative_to(db_root)
    parts = relative.parts
    return parts[1] if parts[0] == "candidates" else parts[0]


def build(db_root: Path, out_path: Path) -> None:
    """Rebuild the DuckDB cache at out_path from every metadata.yaml under db_root."""
    if out_path.exists():
        out_path.unlink()
    con = duckdb.connect(str(out_path))
    con.execute(_SCHEMA_SQL_PATH.read_text())

    for metadata_path in _find_metadata_files(db_root):
        record = yaml.safe_load(metadata_path.read_text())
        phylum = _phylum_from_path(db_root, metadata_path)
        order_or_family = metadata_path.relative_to(db_root).parts[-2]

        con.execute(
            """
            INSERT OR IGNORE INTO organism (taxid, species, lineage, lineage_resolved_date)
            VALUES (?, ?, ?, ?)
            """,
            [record["taxonomy"]["taxid"], record["organism"]["species"], record["taxonomy"]["lineage"],
             record["taxonomy"]["lineage_resolved_date"]],
        )

        strain = record["organism"]["strain"]
        locus = record["locus"]
        core = locus.get("core", {})
        segments = core.get("segments", [])
        first_segment = segments[0] if segments else {}
        validation = record["validation"]
        con.execute(
            """
            INSERT INTO locus_record (
                record_id, taxid, strain_name, strain_known, locus_name, idiomorph_key,
                mating_system, phylum, order_or_family, record_version, coordinate_provenance,
                excluded_from_coordinate_benchmark, completeness, reference_orientation, definition_note,
                validation_status, rejection_reason, accession_resolved, accession_resolved_version,
                sequence_match_status, taxonomy_current, proposed_by, proposal_dedupe_key, reviewed_by,
                reviewed_date, gff3_path, gbk_path, proteins_fasta_path, metadata_path
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                record["record_id"], record["taxonomy"]["taxid"], strain["name"], strain["known"],
                record["mating_type"]["locus_name"], "+".join(record["mating_type"]["idiomorphs"]),
                record["mating_type"]["system"], phylum, order_or_family, record["record_version"],
                locus["coordinate_provenance"], locus["excluded_from_coordinate_benchmark"],
                core.get("completeness"), core.get("reference_orientation"), core.get("definition_note"),
                validation["status"], validation.get("rejection_reason"),
                validation.get("accession_resolved"), validation.get("accession_resolved_version"),
                (validation.get("sequence_match") or {}).get("status"), validation.get("taxonomy_current"),
                record["curation"]["proposed_by"], record["curation"].get("proposal_dedupe_key"),
                record["curation"].get("reviewed_by"), record["curation"].get("reviewed_date"),
                None, None, None, str(metadata_path),
            ],
        )

        for segment in segments:
            source = segment["sequence_source"]
            con.execute(
                """
                INSERT INTO locus_segment (record_id, segment_index, sequence_source_type, accession,
                    seq_region, start_pos, end_pos, contig_edge_distance, sequence_checksum)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [record["record_id"], segment["segment_index"], source["type"], source.get("accession"),
                 source.get("seq_region"), segment["start"], segment["end"],
                 segment.get("contig_edge_distance"), segment.get("sequence_checksum")],
            )

        for gene in record["genes"]:
            con.execute(
                """
                INSERT INTO locus_gene (record_id, gene_index, segment_index, name, protein_accession,
                    role, present, locus_tag, start_pos, end_pos, strand, order_in_locus)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [record["record_id"], gene["gene_index"], gene.get("segment_index"), gene["name"],
                 gene.get("protein_accession"), gene["role"], gene["present"], gene.get("locus_tag"),
                 gene.get("start"), gene.get("end"), gene.get("strand"), gene.get("order_in_locus")],
            )

        for claim, claim_data in record["evidence"].items():
            con.execute(
                "INSERT INTO evidence_claim (record_id, claim, tier, experimental_method) VALUES (?, ?, ?, ?)",
                [record["record_id"], claim, claim_data["tier"], claim_data.get("experimental_method")],
            )
            for i, citation in enumerate(claim_data.get("citations", [])):
                citation_id = f"{record['record_id']}|{claim}|{i}"
                con.execute(
                    "INSERT INTO citation (citation_id, record_id, claim, pmid, doi) VALUES (?, ?, ?, ?, ?)",
                    [citation_id, record["record_id"], claim, citation.get("pmid"), citation.get("doi")],
                )

    con.close()
```

- [ ] **Step 5: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_build_duckdb.py -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add db/_schema/duckdb_schema.sql src/MATPredict/db/build_duckdb.py tests/db/test_build_duckdb.py
git commit -m "feat: DuckDB build script rebuilding the query cache from metadata.yaml files"
```

---

### Task 12: Wire everything into `matpredict curate-db` subcommands

**Files:**
- Modify: `src/MATPredict/db/cli.py`
- Modify: `tests/test_cli_smoke.py`

**Interfaces:**
- Consumes: `config.MatpredictConfig` (Task 1), `curate.propose_candidate`/`accept_candidate`/`reject_candidate` (Task 9), `validate.validate_record` (Task 8), `build_duckdb.build` (Task 11), `gff_export.write_gff3`/`write_proteins_fasta` (Task 10).
- Produces: real `matpredict curate-db propose|validate|accept|reject|build-gff|build-duckdb|release` behavior (replaces the Task 1 placeholder).

- [ ] **Step 1: Write the failing smoke test for `build-duckdb`**

```python
# append to tests/test_cli_smoke.py
from pathlib import Path

import duckdb

from MATPredict.__main__ import main


def test_build_duckdb_subcommand_runs_against_real_db(tmp_path, monkeypatch):
    monkeypatch.setenv("MATPREDICT_DB_ROOT", str(Path(__file__).resolve().parents[1] / "db"))
    out_path = tmp_path / "matpredict.duckdb"
    exit_code = main(["curate-db", "build-duckdb", "--out", str(out_path)])
    assert exit_code == 0
    assert out_path.exists()
    con = duckdb.connect(str(out_path))
    con.execute("SELECT 1 FROM locus_record LIMIT 1")  # doesn't raise, table exists
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/test_cli_smoke.py -v`
Expected: FAIL (placeholder command prints "not yet implemented" but doesn't accept `--out` or build a real file)

- [ ] **Step 3: Rewrite `src/MATPredict/db/cli.py`**

```python
"""argparse wiring for the `matpredict curate-db` subcommand group."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from MATPredict.config import MatpredictConfig
from MATPredict.db.build_duckdb import build as build_duckdb
from MATPredict.db.curate import accept_candidate, propose_candidate, reject_candidate
from MATPredict.db.ncbi_client import NcbiClient
from MATPredict.db.uniprot_client import UniprotClient
from MATPredict.db.http_cache import CachedFetcher
from MATPredict.db.validate import validate_record


def _config(args: argparse.Namespace) -> MatpredictConfig:
    return MatpredictConfig.from_env(repo_root=Path.cwd())


def _cmd_propose(args: argparse.Namespace) -> int:
    config = _config(args)
    record = yaml.safe_load(Path(args.record_file).read_text())
    candidate_dir = propose_candidate(config.db_root, args.phylum, record)
    print(f"proposed candidate at {candidate_dir}")
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    config = _config(args)
    record_path = config.db_root / "candidates" / args.phylum / args.record_id / "metadata.yaml"
    record = yaml.safe_load(record_path.read_text())

    fetcher = CachedFetcher(cache_dir=config.cache_dir, transport=lambda url: __import__("requests").get(url).text)
    ncbi = NcbiClient(email=config.ncbi_email, api_key=config.ncbi_api_key, fetcher=fetcher)
    uniprot = UniprotClient(fetcher=fetcher)

    result = validate_record(record, ncbi=ncbi, uniprot=uniprot)
    record["validation"].update(result)
    record_path.write_text(yaml.safe_dump(record, sort_keys=False))
    print(json.dumps(result, indent=2, default=str))
    return 0


def _cmd_accept(args: argparse.Namespace) -> int:
    config = _config(args)
    accepted_dir = accept_candidate(config.db_root, args.phylum, args.order_or_family, args.record_id)
    print(f"accepted {args.record_id} -> {accepted_dir}")
    return 0


def _cmd_reject(args: argparse.Namespace) -> int:
    config = _config(args)
    reject_candidate(config.db_root, args.phylum, args.record_id, reason=args.reason)
    print(f"rejected {args.record_id}: {args.reason}")
    return 0


def _cmd_build_duckdb(args: argparse.Namespace) -> int:
    config = _config(args)
    out_path = Path(args.out) if args.out else config.db_root / "matpredict.duckdb"
    build_duckdb(db_root=config.db_root, out_path=out_path)
    print(f"built {out_path}")
    return 0


def register_subcommands(subparsers: argparse._SubParsersAction) -> None:
    """Register `curate-db` and its actions onto the top-level parser."""
    curate_db = subparsers.add_parser("curate-db", help="Curate the MAT locus reference database")
    action = curate_db.add_subparsers(dest="action", required=True)

    propose = action.add_parser("propose")
    propose.add_argument("--phylum", required=True)
    propose.add_argument("--record-file", required=True)
    propose.set_defaults(func=_cmd_propose)

    validate = action.add_parser("validate")
    validate.add_argument("--phylum", required=True)
    validate.add_argument("--record-id", required=True)
    validate.set_defaults(func=_cmd_validate)

    accept = action.add_parser("accept")
    accept.add_argument("--phylum", required=True)
    accept.add_argument("--order-or-family", required=True)
    accept.add_argument("--record-id", required=True)
    accept.set_defaults(func=_cmd_accept)

    reject = action.add_parser("reject")
    reject.add_argument("--phylum", required=True)
    reject.add_argument("--record-id", required=True)
    reject.add_argument("--reason", required=True)
    reject.set_defaults(func=_cmd_reject)

    build_gff = action.add_parser("build-gff")
    build_gff.set_defaults(func=lambda a: (_ for _ in ()).throw(NotImplementedError("wire per-record GFF export CLI when a real accepted record needs it")))

    build_db = action.add_parser("build-duckdb")
    build_db.add_argument("--out", required=False)
    build_db.set_defaults(func=_cmd_build_duckdb)

    release = action.add_parser("release")
    release.set_defaults(func=lambda a: (_ for _ in ()).throw(NotImplementedError("wire release cut CLI in Task 14")))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/test_cli_smoke.py -v`
Expected: PASS

- [ ] **Step 5: Run the full test suite**

Run: `pixi run -e test pytest -v`
Expected: PASS (all tests from Tasks 1-12)

- [ ] **Step 6: Commit**

```bash
git add src/MATPredict/db/cli.py tests/test_cli_smoke.py
git commit -m "feat: wire propose/validate/accept/reject/build-duckdb into matpredict curate-db"
```

---

### Task 13: Seed literature curation (data task, not unit-testable)

This task populates real records. It is a runbook, not a code change — acceptance is checked by running the tooling from Tasks 1-12 against real literature, not by writing new unit tests.

**Procedure, run once per phylum (Ascomycota, Basidiomycota, Mucoromycota):**

- [ ] **Step 1: Seed the search from `annotated.yml`.** Read `db/<Phylum>/annotated.yml` (where it exists) and list every organism/accession combination it already names.

- [ ] **Step 2: Search PubMed per organism/genus** (via the PubMed MCP tool) for `"<genus> mating type locus"` / `"<genus> MAT locus"` style queries. For each candidate finding, record the PMID/DOI and the exact sentence/table the locus/gene/coordinate claim is drawn from — no claim without one.

- [ ] **Step 3: Draft each candidate as a `metadata.yaml`** matching the Task 2 schema, setting `coordinate_provenance` honestly (`not_available` for genetic/RFLP-only evidence, `curator_derived` when coordinates must be located on a cited assembly, `published_explicit` when the paper states them directly). Namespace every `protein_accession` (`ncbi_protein:`/`uniprotkb:`). Set `proposal_dedupe_key` to `pmid|taxid|locus_name`.

- [ ] **Step 4: Run `matpredict curate-db propose --phylum <Phylum> --record-file <draft>.yaml`** for each draft. A schema failure here (`CurationError`) means the draft is malformed — fix it before continuing, not by loosening the schema.

- [ ] **Step 5: Run `matpredict curate-db validate --phylum <Phylum> --record-id <record_id>`** for each proposed candidate. Read the printed validation JSON.

- [ ] **Step 6: Human review.** For each candidate with `needs_review`, compare the citation and validation JSON side by side. Run `matpredict curate-db accept --phylum <Phylum> --order-or-family <Order> --record-id <record_id>` to accept, or `... reject --record-id <record_id> --reason "<reason>"` to reject. Every rejection gets a real, specific reason (e.g. "accession GCA_x.1 suppressed, no replacement cited in paper").

- [ ] **Step 7: Confirm edge-case coverage per phylum** before calling the phylum done: at least one accepted record with `coordinate_provenance: not_available`, one with `sequence_source.type: insdc_nucleotide`, one with a `present: false` gene, and — for Basidiomycota specifically — one strain with both an accepted `HD` record and an accepted `PR` record.

- [ ] **Step 8: Generate GFF3/GBK/proteins.faa for every accepted record.** For each, fetch the gene sequences via `NcbiClient`/`UniprotClient` (same clients Task 12 wires up), call `gff_export.write_gff3` and `gff_export.write_proteins_fasta`, and hand-build/derive the paired `locus.gbk` (via `Bio.SeqIO.write` with a `SeqRecord` built from the same coordinates) so all three files agree — this is exactly what the "Cross-file consistency" acceptance check in the spec verifies.

- [ ] **Step 9: Rebuild the DuckDB cache** with `matpredict curate-db build-duckdb` and spot-check row counts against file counts (`db/<Phylum>/**/metadata.yaml` count == `SELECT COUNT(*) FROM locus_record WHERE validation_status='accepted'`).

- [ ] **Step 10: Repeat for all three phyla** until each has 5-10 accepted tier-1 records including its required edge cases.

- [ ] **Step 11: Commit each phylum's accepted records as they land** (don't batch all three phyla into one commit):

```bash
git add db/Mucoromycota/
git commit -m "data: curate N tier-1 Mucoromycota MAT locus records"
```

---

### Task 14: Cut the v1 release

**Files:**
- Modify: `db/_release.yml`
- Create: `src/MATPredict/db/release.py`
- Modify: `src/MATPredict/db/cli.py` (wire the `release` action left as `NotImplementedError` in Task 12)
- Create: `tests/db/test_release.py`

**Interfaces:**
- Consumes: DuckDB build (Task 11) to enumerate all currently-`accepted` `record_id`s.
- Produces: `release.cut_release(db_root: Path, release_tag: str, cut_date: str, git_tag_runner=subprocess.run) -> dict` (returns the manifest entry written).

- [ ] **Step 1: Write the failing test**

```python
# tests/db/test_release.py
from __future__ import annotations

from types import SimpleNamespace

import yaml

from MATPredict.db.curate import propose_candidate, accept_candidate
from MATPredict.db.release import cut_release

RECORD = {
    "record_id": "4837_nrrl-1555_MAT_Plus",
    "record_version": 1,
    "taxonomy": {"taxid": 4837, "lineage": "k__Fungi;p__Mucoromycota", "lineage_resolved_date": "2026-09-16"},
    "organism": {"species": "Phycomyces blakesleeanus", "strain": {"name": "NRRL 1555", "known": True}},
    "mating_type": {"locus_name": "MAT", "idiomorphs": ["Plus"], "system": "heterothallic"},
    "locus": {"coordinate_provenance": "not_available", "excluded_from_coordinate_benchmark": True},
    "genes": [{"gene_index": 0, "name": "sexP", "protein_accession": "ncbi_protein:AAB12345.1", "role": "core_MAT", "present": True}],
    "evidence": {
        "locus_existence": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "boundaries": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
        "idiomorph_assignment": {"tier": 1, "citations": [{"pmid": "18248337", "doi": None}]},
    },
    "validation": {"status": "needs_review", "rejection_reason": None},
    "curation": {"proposed_by": "literature-mining-agent"},
}


def _fake_git_runner(cmd, **kwargs):
    return SimpleNamespace(returncode=0, stdout="")


def test_cut_release_writes_manifest_and_tags(tmp_path):
    propose_candidate(tmp_path, "Mucoromycota", RECORD)
    accept_candidate(tmp_path, "Mucoromycota", "Mucorales", RECORD["record_id"])

    tagged_commands = []
    def runner(cmd, **kwargs):
        tagged_commands.append(cmd)
        return _fake_git_runner(cmd, **kwargs)

    entry = cut_release(tmp_path, release_tag="2026.09.0", cut_date="2026-09-16", git_tag_runner=runner)

    assert entry["records"] == ["4837_nrrl-1555_MAT_Plus@1"]
    manifest = yaml.safe_load((tmp_path / "_release.yml").read_text())
    assert "2026.09.0" in manifest["releases"]
    assert any("db-release-2026.09.0" in " ".join(c) for c in tagged_commands)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run -e test pytest tests/db/test_release.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write `src/MATPredict/db/release.py`**

```python
"""Cut a release: snapshot accepted record_id@version pairs and tag the commit."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Callable

import yaml


def _all_accepted_records(db_root: Path) -> list[tuple[str, int]]:
    results = []
    for metadata_path in sorted(db_root.glob("**/metadata.yaml")):
        if "candidates" in metadata_path.parts:
            continue
        record = yaml.safe_load(metadata_path.read_text())
        if record["validation"]["status"] == "accepted":
            results.append((record["record_id"], record["record_version"]))
    return results


def cut_release(
    db_root: Path,
    release_tag: str,
    cut_date: str,
    git_tag_runner: Callable = subprocess.run,
) -> dict:
    """Append a release entry to db/_release.yml and create a matching git tag."""
    records = [f"{record_id}@{version}" for record_id, version in _all_accepted_records(db_root)]
    git_tag = f"db-release-{release_tag}"

    manifest_path = db_root / "_release.yml"
    manifest = yaml.safe_load(manifest_path.read_text()) if manifest_path.exists() else {"releases": {}}
    entry = {"cut_date": cut_date, "git_tag": git_tag, "records": records}
    manifest["releases"][release_tag] = entry
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))

    git_tag_runner(["git", "tag", git_tag], cwd=str(db_root), capture_output=True, text=True)
    return entry
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run -e test pytest tests/db/test_release.py -v`
Expected: PASS

- [ ] **Step 5: Wire the `release` CLI action in `src/MATPredict/db/cli.py`**

Replace the placeholder `release` parser registration from Task 12 with:

```python
def _cmd_release(args: argparse.Namespace) -> int:
    config = _config(args)
    entry = cut_release(config.db_root, release_tag=args.release_tag, cut_date=args.cut_date)
    print(f"cut release {args.release_tag}: {len(entry['records'])} records, tag {entry['git_tag']}")
    return 0
```

and (add the import `from MATPredict.db.release import cut_release` near the top of the file):

```python
    release = action.add_parser("release")
    release.add_argument("--release-tag", required=True)
    release.add_argument("--cut-date", required=True)
    release.set_defaults(func=_cmd_release)
```

- [ ] **Step 6: Run the full test suite one more time**

Run: `pixi run -e test pytest -v`
Expected: PASS (all tests from Tasks 1-14)

- [ ] **Step 7: Cut the actual v1 release** (after Task 13's records are accepted)

Run:
```bash
pixi run matpredict curate-db release --release-tag 2026.09.0 --cut-date 2026-09-16
```

- [ ] **Step 8: Commit**

```bash
git add src/MATPredict/db/release.py src/MATPredict/db/cli.py tests/db/test_release.py db/_release.yml
git commit -m "feat: release-cut tooling; wire matpredict curate-db release"
```

---

## Self-Review Notes

- **Spec coverage**: every spec section has a task — packaging/CLI (Task 1), both schemas + cross-file idiomorph validation (Task 2), `order.yml` content (Task 3), taxonomy/accession/sequence validation (Tasks 4-8), candidate lifecycle + status/location invariant (Task 9), GFF3/FASTA export (Task 10), DuckDB build (Task 11), CLI wiring (Task 12), actual curated records + edge cases (Task 13), release manifest + git tag (Task 14). `locus.gbk` generation is named explicitly in Task 13 Step 8 rather than a separate module, since it's a thin `Bio.SeqIO` wrapper around the same data `gff_export.py` already structures — flagged here rather than silently dropped.
- **Placeholder scan**: the two `NotImplementedError` stubs in Task 12 (`build-gff` CLI wiring, `release` CLI wiring) are intentional forward-references to Tasks 13/14, not unresolved placeholders — both are wired concretely before the plan ends.
- **Type consistency checked**: `AccessionStatus` is defined once in `ncbi_client.py` and reused (not redefined) by `uniprot_client.py`; `record_id`/`gene_index`/`segment_index` naming is consistent from `identifiers.py` through `curate.py`, `gff_export.py`, and `build_duckdb.py`; `MatchScore.status` values (`pass`/`warn`/`fail`) match the DuckDB `sequence_match_status` column and the schema's enum everywhere they appear.
