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
