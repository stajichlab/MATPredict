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


def validate_gene_vocabulary(record: dict, order_doc: dict) -> list[str]:
    """Cross-check every present gene name in a record against its locus's declared `genes`.

    `detect.search._attribute` treats the family's `order.yml` gene list as a hard
    filter: a hit whose curated gene name is not declared there is dropped before
    clustering, scoring and reporting. A curated gene whose name is absent from the
    declaration is therefore silently unreachable -- it can never contribute a hit,
    a cluster membership or a report line in any genome.

    This check makes that failure loud at curation time. It is deliberately
    one-directional: it asserts every curated name is declared, and does NOT assert
    the converse (that every declared name has a curated protein). The converse is
    legitimately false today -- a family may declare genes belonging to an idiomorph
    that has not been curated yet (Ascomycota MATsc's MATA1/MATA2), or a
    `flanking_variable` gene that only some records carry (Mucoromycota algA/glrA).

    A gene marked `present: false` is skipped: it records a documented absence, not
    a sequence that search will ever match.
    """
    locus_name = record["mating_type"]["locus_name"]
    locus_entry = next((l for l in order_doc["loci"] if l["locus_name"] == locus_name), None)
    if locus_entry is None:
        return [f"no order.yml locus entry named '{locus_name}'"]

    # Aliases count as declared. A roster gene may collapse several curated
    # names onto one canonical name (`Family.gene_aliases`) when their proteins
    # are byte-identical and no homology score could ever separate them --
    # MFa1/MFa2/MFa3 -> MFa. The curated record keeps the gene name its
    # publication deposited, and `reference_fasta.build_reference_fasta`
    # rewrites it to the canonical name on the way into the search, so the hit
    # IS attributable and this check must not call it orphaned.
    declared = {
        name
        for gene in locus_entry.get("genes", [])
        for name in [gene["name"], *(gene.get("aliases") or [])]
    }
    errors: list[str] = []
    for gene in record.get("genes", []):
        if not gene.get("present", True):
            continue
        name = gene["name"]
        if name not in declared:
            errors.append(
                f"gene '{name}' (gene_index {gene.get('gene_index')}) is not declared in "
                f"order.yml locus '{locus_name}' genes {sorted(declared)}; "
                "detect.search._attribute would silently drop every hit to it"
            )
    return errors


def validate_gene_status(record: dict) -> list[str]:
    """Check the semantics of the optional gene-level `status` (`remnant|pseudogene`).

    `status` says what a recorded gene IS; it is orthogonal to `present`:

    * `status` + `present: false`  -- no model is recorded: no coordinates and no
      protein (a remnant known only from prose or from a genome comparison).
    * `status` + `present: true`   -- a model IS recorded: a protein accession or
      coordinates must be given, or the status describes nothing.
    * `remnant_note` without a `status` is an orphan note.

    The enum itself (and the rejection of `partial`) is the JSON schema's job;
    `status` is never for assembly or annotation truncation (that is
    `locus.core.completeness`). Returns one message per problem.
    """
    errors: list[str] = []
    for gene in record.get("genes", []):
        label = f"gene '{gene.get('name')}' (gene_index {gene.get('gene_index')})"
        status = gene.get("status")
        if status is None:
            if gene.get("remnant_note"):
                errors.append(f"{label} has a remnant_note but no status")
            continue
        has_model = bool(gene.get("protein_accession")) or (
            gene.get("start") is not None and gene.get("end") is not None
        )
        if gene.get("present", True):
            if not has_model:
                errors.append(
                    f"{label} has status '{status}' and present: true but records no model "
                    "(no protein_accession and no start/end); a gene with no model is present: false"
                )
        elif has_model:
            errors.append(
                f"{label} has status '{status}' and present: false but records a model "
                "(protein_accession or coordinates); present: false means no model is recorded"
            )
    return errors


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
