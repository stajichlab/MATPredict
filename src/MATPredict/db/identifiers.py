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
