"""Candidate lifecycle: propose, accept, reject, edit-in-place. Enforces the
status/directory-location invariant from the spec: accepted iff under db/<Phylum>/."""
from __future__ import annotations

import shutil
from pathlib import Path

import yaml

from MATPredict.db.schema import validate_idiomorphs, validate_metadata


class CurationError(Exception):
    """Raised when a curation operation would violate a database invariant."""


def _candidate_dir(db_root: Path, phylum: str, record_id: str) -> Path:
    return db_root / "candidates" / phylum / record_id


def _load_order_doc(db_root: Path, phylum: str) -> dict:
    return yaml.safe_load((db_root / phylum / "order.yml").read_text())


def _validate_against_schema_and_order(db_root: Path, phylum: str, record: dict) -> None:
    errors = validate_metadata(record)
    order_doc = _load_order_doc(db_root, phylum)
    errors += validate_idiomorphs(record, order_doc)
    if errors:
        raise CurationError(f"invalid candidate record: {'; '.join(errors)}")


def _check_dedupe_key_not_rejected(db_root: Path, phylum: str, record: dict) -> None:
    dedupe_key = record.get("curation", {}).get("proposal_dedupe_key")
    if not dedupe_key:
        return

    candidates_root = db_root / "candidates" / phylum
    for metadata_path in candidates_root.glob("*/metadata.yaml"):
        existing = yaml.safe_load(metadata_path.read_text())
        existing_key = existing.get("curation", {}).get("proposal_dedupe_key")
        existing_status = existing.get("validation", {}).get("status")
        if existing_key == dedupe_key and existing_status == "rejected":
            reason = existing.get("validation", {}).get("rejection_reason")
            raise CurationError(
                f"a candidate with dedupe key '{dedupe_key}' was already rejected: {reason}"
            )


def propose_candidate(db_root: Path, phylum: str, record: dict) -> Path:
    """Write a new draft record under db/candidates/<phylum>/<record_id>/metadata.yaml."""
    _validate_against_schema_and_order(db_root, phylum, record)
    _check_dedupe_key_not_rejected(db_root, phylum, record)

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
    _validate_against_schema_and_order(db_root, phylum, record)

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
