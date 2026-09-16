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
