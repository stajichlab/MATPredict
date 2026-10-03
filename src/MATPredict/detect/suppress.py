"""Genomes the curator has suppressed: never run, and counted when skipped.

Curator's ruling 2026-09-26. Two lists are merged:

* the BFD list, `Fungi_BFD/data/curation/suppress.txt`, shared with the BFD,
  funannotate and ANI pipelines -- one `ASMID,REASON` row per genome, with a
  header. Too-small assemblies, contaminated MAGs, Microsporidia, and the nine
  PRJEB104476 records that are single Sanger amplicons deposited as "Complete
  Genome" (GCA_986280975.1 is 655 bp and reached a panel as an empty FASTA);
* a MATPredict-only list, `suppress.txt` in the curated database root, for
  genomes wrong for MAT detection but fine for BFD. It lives beside `db/` so a
  frozen `run-<sha>` worktree carries the list it was launched with, exactly
  as `MATPREDICT_DB_ROOT` already pins the database.

Both formats are read by one parser: blank lines and `#` comments are
skipped, the first comma- or tab-separated field is the ASMID, the rest is the
reason, and a header row whose first field is `ASMID` is ignored.

A row is matched by its full ASMID (`GCA_986280975.1_DAH1005FM`) or by the
bare accession (`GCA_986280975.1`), since some runners carry only the latter.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

from MATPredict import logger

#: The shared BFD list. `MATPREDICT_BFD_SUPPRESS` overrides it.
BFD_SUPPRESS_PATH = Path(
    "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt"
)
#: The MATPredict-only list's name inside the curated database root.
LOCAL_SUPPRESS_FILENAME = "suppress.txt"

_ACCESSION = re.compile(r"^(GC[AF]_\d+\.\d+)")


def default_suppress_paths(db_root: Path) -> list[Path]:
    """The BFD list, then the MATPredict-only list in `db_root`."""
    bfd = Path(os.environ.get("MATPREDICT_BFD_SUPPRESS", BFD_SUPPRESS_PATH))
    return [bfd, Path(db_root) / LOCAL_SUPPRESS_FILENAME]


def load_suppress_list(paths: list[Path]) -> dict[str, str]:
    """{ASMID: reason} merged over `paths`. A missing file is logged and
    skipped: the BFD list is not reachable off /bigdata, and a run there
    must not fail for it."""
    suppressed: dict[str, str] = {}
    for path in paths:
        path = Path(path)
        if not path.is_file():
            logger.warning("suppress list %s not found; skipped", path)
            continue
        for line in path.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            fields = [f.strip() for f in re.split(r"[,\t]", line, maxsplit=1)]
            asmid = fields[0]
            if not asmid or asmid.upper() == "ASMID":
                continue
            suppressed[asmid] = fields[1] if len(fields) > 1 else ""
    return suppressed


def _keys(asmid: str) -> set[str]:
    match = _ACCESSION.match(asmid)
    return {asmid, match.group(1)} if match else {asmid}


def filter_rows(rows: list[str], suppressed: dict[str, str]) -> tuple[list[str], list[str]]:
    """Split panel rows (`ASMID<TAB>...`) into (kept, skipped)."""
    blocked: set[str] = set()
    for asmid in suppressed:
        blocked |= _keys(asmid)
    kept, skipped = [], []
    for row in rows:
        asmid = re.split(r"[\t ,]", row.strip(), maxsplit=1)[0]
        (skipped if asmid and _keys(asmid) & blocked else kept).append(row)
    return kept, skipped
