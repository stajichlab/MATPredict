"""Offline NCBI taxonomy: a slim table built from an NCBI taxdump.

`detect` needs three facts per taxid: the ancestor lineage (for
`taxonomic_scope` routing), the phylum name (phylum fallback) and the nuclear
genetic code. All three are in the taxdump's `nodes.dmp` (parent, rank,
genetic code id) and `names.dmp` (scientific name). This module writes them,
for every taxon, to one zstd-compressed text table and answers the three
questions from it, so a run needs no NCBI E-utilities call. Measured on the
2026-05 taxdump: 2.83 M taxa, 25 MB compressed.

Table format (UTF-8 text, zstd-compressed):

    #matpredict-taxonomy<TAB>v1<TAB>snapshot=<date><TAB>source=<file>
    N<TAB>taxid<TAB>parent<TAB>rank<TAB>genetic_code<TAB>scientific_name
    M<TAB>old_taxid<TAB>new_taxid        (merged.dmp)

Lookups mirror `NcbiClient` (efetch db=taxonomy) exactly:
- lineage: ancestor taxids, root-first, without the root (1) and without the
  taxid itself (efetch `LineageEx`);
- phylum: scientific name of the phylum-rank ANCESTOR (not the taxid itself);
- genetic code: the taxid's own nuclear code (`GCId`).
A merged taxid is answered for the taxid it was merged into, as efetch does.

Where the table is: `$MATPREDICT_TAXONOMY`. With `$MATPREDICT_OFFLINE=1` a
taxid missing from the table raises `TaxidNotInSnapshot` (the router records
it in `routing_error`); otherwise the caller falls back to E-utilities.
"""
from __future__ import annotations

import gzip
import io
import os
import subprocess
import tarfile
import zipfile
from array import array
from pathlib import Path

FORMAT_TAG = "#matpredict-taxonomy"
FORMAT_VERSION = "v1"


class TaxidNotInSnapshot(LookupError):
    """The taxid is not in the local taxonomy snapshot (newer, deleted or wrong)."""


# --------------------------------------------------------------------- zstd I/O

def _zstd_decompress(data: bytes) -> bytes:
    try:
        from compression import zstd  # Python >= 3.14
        return zstd.decompress(data)
    except ImportError:
        pass
    try:
        import zstandard
        return zstandard.ZstdDecompressor().decompressobj().decompress(data)
    except ImportError:
        pass
    return subprocess.run(["zstd", "-dc"], input=data, capture_output=True, check=True).stdout


def _zstd_compress(data: bytes, level: int = 19) -> bytes:
    try:
        from compression import zstd  # Python >= 3.14
        return zstd.compress(data, level=level)
    except ImportError:
        pass
    try:
        import zstandard
        return zstandard.ZstdCompressor(level=level).compress(data)
    except ImportError:
        pass
    return subprocess.run(["zstd", "-q", f"-{level}", "-c"], input=data, capture_output=True,
                          check=True).stdout


def _open_text(path: Path):
    """Open a table (plain or zstd) for streamed text reading."""
    with open(path, "rb") as fh:
        magic = fh.read(4)
    if magic != b"\x28\xb5\x2f\xfd":
        return open(path, encoding="utf-8")
    try:
        from compression import zstd  # Python >= 3.14: streamed
        return zstd.open(path, "rt", encoding="utf-8")
    except ImportError:
        pass
    try:
        import zstandard  # Python < 3.14 (a dependency there): streamed
        raw = open(path, "rb")
        return io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(raw, closefd=True),
                                encoding="utf-8")
    except ImportError:
        pass
    return io.StringIO(_zstd_decompress(Path(path).read_bytes()).decode("utf-8"))


def _decode(name: str, data: bytes) -> str:
    """Bytes of a dump file, plain or compressed (detected from the magic bytes)."""
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    elif data[:4] == b"\x28\xb5\x2f\xfd":
        data = _zstd_decompress(data)
    return data.decode("utf-8")


def _read_dump_files(taxdump: Path) -> dict[str, str]:
    """`nodes.dmp`, `names.dmp` and `merged.dmp` (optional) from a taxdump.

    `taxdump` is a directory (files plain, `.gz` or `.zst`), an NCBI archive
    `.zip` (`taxdmp_YYYY-MM-DD.zip`) or a `.tar.gz` (`taxdump.tar.gz`).
    """
    wanted = ("nodes.dmp", "names.dmp", "merged.dmp")
    out: dict[str, str] = {}
    taxdump = Path(taxdump)
    if taxdump.is_dir():
        for name in wanted:
            for cand in (name, name + ".gz", name + ".zst"):
                p = taxdump / cand
                if p.exists():
                    out[name] = _decode(cand, p.read_bytes())
                    break
    elif zipfile.is_zipfile(taxdump):
        with zipfile.ZipFile(taxdump) as zf:
            for name in wanted:
                if name in zf.namelist():
                    out[name] = _decode(name, zf.read(name))
    elif tarfile.is_tarfile(taxdump):
        with tarfile.open(taxdump) as tf:
            for name in wanted:
                try:
                    member = tf.getmember(name)
                except KeyError:
                    continue
                out[name] = _decode(name, tf.extractfile(member).read())
    else:
        raise ValueError(f"not a taxdump directory, .zip or .tar.gz: {taxdump}")
    for name in ("nodes.dmp", "names.dmp"):
        if name not in out:
            raise FileNotFoundError(f"{name} not found in {taxdump}")
    return out


def _dmp_rows(text: str):
    """Fields of NCBI `.dmp` rows (`\\t|\\t` separated, `\\t|` at line end)."""
    for line in text.splitlines():
        if line:
            yield line.rstrip("\t|").split("\t|\t")


def build_table(taxdump: Path, out_path: Path, snapshot: str, source: str | None = None) -> dict:
    """Write the slim table for every taxon in `taxdump`. Returns counts."""
    files = _read_dump_files(taxdump)
    names: dict[str, str] = {}
    for f in _dmp_rows(files["names.dmp"]):
        if f[3] == "scientific name":
            names[f[0]] = f[1]
    buf = io.StringIO()
    max_taxid = max(int(f[0]) for f in _dmp_rows(files["nodes.dmp"]))
    buf.write(f"{FORMAT_TAG}\t{FORMAT_VERSION}\tsnapshot={snapshot}\t"
              f"source={source or Path(taxdump).name}\tmax_taxid={max_taxid}\n")
    n_nodes = 0
    for f in _dmp_rows(files["nodes.dmp"]):
        # nodes.dmp: tax_id, parent tax_id, rank, embl code, division id,
        # inherited div flag, genetic code id, ...
        buf.write(f"N\t{f[0]}\t{f[1]}\t{f[2]}\t{f[6]}\t{names.get(f[0], '')}\n")
        n_nodes += 1
    n_merged = 0
    for f in _dmp_rows(files.get("merged.dmp", "")):
        buf.write(f"M\t{f[0]}\t{f[1]}\n")
        n_merged += 1
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_name(out_path.name + ".tmp")
    tmp.write_bytes(_zstd_compress(buf.getvalue().encode("utf-8")))
    tmp.replace(out_path)
    return {"taxa": n_nodes, "merged": n_merged, "snapshot": snapshot}


# ----------------------------------------------------------------------- lookups

class LocalTaxonomy:
    """Lineage, phylum name and genetic code from a slim table.

    Memory: parent and genetic code are int arrays indexed by taxid (about 4
    bytes x max taxid each); names are kept only for phylum-rank taxa.
    """

    def __init__(self, path: Path):
        self.path = Path(path)
        self._is_phylum: set[int] = set()
        self._phylum_name: dict[int, str] = {}
        self._merged: dict[int, int] = {}
        self.n_taxa = 0
        with _open_text(self.path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if header[0] != FORMAT_TAG or header[1] != FORMAT_VERSION:
                raise ValueError(f"not a {FORMAT_TAG} {FORMAT_VERSION} table: {self.path}")
            meta = dict(h.split("=", 1) for h in header[2:] if "=" in h)
            self.snapshot = meta.get("snapshot", "unknown")
            self.source = meta.get("source", "")
            size = int(meta.get("max_taxid", 0)) + 1
            # Streamed line by line into int arrays indexed by taxid, so memory
            # stays near 2 x 4 bytes x max taxid (measured below in the tests'
            # docstring) instead of a list of 3 M parsed rows.
            self._parent = array("i", [-1]) * size
            self._code = array("i", [0]) * size
            for line in fh:
                kind, _, rest = line.partition("\t")
                if kind == "N":
                    taxid, parent, rank, code, name = rest.rstrip("\n").split("\t", 4)
                    t = int(taxid)
                    if t >= len(self._parent):  # table without max_taxid
                        grow = t + 1 - len(self._parent)
                        self._parent.extend([-1] * grow)
                        self._code.extend([0] * grow)
                    self._parent[t] = int(parent)
                    self._code[t] = int(code) if code else 0
                    if rank == "phylum":
                        self._is_phylum.add(t)
                        self._phylum_name[t] = name
                    self.n_taxa += 1
                elif kind == "M":
                    old, new = rest.rstrip("\n").split("\t")
                    self._merged[int(old)] = int(new)

    def _resolve(self, taxid: int) -> int:
        taxid = self._merged.get(int(taxid), int(taxid))
        if taxid <= 0 or taxid >= len(self._parent) or self._parent[taxid] < 0:
            raise TaxidNotInSnapshot(
                f"taxid {taxid} is not in the local NCBI taxonomy snapshot {self.snapshot}")
        return taxid

    def __contains__(self, taxid: int) -> bool:
        try:
            self._resolve(taxid)
            return True
        except TaxidNotInSnapshot:
            return False

    def lineage(self, taxid: int) -> list[int]:
        """Ancestors root-first, without the root (1) and without `taxid`."""
        t = self._resolve(taxid)
        out = []
        x = self._parent[t]
        while x > 1:
            out.append(x)
            nxt = self._parent[x]
            if nxt == x:
                break
            x = nxt
        return out[::-1]

    def phylum_name(self, taxid: int) -> str | None:
        for a in self.lineage(taxid):
            if a in self._is_phylum:
                return self._phylum_name.get(a) or None
        return None

    def genetic_code(self, taxid: int) -> int | None:
        code = self._code[self._resolve(taxid)]
        return code or None


_loaded: LocalTaxonomy | None = None
_loaded_path: str | None = None


def configured_path() -> Path | None:
    """`$MATPREDICT_TAXONOMY` when set and the file exists, else None."""
    p = os.environ.get("MATPREDICT_TAXONOMY")
    if p and Path(p).is_file():
        return Path(p)
    return None


def offline() -> bool:
    """True when `$MATPREDICT_OFFLINE` is 1/true/yes: never call E-utilities."""
    return os.environ.get("MATPREDICT_OFFLINE", "").strip().lower() in ("1", "true", "yes")


def get_local() -> LocalTaxonomy | None:
    """The configured table, loaded once per process (None if not configured)."""
    global _loaded, _loaded_path
    p = configured_path()
    if p is None:
        return None
    if _loaded is None or _loaded_path != str(p):
        _loaded, _loaded_path = LocalTaxonomy(p), str(p)
    return _loaded


def source_label() -> str:
    """What answers taxonomy lookups in this process, for the report."""
    local = get_local()
    if local is not None:
        tail = "" if offline() else "; NCBI E-utilities for taxids not in it"
        return f"local NCBI taxonomy snapshot {local.snapshot}{tail}"
    if offline():
        return "none (offline, no local taxonomy table)"
    return "NCBI E-utilities"
