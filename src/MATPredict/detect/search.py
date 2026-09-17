"""diamond (fast path) and exonerate (fallback path) search wrappers.

Reference-protein headers are `{record_id}|gene{gene_index}|{gene_name}`
(the form `reference_fasta.build_reference_fasta` rewrites
`gff_export.write_proteins_fasta`'s real headers into) -- parsed back out
here to recover which curated record and gene a hit corresponds to.

Attribution (which `(phylum, locus_name)` family a hit belongs to) is keyed
on that `record_id`, never on the bare gene name. Gene names are NOT unique
across families in the real database: `pheromone`/`pheromone_receptor` occur
in `Basidiomycota:PR`, `Basidiomycota:Balpha` and `Basidiomycota:Bbeta`;
`Y`/`Z` in `Basidiomycota:Aalpha` and `Basidiomycota:Abeta`;
`matPc`/`matPi`/`matMc`/`matMi` across `Ascomycota:PM`, `Ascomycota:mat2`
and `Ascomycota:mat3`; `sla2` in `Ascomycota:MATsc` and `Ascomycota:MATyl`.
A bare-gene-name lookup keeps only the last family written for each of those
names, so every other family sharing the name becomes unreachable and its
evidence is misattributed. A curated record, by contrast, belongs to exactly
one family (`family_registry.load_record_families`), so the record the hit's
reference sequence came from identifies the family exactly.

--outfmt design (fast path): diamond blastp is a protein-vs-protein search,
so it has no native genomic-coordinate columns. The predicted proteome
passed as `proteome_fasta` is expected to carry its gene's genomic location
on the FASTA defline after the id, as `contig:start-end:strand`
(a convention already used by several gene-prediction pipelines' protein
FASTA output). Requesting diamond's `qtitle` field (everything on the query
defline) recovers that location string, which we then pick apart ourselves.
Verified against diamond v2.2.6: `qtitle` is the WHOLE query defline,
including the query id itself (`q1 contigA:100-400:+`), not only the part
after it -- so the id is stripped and the location token is searched for
among the remaining whitespace-separated fields. `scovhsp` (percentage of the *subject* -- i.e. the curated
reference protein -- covered by the alignment) is requested alongside
`pident` because the spec's report section asks for per-gene identity *and*
coverage against the matched reference protein. The outfmt is therefore:

    6 qseqid sseqid pident scovhsp qtitle

diamond database: `search_fast_path` runs `diamond makedb` on the reference
FASTA into a temporary directory and searches against the resulting `.dmnd`.
diamond v2.2.6 was observed to also accept a plain FASTA as `--db` (building
an index on the fly), but that is a convenience of recent versions, not the
documented contract -- `makedb` is the interface diamond actually specifies
for `blastp --db`, and building it once up front is also what makes the
reference set reusable across queries.

Window-restriction (genomic fallback): rather than relying on an exonerate
flag to restrict the search region (exonerate has no first-class "search
only this coordinate range of this multi-contig file" option), we actually
slice the requested contig region out of `genome_fasta` into a small
temporary FASTA file and pass *that* as `--target`. Exonerate's reported
hit coordinates are then local to the slice (1-based from the slice start),
so we add the window's start offset back on to recover real genome
coordinates before returning `SearchHit`s.

Relaxation (genomic second pass): exonerate's `--score` is an absolute
score threshold (its protein2genome default is 100). Lowering it makes the
search MORE permissive. `--percent` is the opposite knob -- it discards
matches scoring below that percentage of the query's maximal score, and its
default is 0 -- so adding `--percent 50` to the "relaxed" pass (as this
module previously did) made the second pass *stricter* than the standard
pass, inverting the spec's intent. The relaxed pass therefore lowers
`--score` instead, and the standard pass states exonerate's own default
explicitly so the difference between the two is visible in the command.
"""
from __future__ import annotations

import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from Bio import SeqIO

from MATPredict.detect.family_registry import Family, FamilyKey

_DIAMOND_OUTFMT = ["6", "qseqid", "sseqid", "pident", "scovhsp", "qtitle"]

#: exonerate protein2genome's own default score threshold.
STANDARD_EXONERATE_SCORE = 100
#: Lower (more permissive) threshold used by the relaxed genomic second pass.
RELAXED_EXONERATE_SCORE = 50

PROTEOME_DEFLINE_FORMAT = "contig:start-end:strand"
_LOCATION_RE = re.compile(r"^(?P<contig>.+):(?P<start>\d+)-(?P<end>\d+):(?P<strand>[+-])$")


class SearchToolError(RuntimeError):
    """A search binary exited non-zero, so its empty output is not a real negative."""


class ProteomeDeflineError(ValueError):
    """A predicted-proteome defline did not carry the required genomic location."""


@dataclass(frozen=True)
class SearchHit:
    family_key: FamilyKey
    gene_name: str
    role: str  # core_MAT | flanking_conserved | flanking_variable
    contig: str
    start: int  # 1-based
    end: int  # 1-based, inclusive
    strand: str  # "+" | "-"
    identity: float
    reference_record_id: str  # which curated record's protein this matched
    method: str  # "diamond_proteome" | "exonerate_genome" | "exonerate_genome_relaxed"
    coverage: float | None = None  # % of the matched reference protein covered; None when unknown


def _roles_by_family(families: list[Family]) -> dict[FamilyKey, dict[str, str]]:
    """family_key -> {gene name: role} for every attempted family.

    Scoped per family on purpose: the same gene name legitimately carries a
    role in more than one family, and each family's own table is the only
    correct source for its own genes.
    """
    return {family.key: {g["name"]: g["role"] for g in family.genes} for family in families}


def _attribute(
    record_id: str,
    gene_name: str,
    record_families: dict[str, FamilyKey],
    roles_by_family: dict[FamilyKey, dict[str, str]],
) -> tuple[FamilyKey, str] | None:
    """Resolve a reference hit to (family_key, role), or None if it should be dropped.

    A hit is dropped when its curated record is unknown, when that record's
    family is not among the attempted families, or when the gene name is not
    one that family actually expects.
    """
    family_key = record_families.get(record_id)
    if family_key is None:
        return None
    roles = roles_by_family.get(family_key)
    if roles is None:
        return None
    role = roles.get(gene_name)
    if role is None:
        return None
    return family_key, role


def _parse_reference_header(sseqid: str) -> tuple[str, str]:
    """`{record_id}|gene{N}|{gene_name}` -> (record_id, gene_name)."""
    record_id, _, gene_name = sseqid.split("|")
    return record_id, gene_name


def _parse_proteome_location(qtitle: str, qseqid: str | None = None) -> tuple[str, int, int, str]:
    """Pull `contig:start-end:strand` out of a query defline.

    `qtitle` is diamond's whole query defline, id included, so the id is
    dropped and every remaining whitespace-separated token is tested against
    the location format (a proteome may carry a free-text description after
    the location). The contig name itself may contain ':', so only the last
    two ':' separators are structural.

    Raises ProteomeDeflineError naming both the required format and the
    offending defline. Predicted-proteome FASTA headers from NCBI, AUGUSTUS
    and BRAKER do NOT carry this location by default, so a bare tuple-unpack
    ValueError here would surface as an unexplained crash on the very first
    real proteome anyone supplies. Failing loudly (rather than skipping the
    line) is deliberate: silently dropping unparseable deflines would turn a
    malformed input into a confident "no MAT locus found", which is the exact
    class of false negative this pipeline exists to avoid.
    """
    tokens = qtitle.strip().split()
    if qseqid and tokens and tokens[0] == qseqid:
        tokens = tokens[1:]
    for token in tokens:
        match = _LOCATION_RE.match(token)
        if match:
            return (
                match["contig"],
                int(match["start"]),
                int(match["end"]),
                match["strand"],
            )
    raise ProteomeDeflineError(
        f"predicted-proteome defline must carry the gene's genomic location as "
        f"'{PROTEOME_DEFLINE_FORMAT}' (e.g. '>gene1 contig_3:10500-11200:+'); "
        f"no such field in: {qtitle!r}"
    )


def _run_checked(runner: Callable, cmd: list[str]):
    """Run cmd and raise SearchToolError if it exits non-zero.

    Without this check a missing binary, a malformed database or a crash all
    produce empty stdout, which the parsers below read as "zero hits" -- a
    tool failure and a genuine negative result would be indistinguishable.
    """
    result = runner(cmd, capture_output=True, text=True)
    returncode = getattr(result, "returncode", 0)
    if returncode != 0:
        stderr = (getattr(result, "stderr", "") or "").strip()
        raise SearchToolError(
            f"`{' '.join(cmd)}` exited {returncode}"
            + (f": {stderr}" if stderr else "")
        )
    return result


def build_diamond_db(reference_fasta: Path, db_path: Path, runner: Callable = subprocess.run) -> Path:
    """Run `diamond makedb` on reference_fasta, returning the `.dmnd` path.

    `diamond blastp --db` requires diamond's own binary database format; it
    does not accept a plain FASTA. diamond appends `.dmnd` to the `--db`
    prefix it is given, so the prefix passed here is stripped of that suffix
    and the suffixed path is what gets returned.
    """
    prefix = db_path.with_suffix("") if db_path.suffix == ".dmnd" else db_path
    _run_checked(runner, ["diamond", "makedb", "--in", str(reference_fasta), "--db", str(prefix)])
    return prefix.with_suffix(".dmnd")


def search_fast_path(
    proteome_fasta: Path,
    families: list[Family],
    reference_fasta: Path,
    record_families: dict[str, FamilyKey],
    runner: Callable = subprocess.run,
) -> list[SearchHit]:
    """Search an existing predicted proteome against curated reference proteins via diamond blastp."""
    roles_by_family = _roles_by_family(families)

    with tempfile.TemporaryDirectory() as tmp_dir_name:
        diamond_db = build_diamond_db(
            reference_fasta, Path(tmp_dir_name) / "reference.dmnd", runner=runner
        )
        cmd = [
            "diamond", "blastp",
            "--query", str(proteome_fasta),
            "--db", str(diamond_db),
            "--outfmt", *_DIAMOND_OUTFMT,
        ]
        result = _run_checked(runner, cmd)

        hits: list[SearchHit] = []
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            qseqid, sseqid, pident, scovhsp, qtitle = line.split("\t")
            record_id, gene_name = _parse_reference_header(sseqid)
            attribution = _attribute(record_id, gene_name, record_families, roles_by_family)
            if attribution is None:
                continue
            family_key, role = attribution
            contig, start, end, strand = _parse_proteome_location(qtitle, qseqid)
            hits.append(
                SearchHit(
                    family_key=family_key,
                    gene_name=gene_name,
                    role=role,
                    contig=contig,
                    start=start,
                    end=end,
                    strand=strand,
                    identity=float(pident),
                    reference_record_id=record_id,
                    method="diamond_proteome",
                    coverage=float(scovhsp),
                )
            )
        return hits


def _extract_window(genome_fasta: Path, window: tuple[str, int, int], tmp_dir: Path) -> Path:
    """Slice (contig, start, end) (1-based, inclusive) out of genome_fasta into its own FASTA file.

    The slice keeps the original contig name as its record id, so exonerate's
    GFF output still reports the real contig name -- only start/end need
    re-basing back onto genome coordinates afterward.
    """
    contig, start, end = window
    index = SeqIO.index(str(genome_fasta), "fasta")
    try:
        record = index[contig]
    finally:
        index.close()
    sliced = record[start - 1 : end]
    sliced.id = contig
    sliced.description = ""
    window_fasta = tmp_dir / f"window_{contig}_{start}_{end}.fasta"
    SeqIO.write(sliced, window_fasta, "fasta")
    return window_fasta


def search_genomic(
    genome_fasta: Path,
    families: list[Family],
    reference_fasta: Path,
    record_families: dict[str, FamilyKey],
    relaxed: bool = False,
    window: tuple[str, int, int] | None = None,
    runner: Callable = subprocess.run,
) -> list[SearchHit]:
    """Spliced protein-to-genome search via exonerate --model protein2genome.

    `window` restricts the search to (contig, start, end) (1-based, inclusive)
    -- used for the flanking-anchored second pass -- by slicing that region
    out of `genome_fasta` into a temporary FASTA and searching against just
    that slice; hit coordinates are re-based back onto genome coordinates
    before being returned. `relaxed=True` lowers exonerate's absolute score
    threshold (`--score`) for that same second pass, making it genuinely more
    permissive than the standard pass.
    """
    roles_by_family = _roles_by_family(families)

    with tempfile.TemporaryDirectory() as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        offset = 0
        target_fasta = genome_fasta
        if window is not None:
            target_fasta = _extract_window(genome_fasta, window, tmp_dir)
            offset = window[1] - 1  # window start is 1-based; local coords start at 1

        score = RELAXED_EXONERATE_SCORE if relaxed else STANDARD_EXONERATE_SCORE
        cmd = [
            "exonerate", "--model", "protein2genome",
            "--query", str(reference_fasta),
            "--target", str(target_fasta),
            "--score", str(score),
            "--showtargetgff", "yes",
            "--showalignment", "no",
        ]

        result = _run_checked(runner, cmd)

        hits: list[SearchHit] = []
        for line in result.stdout.splitlines():
            if "\tgene\t" not in line:
                continue
            fields = line.split("\t")
            contig, _src, _feat, start, end, _score, strand, _frame, attrs = fields
            # exonerate GFF attrs carry the query id under "sequence <id>" and
            # the real percent identity under "identity <value>", e.g.:
            #   gene_id 1 ; sequence rec1|gene0|mfa1 ; gene_orientation . ; identity 100.00 ; similarity 100.00
            attr_parts = attrs.split(" ; ")
            query_id = next(
                part.split(" ")[1] for part in attr_parts if part.startswith("sequence ")
            )
            identity = 0.0
            for part in attr_parts:
                if part.startswith("identity "):
                    identity = float(part.split(" ")[1])
                    break
            record_id, gene_name = _parse_reference_header(query_id)
            attribution = _attribute(record_id, gene_name, record_families, roles_by_family)
            if attribution is None:
                continue
            family_key, role = attribution
            method = "exonerate_genome_relaxed" if relaxed else "exonerate_genome"
            hits.append(
                SearchHit(
                    family_key=family_key,
                    gene_name=gene_name,
                    role=role,
                    contig=contig,
                    start=int(start) + offset,
                    end=int(end) + offset,
                    strand=strand,
                    identity=identity,
                    reference_record_id=record_id,
                    method=method,
                    coverage=None,
                )
            )
        return hits
