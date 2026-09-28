"""Pheromone precursors found by a strict-CAAX ORF scan beside a receptor.

Curator's ruling 2026-09-27 (receptor queue, option ii). A B/PR mating locus is
often several pheromone-receptor copies plus tiny pheromone precursors
(25-80 aa) that tblastn cannot find. Every receptor copy carries the family's
one receptor gene name, so such a locus shows ONE distinct gene and never
clears the 2-distinct-gene admission floor: the Russula nobilis record
(three tandem receptors + a 25 aa CVVG precursor) called nothing, not even its
own genome (results/2026-09-27_russulaceae_receptor/).

The fix follows the measured signal. Within `window_bp` of each receptor hit,
a stop-anchored six-frame scan looks for a short ORF (in-frame Met
`min_codons`-`max_codons` upstream of the stop) whose last four residues match
the family's CAAX motif. A found ORF becomes a hit of the family's precursor
gene (method `caax_scan`), which the admission floor and the modelled-gene bar
count as a second gene.

Measured basis (results/2026-09-27_pheromone_positional/NOTE.md): the strict
motif C[VI][IV][AVMG] within 10 kb flagged 6/9 known mating receptors and 0/25
other STE3 copies, against 2.5% of random windows of the same size; the
textbook CAAX alphabet C[AVLIM][AVLIM]X hit 42% of random windows and is
useless. Sporidiobolales/Rhodotorula precursors end CTxA (e.g. CTIA/CTVA), so
the motif is a per-family roster parameter; the scan is NOT enabled for redPR
until a lineage motif is measured.

Confidence (`pipeline._build`): a scan precursor satisfies the gene-count
rules but cannot by itself raise a call to `high` -- it is a motif match, not
a homology model. A call whose modelled genes reach the bar only through the
scan is capped at medium.
"""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

from Bio.Seq import Seq

from MATPredict.detect.search import SearchHit

#: `SearchHit.method` / `GeneEvidence.method` of a scan precursor.
CAAX_METHOD = "caax_scan"
#: `GeneEvidence.status` of a scan precursor: a distinct evidence type, never
#: a polished homology model and never `unpolished`.
STATUS_CAAX_ORF = "caax_orf"

DEFAULT_MOTIF = "C[VI][IV][AVMG]"
DEFAULT_WINDOW_BP = 10_000
DEFAULT_MIN_CODONS = 20
DEFAULT_MAX_CODONS = 130

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def scan_config(spec: dict | None) -> dict | None:
    """A roster `pheromone_precursor_scan` block with defaults filled in, or None."""
    if not spec:
        return None
    return {
        "gene": spec["gene"],
        "receptor_genes": list(spec["receptor_genes"]),
        "motif": spec.get("motif", DEFAULT_MOTIF),
        "window_bp": int(spec.get("window_bp", DEFAULT_WINDOW_BP)),
        "min_codons": int(spec.get("min_codons", DEFAULT_MIN_CODONS)),
        "max_codons": int(spec.get("max_codons", DEFAULT_MAX_CODONS)),
    }


@dataclass(frozen=True)
class CaaxOrf:
    start: int  # 1-based, inclusive, in the scanned sequence (Met .. stop codon)
    end: int
    strand: str
    length_aa: int  # Met .. last CAAX residue, stop not counted
    motif: str


def scan_caax_orfs(
    seq: str,
    motif: str = DEFAULT_MOTIF,
    min_codons: int = DEFAULT_MIN_CODONS,
    max_codons: int = DEFAULT_MAX_CODONS,
    genetic_code: int = 1,
) -> list[CaaxOrf]:
    """Every stop-anchored ORF in `seq` (both strands) that ends in `motif`.

    A candidate is a stop codon whose upstream in-frame segment ends in the
    motif and holds a Met `min_codons`-`max_codons` codons before the stop
    (the most upstream such Met starts the ORF), exactly the rule of
    results/2026-09-27_pheromone_positional/scan_genome.py.
    """
    tail = re.compile(f"(?:{motif})$")
    out: list[CaaxOrf] = []
    n = len(seq)
    rc = seq.translate(_COMP)[::-1]
    for strand, s in (("+", seq), ("-", rc)):
        for frame in range(3):
            sub = s[frame:]
            sub = sub[: len(sub) - len(sub) % 3]
            if not sub:
                continue
            prot = str(Seq(sub).translate(table=genetic_code))
            seg_start = 0
            for m in re.finditer(r"\*", prot):
                seg = prot[seg_start:m.start()]
                seg_begin = seg_start
                seg_start = m.end()
                if len(seg) < min_codons or not tail.search(seg):
                    continue
                lo = max(0, len(seg) - max_codons)
                hi = len(seg) - min_codons + 1
                met = seg.find("M", lo, hi)
                if met < 0:
                    continue
                aa_begin = seg_begin + met          # index of Met in prot
                aa_stop = m.start()                 # index of the stop in prot
                nt_first = frame + 3 * aa_begin     # 0-based in s
                nt_last = frame + 3 * aa_stop + 2   # 0-based, last base of stop
                if strand == "+":
                    start, end = nt_first + 1, nt_last + 1
                else:
                    start, end = n - nt_last, n - nt_first
                out.append(CaaxOrf(start, end, strand, aa_stop - aa_begin, seg[-4:]))
    return sorted(out, key=lambda o: (o.start, o.end, o.strand))


def _read_contigs(genome_fasta: Path, wanted: set[str]) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    opener = gzip.open if str(genome_fasta).endswith(".gz") else open
    try:
        with opener(genome_fasta, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    if name in wanted:
                        seqs[name] = []
                elif name in seqs:
                    seqs[name].append(line.strip())
    except OSError:
        return {}
    return {k: "".join(v).upper() for k, v in seqs.items()}


def caax_precursor_hits(
    genome_fasta: Path,
    families: list,
    hits: list[SearchHit],
    genetic_code: int = 1,
) -> list[SearchHit]:
    """Scan precursors for every family with a `pheromone_precursor_scan`.

    One hit per distinct ORF (an ORF inside two receptors' windows is one
    hit), each carrying the family's precursor gene name and role, identity
    0.0 (there is no alignment), `align_length_aa` = ORF length, and the motif
    in `reference_record_id` ("caax_scan:CVIA").
    """
    jobs = []
    for family in families:
        cfg = scan_config(getattr(family, "pheromone_precursor_scan", None))
        if cfg is None:
            continue
        receptors = [
            h for h in hits
            if h.family_key == family.key and h.gene_name in cfg["receptor_genes"]
            and h.superseded_by is None
        ]
        if receptors:
            role = next((g["role"] for g in family.genes if g["name"] == cfg["gene"]), "core_MAT")
            jobs.append((family, cfg, role, receptors))
    if not jobs:
        return []
    contigs = _read_contigs(genome_fasta, {h.contig for _, _, _, rs in jobs for h in rs})
    found: list[SearchHit] = []
    for family, cfg, role, receptors in jobs:
        seen: set[tuple[str, int, int, str]] = set()
        for rec in receptors:
            seq = contigs.get(rec.contig)
            if not seq:
                continue
            lo = max(1, min(rec.start, rec.end) - cfg["window_bp"])
            hi = min(len(seq), max(rec.start, rec.end) + cfg["window_bp"])
            for orf in scan_caax_orfs(seq[lo - 1:hi], cfg["motif"], cfg["min_codons"],
                                      cfg["max_codons"], genetic_code):
                start, end = orf.start + lo - 1, orf.end + lo - 1
                key = (rec.contig, start, end, orf.strand)
                if key in seen:
                    continue
                seen.add(key)
                found.append(SearchHit(
                    family.key, cfg["gene"], role, rec.contig, start, end, orf.strand,
                    0.0, f"{CAAX_METHOD}:{orf.motif}", CAAX_METHOD,
                    align_length_aa=orf.length_aa,
                ))
    return found


def admitted_only_through_scan(
    distinct_genes: set[str], scan_names: set[str], homology_modelled: int, modelled_bar: int,
) -> bool:
    """Does a call reach the admission bar only by counting scan precursors?

    True when it holds a scan precursor and, without it, either fewer than two
    distinct genes remain (the evidence floor) or fewer homology-modelled genes
    than the modelled-gene bar requires. Read by the CAAX unverified label
    (curator's ruling 2026-09-27; `verification.label_caax_unverified`).
    """
    if not scan_names:
        return False
    return len(set(distinct_genes) - set(scan_names)) < 2 or homology_modelled < modelled_bar
