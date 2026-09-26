#!/usr/bin/env python3
"""Harvest assembly-mapped MAT gene coordinates from NCBI Gene, for the
detection coordinate benchmark.

WHY THIS IS SEPARATE FROM THE CURATED RECORDS. A curated record is the QUERY
the pipeline searches with, and most of them are classic locus-specific GenBank
deposits (M33876.1, AF318048.1) whose coordinates are positions within that
deposit, not within any assembly. Scoring detection against a record's own
coordinates would also be circular. What a coordinate benchmark needs is an
INDEPENDENT statement of where the locus sits on a real assembly, and NCBI's
own annotation of that assembly is exactly that.

WHAT IT CAN AND CANNOT MEASURE. This gives "did the pipeline put the locus in
the right place on this assembly". It is NOT a held-out recall estimate: for
most of these species the curated reference is the same species, sometimes the
same strain, so a hit here says nothing about a novel genome. Held-out recall
needs `run_pipeline` to accept a holdout-filtered reference set, which it does
not yet support -- see `detect/benchmark.py`, which keeps that distinction
carefully and reports sensitivity=None rather than pretending.

COORDINATES. NCBI esummary `genomicinfo` reports chrstart/chrstop 0-BASED and
orientation-carrying: chrstart > chrstop means the minus strand. This script
converts to 1-based inclusive start<=end with an explicit strand, matching the
convention the curated records use. Verified against the one record whose
assembly coordinates the curator had already derived by hand: matA-1 came back
1860691..1862193 (0-based) against the curator's 1860692..1862194 (1-based).

Usage:
    python scripts/harvest_coordinate_truth.py --db-root db \
        --out testset/Pezizomycotina/coordinate_truth.yaml
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import yaml

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
#: NCBI asks for <=3 requests/second without an API key. Be a good citizen.
MIN_INTERVAL_S = 0.34
_last = [0.0]


def _get(endpoint: str, **params) -> dict:
    wait = MIN_INTERVAL_S - (time.monotonic() - _last[0])
    if wait > 0:
        time.sleep(wait)
    params.setdefault("retmode", "json")
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=60) as fh:
        body = fh.read().decode()
    _last[0] = time.monotonic()
    return json.loads(body)


def gene_ids(locus_tag: str, organism: str | None) -> list[str]:
    """Gene IDs for a locus tag. The organism term disambiguates a tag that
    several assemblies of one species reuse."""
    term = f"{locus_tag}[Gene Name] OR {locus_tag}[All Fields]"
    if organism:
        term = f"({term}) AND \"{organism}\"[Organism]"
    res = _get("esearch.fcgi", db="gene", term=term)
    return res.get("esearchresult", {}).get("idlist", []) or []


def _names_this_gene(g: dict, locus_tag: str) -> bool:
    """Does this NCBI Gene record actually BELONG to `locus_tag`?

    esearch's `[All Fields]` matches substrings, so a tag is also returned for
    genes whose own name merely contains it. Measured: querying the
    S. cerevisiae tag `YCR097W` returns both HMRA1 (correct) and `YCR097W-A`,
    a different gene 111 bp downstream, which then appeared as a second,
    spurious placement. Require the tag to be the gene's own name or one of
    its aliases, compared exactly.
    """
    fields = [g.get("name") or ""]
    fields += (g.get("otheraliases") or "").split(",")
    fields += (g.get("otherdesignations") or "").split("|")
    return any(f.strip() == locus_tag for f in fields)


def genomic_placements(gene_id: str, locus_tag: str | None = None) -> list[dict]:
    """Every assembly placement NCBI records for this gene, 1-based inclusive."""
    res = _get("esummary.fcgi", db="gene", id=gene_id)
    result = res.get("result", {})
    out = []
    for uid in result.get("uids", []):
        g = result[uid]
        if locus_tag and not _names_this_gene(g, locus_tag):
            continue
        for gi in g.get("genomicinfo", []) or []:
            a, b = int(gi["chrstart"]), int(gi["chrstop"])
            # 0-based, and chrstart > chrstop encodes the minus strand.
            strand = "+" if a <= b else "-"
            start, end = (a, b) if a <= b else (b, a)
            out.append({
                "gene_id": uid,
                "ncbi_gene_name": g.get("name"),
                "sequence_accession": gi.get("chraccver"),
                "start": start + 1,
                "end": end + 1,
                "strand": strand,
                "exon_count": gi.get("exoncount"),
            })
    return out



def _efetch_text(db: str, uid: str, rettype: str) -> str:
    wait = MIN_INTERVAL_S - (time.monotonic() - _last[0])
    if wait > 0:
        time.sleep(wait)
    url = (f"{EUTILS}/efetch.fcgi?"
           + urllib.parse.urlencode({"db": db, "id": uid, "rettype": rettype, "retmode": "text"}))
    with urllib.request.urlopen(url, timeout=120) as fh:
        body = fh.read().decode(errors="replace")
    _last[0] = time.monotonic()
    return body


#: `/coded_by="complement(join(DS499596.1:2521591..2521971,DS499596.1:2522..."`
_CODED_BY = re.compile(r'/coded_by="([^"]+)"', re.S)
_SPAN = re.compile(r"([A-Z]{1,2}[_A-Z0-9]*\.\d+):(\d+)\.\.(\d+)")


def placement_from_protein(protein_acc: str) -> dict | None:
    """Genomic placement read off a protein's GenPept `/coded_by`.

    Works for GenBank/WGS proteins (EDP..., EAS...), whose `coded_by` names a
    genomic contig directly. RefSeq proteins (XP_...) point at an mRNA instead,
    which carries no genomic coordinates -- those fall through to the feature
    table route.
    """
    gp = _efetch_text("protein", protein_acc, "gp")
    m = _CODED_BY.search(gp)
    if not m:
        return None
    expr = " ".join(m.group(1).split())
    spans = _SPAN.findall(expr)
    if not spans:
        return None
    accs = {a for a, _, _ in spans}
    if len(accs) != 1:
        return None
    acc = accs.pop()
    if acc.startswith(("XM_", "NM_")):          # an mRNA, not the genome
        return None
    lo = min(int(b) for _, b, _ in spans)
    hi = max(int(c) for _, _, c in spans)
    return {
        "sequence_accession": acc,
        "start": lo,
        "end": hi,
        "strand": "-" if expr.startswith("complement") else "+",
        "via": "protein /coded_by",
    }


_FT_CACHE: dict[str, list[str]] = {}


def placement_from_feature_table(genomic_acc: str, locus_tag: str) -> dict | None:
    """Genomic placement read off a nuccore FEATURE TABLE `gene` feature.

    The route for RefSeq genes whose NCBI Gene record has an empty
    `genomicinfo` -- which happens when an assembly's annotation has been
    retired, as it has for the Coccidioides immitis RS build. One fetch per
    accession, cached, because a scaffold's table runs to tens of thousands
    of lines.
    """
    lines = _FT_CACHE.get(genomic_acc)
    if lines is None:
        lines = _efetch_text("nuccore", genomic_acc, "ft").splitlines()
        _FT_CACHE[genomic_acc] = lines
    for i, line in enumerate(lines):
        if "locus_tag" not in line or locus_tag not in line:
            continue
        # Walk back to the owning `gene` feature line: "<start>\t<end>\tgene".
        for j in range(i, max(-1, i - 40), -1):
            parts = lines[j].split("\t")
            if len(parts) >= 3 and parts[2] == "gene":
                a, b = parts[0].lstrip("<>"), parts[1].lstrip("<>")
                if not (a.isdigit() and b.isdigit()):
                    break
                a, b = int(a), int(b)
                return {
                    "sequence_accession": genomic_acc,
                    "start": min(a, b),
                    "end": max(a, b),
                    "strand": "+" if a <= b else "-",
                    "via": "nuccore feature table",
                }
        break
    return None


def genomic_accessions_for_gene(gene_id: str) -> list[str]:
    """Genomic nuccore accessions linked to a gene, for the feature-table route."""
    # NOT linkname="gene_nuccore_pos": that link is absent for genes whose
    # assembly annotation has been retired, which is exactly the case this
    # route exists to serve. Take the general link and keep the genomic
    # accessions, skipping transcripts.
    res = _get("elink.fcgi", dbfrom="gene", db="nuccore", id=gene_id)
    ids: list = []
    for ls in res.get("linksets", []):
        for db in ls.get("linksetdbs", []):
            if db.get("linkname") in (None, "gene_nuccore", "gene_nuccore_pos"):
                ids.extend(db.get("links", []) or [])
    out = []
    for uid in ids[:4]:
        sm = _get("esummary.fcgi", db="nuccore", id=str(uid))
        r = sm.get("result", {})
        for u in r.get("uids", []):
            acc = r[u].get("accessionversion") or ""
            if acc and not acc.startswith(("XM_", "NM_", "XR_", "NR_")):
                out.append(acc)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db-root", type=Path, default=Path("db"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--phylum", default="Ascomycota")
    args = ap.parse_args()

    records = []
    for meta in sorted((args.db_root / args.phylum).glob("*/*/metadata.yaml")):
        d = yaml.safe_load(meta.read_text())
        genes = [g for g in (d.get("genes") or []) if g.get("locus_tag")]
        if not genes:
            continue
        species = d["organism"]["species"]
        entry = {
            "record_id": d["record_id"],
            "order": meta.parent.parent.name,
            "species": species,
            "strain": (d["organism"].get("strain") or {}).get("name"),
            "idiomorphs": d["mating_type"]["idiomorphs"],
            "genes": [],
            "unresolved": [],
        }
        for g in genes:
            tag = g["locus_tag"]
            try:
                ids = gene_ids(tag, species)
            except Exception as exc:                      # network/NCBI hiccup
                entry["unresolved"].append({"locus_tag": tag, "error": str(exc)})
                continue
            places = []
            for gid in ids[:3]:
                try:
                    places.extend(genomic_placements(gid, tag))
                except Exception as exc:
                    entry["unresolved"].append({"locus_tag": tag, "error": str(exc)})
            if not places:
                # Fallbacks, in order of directness. NCBI Gene carries an empty
                # `genomicinfo` whenever an assembly's annotation has been
                # retired (Coccidioides immitis RS) and has no Gene record at
                # all for a GenBank-only assembly (Aspergillus fumigatus A1163),
                # so the locus tag alone is not enough for either.
                pacc = (g.get("protein_accession") or "").split(":")[-1]
                got = None
                if pacc:
                    try:
                        got = placement_from_protein(pacc)
                    except Exception as exc:
                        entry["unresolved"].append({"locus_tag": tag, "error": f"coded_by: {exc}"})
                if got is None:
                    for gid in ids[:1]:
                        try:
                            for acc in genomic_accessions_for_gene(gid):
                                got = placement_from_feature_table(acc, tag)
                                if got:
                                    break
                        except Exception as exc:
                            entry["unresolved"].append(
                                {"locus_tag": tag, "error": f"feature table: {exc}"})
                        if got:
                            break
                if got is None:
                    entry["unresolved"].append({"locus_tag": tag, "error": "no genomic placement"})
                    continue
                got.setdefault("gene_id", ids[0] if ids else None)
                places = [got]
            entry["genes"].append({
                "record_gene_name": g.get("name"),
                "role": g.get("role"),
                "locus_tag": tag,
                "placements": places,
            })
            print(f"  {d['record_id'][:34]:36}{tag:18}"
                  f"{places[0]['sequence_accession']} {places[0]['start']}..{places[0]['end']}"
                  f"{places[0]['strand']}", file=sys.stderr)
        records.append(entry)

    doc = {
        "provenance": {
            "source": "NCBI Gene esummary genomicinfo, via E-utilities",
            "harvested_date": date.today().isoformat(),
            "coordinate_convention": "1-based inclusive, start <= end, strand explicit",
            "note": (
                "Independent of the curated records' own coordinates, which are "
                "positions within locus-specific GenBank deposits. Measures "
                "placement on an assembly, NOT held-out recall -- the curated "
                "reference is often the same species or strain."
            ),
        },
        "records": records,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(yaml.safe_dump(doc, sort_keys=False, width=100))
    n_g = sum(len(r["genes"]) for r in records)
    n_u = sum(len(r["unresolved"]) for r in records)
    print(f"\n{len(records)} records, {n_g} genes placed, {n_u} unresolved -> {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
