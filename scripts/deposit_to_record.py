#!/usr/bin/env python3
"""Deterministic generator: one GenBank MAT1-1-1 / MAT1-2-1 deposit -> one core-only curated record.

Pilot for growing the reference roster with annotated, genus-diverse deposits
(docs/notes/2026-10-21_gene-model-extension-design.md sections 14b/15b; curator
approval 2026-10-22). Given a protein accession (for example a row of
results/2026-10-21_gene_model_extension/heldout/selected.tsv), it

1. fetches the GenPept record, the nucleotide record around the CDS, the protein
   FASTA, the NCBI taxonomy record and (when the GenPept cites one) the PubMed
   summary, with public E-utilities, at most 3 requests per second and WITHOUT an
   e-mail address (``tool=`` only); every raw response is cached on disk under
   ``--cache-dir`` so a re-run is byte-identical and offline;
2. derives every field it can mechanically (accession, exons in transcript order
   with the stop codon included, strand, codon_start, transl_table, partial flags,
   organism, strain/isolate, culture collection, taxid, lineage, PubMed id, DOI,
   family by order.yml scope, idiomorph from the gene name);
3. re-checks the result independently (the joined exons translate to the NCBI
   protein exactly; ATG..stop or a recorded partial flag; canonical intron
   dinucleotides; gene span equals the exon bounds unless the deposit's own gene
   feature is wider; no collision with an existing record);
4. writes ``<out>/<record_id>/metadata.yaml`` (a candidate: ``validation.status:
   needs_review``) for every deposit whose verdict is not NEEDS_HUMAN, and a
   report (TSV + JSON) with one verdict per deposit:

   * CONFIRMED             every check passed;
   * CONFIRMED_WITH_CAVEAT every blocking check passed, some informative flag is
                           set (genome-project annotation, third-party annotation,
                           UTR in the gene feature, partial CDS, GC-AG intron, ...);
   * NEEDS_HUMAN           a field could not be derived mechanically or a check
                           failed; no record is written, the reasons are listed.

Two values are NOT derivable from a deposit and are taken from the command line
as curator choices, never guessed: ``--mating-system`` (a single-gene deposit
says nothing about the breeding system) and ``--tier-policy`` (which evidence
tier a depositor's annotation earns). Both are written into the record's
``curation.notes`` as curator-supplied. Nothing in this script calls a model or
uses judgement: the same inputs and cache give the same bytes.

Run (from outside the repository, never with ``uv run`` inside it)::

    uv run --no-project --with biopython --with pyyaml python -I \
        scripts/deposit_to_record.py --repo <worktree> --protein EJP63339.1 \
        --mating-system heterothallic --tier-policy cited2-uncited3 \
        --date 2026-10-22 --cache-dir <cache> --out <dir> --report <prefix>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "MATPredict-deposit-to-record"
MIN_INTERVAL_S = 0.34  # <= 3 requests per second (NCBI's keyless limit)
WINDOW_BP = 20000  # nucleotide window fetched on each side of the CDS
SHORT_INTRON_BP = 35
GENOME_SCALE_BP = 200000  # a nucleotide record this long is a genome contig, not a locus deposit
GENERATOR_VERSION = "1"
LINEAGE_RANKS = [("k", "kingdom"), ("p", "phylum"), ("subphylum", "subphylum"), ("c", "class"),
                 ("o", "order"), ("f", "family"), ("g", "genus"), ("s", "species")]
ROSTER_GENES = {"MAT1-1-1": "MAT1-1", "MAT1-2-1": "MAT1-2"}

CONFIRMED = "CONFIRMED"
CAVEAT = "CONFIRMED_WITH_CAVEAT"
NEEDS_HUMAN = "NEEDS_HUMAN"

TIER_POLICIES = {
    # name -> (tier when the GenPept cites a PubMed article, tier when it cites none)
    "cited1-uncited3": (1, 3),
    "cited2-uncited3": (2, 3),
}


# --------------------------------------------------------------------------- fetching

class Fetcher:
    """GET with an on-disk cache keyed by URL, a 3/s throttle and no e-mail parameter."""

    def __init__(self, cache_dir: Path, offline: bool = False):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self._last = 0.0
        self.network_calls = 0

    def get(self, path: str, **params) -> str:
        params = {**params, "tool": TOOL}
        assert "email" not in params
        url = f"{EUTILS}/{path}?{urllib.parse.urlencode(sorted(params.items()))}"
        key = hashlib.sha256(url.encode()).hexdigest()[:32]
        cached = self.cache_dir / f"{key}.txt"
        if cached.exists():
            return cached.read_text()
        if self.offline:
            raise RuntimeError(f"offline and not cached: {url}")
        last_error = None
        for attempt in range(6):
            wait = MIN_INTERVAL_S - (time.monotonic() - self._last)
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()
            try:
                self.network_calls += 1
                with urllib.request.urlopen(url, timeout=120) as response:
                    body = response.read().decode()
                if body.lstrip().startswith(("Error", '{"error"')) or "<ERROR>" in body[:500]:
                    raise RuntimeError(f"NCBI error body: {body[:200]}")
                cached.write_text(body)
                return body
            except Exception as exc:  # noqa: BLE001 -- retried, then re-raised
                last_error = exc
                time.sleep(min(30, 2 * 2 ** attempt))
        raise RuntimeError(f"fetch failed for {url}: {last_error}")


# --------------------------------------------------------------------------- results

@dataclass
class Check:
    name: str
    status: str  # pass | caveat | needs_human
    detail: str = ""


@dataclass
class Result:
    protein: str
    checks: list[Check] = field(default_factory=list)
    facts: dict = field(default_factory=dict)
    record: dict | None = None
    record_dir: str | None = None

    def add(self, name: str, status: str, detail: str = "") -> None:
        self.checks.append(Check(name, status, detail))

    @property
    def verdict(self) -> str:
        statuses = {c.status for c in self.checks}
        if "needs_human" in statuses:
            return NEEDS_HUMAN
        if "caveat" in statuses:
            return CAVEAT
        return CONFIRMED

    def reasons(self, status: str) -> list[str]:
        return [f"{c.name}: {c.detail}" for c in self.checks if c.status == status]


# --------------------------------------------------------------------------- helpers

_REMNANT_WORDS = re.compile(r"truncat|pseudo|remnant|fragment|disrupted", flags=re.I)
_CODED_BY_PART = re.compile(r"([A-Z]{1,6}_?[0-9]+\.[0-9]+):<?([0-9]+)\.\.>?([0-9]+)")


def normalise_gene_name(value: str | None) -> str | None:
    """MAT1-1-1 / MAT1-2-1 from a /gene value spelled with any case or spacing, else None."""
    if not value:
        return None
    v = re.sub(r"\s+", "", value).upper()
    return v if v in ROSTER_GENES else None


def gene_name_from_product(product: str | None) -> str | None:
    """The single roster name appearing as a token in /product, else None (0 or 2 names)."""
    if not product:
        return None
    found = {normalise_gene_name(t) for t in re.findall(r"MAT\s*1-[12]-1", product, flags=re.I)}
    found.discard(None)
    return found.pop() if len(found) == 1 else None


def lineage_string(lineage_ex: list[dict], self_node: dict) -> tuple[str, dict[str, str], list[str]]:
    """k__..;p__..;subphylum__..;c__..;o__..;f__..;g__..;s__.. from NCBI taxonomy ranks.

    The Ascomycota records in db/ all use this eight-field form. `s__` is the
    species-rank node's name (spaces -> `_`), so an infraspecific taxid (strain,
    forma specialis) still gets its species. A rank NCBI does not assign is left
    empty, as `taxonkit reformat` leaves it.
    """
    nodes = lineage_ex + [self_node]
    by_rank = {}
    for node in nodes:
        by_rank.setdefault(node["rank"], node["name"])
    missing = []
    parts = []
    for prefix, rank in LINEAGE_RANKS:
        name = by_rank.get(rank, "")
        if not name:
            missing.append(rank)
        parts.append(f"{prefix}__{name.replace(' ', '_')}")
    return ";".join(parts), by_rank, missing


def slugify_strain(name: str) -> str:
    """Same rule as MATPredict.db.identifiers.slugify_strain (kept identical, tested)."""
    slug = name.strip().lower().replace("/", "-")
    slug = re.sub(r"\s+", "-", slug)
    slug = re.sub(r"[^a-z0-9-]", "", slug)
    return re.sub(r"-{2,}", "-", slug).strip("-")


def stop_codons(table: int) -> set[str]:
    return set(CodonTable.unambiguous_dna_by_id[table].stop_codons)


# --------------------------------------------------------------------------- the db

def load_db_index(db_root: Path) -> dict:
    """Accepted and candidate records: protein accessions, nucleotide accessions, ids, taxids, species."""
    index = {"proteins": {}, "nucleotides": {}, "record_ids": set(), "taxids": {}, "species": {}}
    for meta in sorted(db_root.glob("**/metadata.yaml")):
        doc = yaml.safe_load(meta.read_text()) or {}
        rid = doc.get("record_id")
        if not rid:
            continue
        index["record_ids"].add(rid)
        index["taxids"].setdefault((doc.get("taxonomy") or {}).get("taxid"), []).append(rid)
        index["species"].setdefault((doc.get("organism") or {}).get("species"), []).append(rid)
        for gene in doc.get("genes") or []:
            acc = gene.get("protein_accession")
            if acc:
                index["proteins"].setdefault(acc.split(":", 1)[-1].split(".")[0], []).append(rid)
        for seg in ((doc.get("locus") or {}).get("core") or {}).get("segments") or []:
            acc = (seg.get("sequence_source") or {}).get("accession")
            if acc:
                index["nucleotides"].setdefault(acc.split(".")[0], []).append(rid)
    return index


def route_family(db_root: Path, phylum: str, lineage_taxids: set[int], gene: str) -> list[tuple[str, dict]]:
    """Every order.yml family of `phylum` whose scope covers the lineage and that declares `gene`."""
    order_path = db_root / phylum / "order.yml"
    if not order_path.exists():
        return []
    doc = yaml.safe_load(order_path.read_text())
    hits = []
    for locus in doc["loci"]:
        if not set(locus.get("taxonomic_scope") or []) & lineage_taxids:
            continue
        names = {n for g in locus.get("genes", []) for n in [g["name"], *(g.get("aliases") or [])]}
        if gene in names:
            hits.append((locus["locus_name"], locus))
    return hits


# --------------------------------------------------------------------------- core

def parse_taxonomy(xml_text: str) -> dict:
    root = ET.fromstring(xml_text)
    taxon = root.find("Taxon")
    if taxon is None:
        raise ValueError("no Taxon in taxonomy efetch")
    lineage_ex = [{"taxid": int(t.findtext("TaxId")), "name": t.findtext("ScientificName"),
                   "rank": t.findtext("Rank")} for t in taxon.findall("LineageEx/Taxon")]
    return {
        "taxid": int(taxon.findtext("TaxId")),
        "name": taxon.findtext("ScientificName"),
        "rank": taxon.findtext("Rank"),
        "lineage_ex": lineage_ex,
        "genetic_code": int(taxon.findtext("GeneticCode/GCId") or 0) or None,
    }


def build(protein: str, fetcher: Fetcher, db_root: Path, db_index: dict, *, date: str,
          mating_system: str, tier_policy: str) -> Result:
    res = Result(protein=protein)
    f = res.facts

    # ---- GenPept
    gp_text = fetcher.get("efetch.fcgi", db="protein", id=protein, rettype="gp", retmode="text")
    gp = SeqIO.read(StringIO(gp_text), "genbank")
    f["protein_version"] = gp.id
    if gp.id != protein:
        res.add("protein_version", "needs_human", f"requested {protein}, NCBI returned {gp.id}")
    gp_source = next((x for x in gp.features if x.type == "source"), None)
    gp_cds = [x for x in gp.features if x.type == "CDS"]
    if len(gp_cds) != 1 or "coded_by" not in gp_cds[0].qualifiers:
        res.add("coded_by", "needs_human", "GenPept has no single CDS feature with /coded_by")
        return res
    coded_by = gp_cds[0].qualifiers["coded_by"][0]
    f["coded_by"] = coded_by
    parts = _CODED_BY_PART.findall(coded_by.replace(" ", ""))
    accessions = {p[0] for p in parts}
    if "order(" in coded_by or len(accessions) != 1 or not parts:
        res.add("coded_by", "needs_human", f"coded_by is not one accession / join: {coded_by}")
        return res
    nuc = accessions.pop()
    f["nucleotide"] = nuc
    lo = min(int(p[1]) for p in parts)
    hi = max(int(p[2]) for p in parts)
    res.add("coded_by", "pass", coded_by)

    gp_gene = (gp_cds[0].qualifiers.get("gene") or [None])[0]
    gp_prot = next((x for x in gp.features if x.type in ("Protein", "Prot")), None)
    product = (gp_prot.qualifiers.get("product") or [None])[0] if gp_prot is not None else None
    if product is None:
        product = (gp_cds[0].qualifiers.get("product") or [None])[0]
    f["gene_qualifier"], f["product"] = gp_gene, product

    # ---- citations (GenPept REFERENCE blocks with a PUBMED line)
    pmids = []
    for ref in gp.annotations.get("references", []):
        if ref.pubmed_id and ref.pubmed_id not in pmids:
            pmids.append(ref.pubmed_id)
    f["pmids"] = pmids
    citations = []
    for pmid in pmids:
        summ = json.loads(fetcher.get("esummary.fcgi", db="pubmed", id=pmid, retmode="json"))
        doc = summ.get("result", {}).get(pmid, {})
        doi = next((a["value"] for a in doc.get("articleids", []) if a.get("idtype") == "doi"), None)
        citations.append({"pmid": pmid, "doi": doi,
                          "_ref": f"{(doc.get('sortfirstauthor') or '').strip()} {doc.get('pubdate', '')[:4]}, "
                                  f"{doc.get('source', '')}"})
    if len(pmids) > 1:
        res.add("citations", "caveat", f"GenPept cites {len(pmids)} PubMed articles: {pmids}")
    elif not pmids:
        res.add("citations", "caveat", "GenPept cites no PubMed article (unpublished or not linked)")
    else:
        res.add("citations", "pass", pmids[0])
    for c in citations:
        if c["doi"] is None:
            res.add("doi", "caveat", f"PubMed {c['pmid']} has no DOI in its summary")

    # ---- protein sequence (independent of the GenPept/nucleotide translation)
    fasta = fetcher.get("efetch.fcgi", db="protein", id=protein, rettype="fasta", retmode="text")
    ncbi_protein = str(SeqIO.read(StringIO(fasta), "fasta-blast").seq)
    f["protein_aa"] = len(ncbi_protein)

    # ---- nucleotide: summary (length, version) and a window around the CDS
    nsum = json.loads(fetcher.get("esummary.fcgi", db="nuccore", id=nuc, retmode="json"))
    uid = (nsum.get("result", {}).get("uids") or [None])[0]
    ndoc = nsum.get("result", {}).get(uid, {}) if uid else {}
    slen = int(ndoc.get("slen") or 0)
    f["nucleotide_length"] = slen
    f["nucleotide_title"] = ndoc.get("title")
    if ndoc.get("accessionversion") != nuc:
        res.add("nucleotide_version", "needs_human",
                f"coded_by names {nuc}; NCBI's current version is {ndoc.get('accessionversion')}")
    w_start = max(1, lo - WINDOW_BP)
    w_end = min(slen, hi + WINDOW_BP) if slen else hi + WINDOW_BP
    gb_text = fetcher.get("efetch.fcgi", db="nuccore", id=nuc, rettype="gbwithparts", retmode="text",
                          seq_start=str(w_start), seq_stop=str(w_end))
    gb = SeqIO.read(StringIO(gb_text), "genbank")
    offset = w_start - 1
    seq = str(gb.seq).upper()
    keywords = [k for k in gb.annotations.get("keywords", []) if k]
    f["keywords"] = keywords
    f["molecule_type"] = gb.annotations.get("molecule_type")
    f["nucleotide_window"] = [w_start, w_end]

    cds_list = [x for x in gb.features if x.type == "CDS"
                and (x.qualifiers.get("protein_id") or [None])[0] == protein]
    if len(cds_list) != 1:
        res.add("cds_found", "needs_human", f"{len(cds_list)} CDS with protein_id {protein} in {nuc}")
        return res
    cds = cds_list[0]
    res.add("cds_found", "pass", f"{nuc}:{lo}..{hi}")
    loc = cds.location
    strand = "-" if loc.strand == -1 else "+"
    if any(p.strand != loc.strand for p in loc.parts):
        res.add("strand", "needs_human", "CDS parts on mixed strands")
        return res
    exons = [(int(p.start) + 1 + offset, int(p.end) + offset) for p in loc.parts]
    partial5 = isinstance(loc.parts[0].start if strand == "+" else loc.parts[0].end,
                          (BeforePosition, AfterPosition))
    partial3 = isinstance(loc.parts[-1].end if strand == "+" else loc.parts[-1].start,
                          (BeforePosition, AfterPosition))
    codon_start = int((cds.qualifiers.get("codon_start") or ["1"])[0])
    transl_table = int((cds.qualifiers.get("transl_table") or ["1"])[0])
    f.update(strand=strand, exons=exons, partial5=partial5, partial3=partial3,
             codon_start=codon_start, transl_table=transl_table,
             locus_tag=(cds.qualifiers.get("locus_tag") or [None])[0])

    # exon order / overlap
    asc = sorted(exons)
    ordered = exons == (asc if strand == "+" else asc[::-1])
    overlap = any(asc[i][1] >= asc[i + 1][0] for i in range(len(asc) - 1))
    if not ordered or overlap:
        res.add("exon_order", "needs_human", f"exons not in transcript order or overlapping: {exons}")
        return res
    res.add("exon_order", "pass", f"{len(exons)} exons, strand {strand}")
    if f["coded_by"] and (min(e[0] for e in exons), max(e[1] for e in exons)) != (lo, hi):
        res.add("coded_by_agrees", "needs_human",
                f"nucleotide CDS bounds {min(e[0] for e in exons)}..{max(e[1] for e in exons)} != coded_by {lo}..{hi}")

    # ---- independent translation from the joined exons
    def sub(s: int, e: int) -> str:
        return seq[s - 1 - offset:e - offset]

    pieces = [sub(s, e) for s, e in exons]
    if strand == "-":
        pieces = [str(Seq(p).reverse_complement()) for p in pieces]
    cds_nt = "".join(pieces)
    framed = cds_nt[codon_start - 1:]
    trimmed = framed[:len(framed) - len(framed) % 3]
    aa = str(Seq(trimmed).translate(table=transl_table))
    stops = stop_codons(transl_table)
    has_stop = bool(trimmed) and trimmed[-3:] in stops
    body = aa[:-1] if has_stop else aa
    if "*" in body:
        res.add("translation", "needs_human", f"internal stop at aa {body.index('*') + 1}")
    elif body != ncbi_protein:
        res.add("translation", "needs_human",
                f"joined exons translate to {len(body)} aa, NCBI protein is {len(ncbi_protein)} aa; "
                f"first difference at aa {next((i + 1 for i, (a, b) in enumerate(zip(body, ncbi_protein)) if a != b), min(len(body), len(ncbi_protein)) + 1)}")
    else:
        res.add("translation", "pass", f"{len(body)} aa identical to {protein}")
    q_translation = (cds.qualifiers.get("translation") or [None])[0]
    if q_translation is not None and q_translation != ncbi_protein:
        res.add("translation_qualifier", "needs_human", "/translation differs from the NCBI protein record")

    # ---- start / stop / frame
    if partial5:
        res.add("start", "caveat", f"5' partial CDS (<), codon_start {codon_start}")
    elif codon_start != 1:
        res.add("start", "needs_human", f"5' complete CDS with codon_start {codon_start}")
    elif cds_nt[:3] != "ATG":
        res.add("start", "needs_human", f"5' complete CDS starts with {cds_nt[:3]}, not ATG")
    else:
        res.add("start", "pass", "ATG")
    if partial3:
        res.add("stop", "caveat", "3' partial CDS (>)")
    elif not has_stop:
        res.add("stop", "needs_human", f"3' complete CDS does not end in a stop codon ({trimmed[-3:]})")
    elif len(framed) % 3:
        res.add("stop", "needs_human", f"3' complete CDS length {len(framed)} is not a multiple of 3")
    else:
        res.add("stop", "pass", f"{trimmed[-3:]} included in the last exon")

    # ---- introns
    introns = []
    for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
        istart, iend = (e1 + 1, s2 - 1) if strand == "+" else (e2 + 1, s1 - 1)
        intron = sub(istart, iend)
        if strand == "-":
            intron = str(Seq(intron).reverse_complement())
        introns.append((istart, iend, intron[:2], intron[-2:], len(intron)))
    f["introns"] = [list(i) for i in introns]
    for istart, iend, donor, acceptor, length in introns:
        tag = f"intron {istart}..{iend} ({length} bp) {donor}..{acceptor}"
        if (donor, acceptor) == ("GT", "AG"):
            status = "pass"
        elif (donor, acceptor) == ("GC", "AG"):
            status = "caveat"
        else:
            status = "needs_human"
        if status == "pass" and length < SHORT_INTRON_BP:
            status = "caveat"
            tag += f" shorter than {SHORT_INTRON_BP} bp"
        res.add("intron", status, tag)

    # ---- gene feature span versus the CDS span
    cds_lo, cds_hi = min(e[0] for e in exons), max(e[1] for e in exons)
    gene_name_q = (cds.qualifiers.get("gene") or [None])[0]
    genes = [g for g in gb.features if g.type == "gene"
             and int(g.location.start) + 1 + offset <= cds_lo and int(g.location.end) + offset >= cds_hi
             and ((g.qualifiers.get("gene") or [None])[0] == gene_name_q
                  or (g.qualifiers.get("locus_tag") or [None])[0] == f["locus_tag"])]
    if not genes:
        res.add("gene_span", "caveat", "no gene feature spans the CDS; record uses the CDS span")
        f["gene_feature"] = None
    else:
        g = genes[0]
        g_lo, g_hi = int(g.location.start) + 1 + offset, int(g.location.end) + offset
        f["gene_feature"] = [g_lo, g_hi]
        if (g_lo, g_hi) == (cds_lo, cds_hi):
            res.add("gene_span", "pass", f"gene feature = exon bounds {cds_lo}..{cds_hi}")
        else:
            five = (cds_lo - g_lo) if strand == "+" else (g_hi - cds_hi)
            three = (g_hi - cds_hi) if strand == "+" else (cds_lo - g_lo)
            res.add("gene_span", "caveat",
                    f"deposit gene feature {g_lo}..{g_hi} is wider than the CDS {cds_lo}..{cds_hi} "
                    f"(5' {five} bp, 3' {three} bp; UTRs or a gene-model span); record uses the CDS span")

    # ---- gene name -> roster gene -> idiomorph
    name = normalise_gene_name(gene_name_q) or normalise_gene_name(gp_gene)
    if name is None:
        from_product = gene_name_from_product(product)
        if from_product:
            name = from_product
            res.add("gene_name", "caveat", f"no /gene MAT1-1-1|MAT1-2-1; named from /product '{product}'")
        else:
            res.add("gene_name", "needs_human", f"/gene {gene_name_q!r}, /product {product!r}: not a roster name")
            return res
    else:
        if (gene_name_q and normalise_gene_name(gene_name_q) is None) or \
                (gene_name_from_product(product) not in (None, name)):
            res.add("gene_name", "needs_human", f"/gene {gene_name_q!r} and /product {product!r} disagree")
            return res
        res.add("gene_name", "pass", name)
    f["gene"] = name
    idiomorph = ROSTER_GENES[name]

    # the other idiomorph's gene on the same accession -> not a heterothallic single-idiomorph deposit
    own_text = " ".join((cds.qualifiers.get("note") or []) + [product or ""])
    if "pseudo" in cds.qualifiers or "pseudogene" in cds.qualifiers or _REMNANT_WORDS.search(own_text):
        res.add("own_status", "needs_human",
                f"the deposit marks this CDS as pseudo/truncated/remnant: '{own_text.strip()}'")
    other_named = []  # (gene, protein_id, is_remnant)
    for x in gb.features:
        if x.type != "CDS" or x is cds:
            continue
        n = normalise_gene_name((x.qualifiers.get("gene") or [None])[0]) or \
            gene_name_from_product((x.qualifiers.get("product") or [None])[0])
        if n:
            text = " ".join((x.qualifiers.get("note") or []) + (x.qualifiers.get("product") or []))
            remnant = ("pseudo" in x.qualifiers or "pseudogene" in x.qualifiers
                       or bool(_REMNANT_WORDS.search(text)))
            other_named.append((n, (x.qualifiers.get("protein_id") or ["-"])[0], remnant))
    f["other_mat_cds"] = other_named
    intact_other = [p for n, p, rem in other_named if n != name and not rem]
    remnant_other = [p for n, p, rem in other_named if n != name and rem]
    if intact_other:
        res.add("other_idiomorph", "needs_human",
                f"the deposit also carries an intact CDS of the other idiomorph's gene ({intact_other}) "
                f"within {WINDOW_BP} bp: idiomorph label and mating system need a curator")
    if remnant_other:
        res.add("other_idiomorph", "caveat",
                f"the deposit also carries a pseudo/truncated CDS of the other idiomorph's gene "
                f"({remnant_other}); not recorded (a curator may record it as status: remnant)")
    if any(n == name for n, _, _ in other_named):
        res.add("other_idiomorph", "caveat", f"another {name} CDS in the window: {other_named}")

    # ---- annotation-source flags (informative)
    inference = cds.qualifiers.get("inference") or []
    note = " ".join(cds.qualifiers.get("note") or [])
    structured = gb.annotations.get("structured_comment") or {}
    genome_scale = ("WGS" in keywords or "Genome-Annotation-Data" in structured
                    or "Genome-Assembly-Data" in structured or slen >= GENOME_SCALE_BP)
    f["genome_scale_record"] = genome_scale
    if genome_scale:
        res.add("annotation_source", "caveat",
                f"genome-project record (WGS keyword, genome structured comment or >= {GENOME_SCALE_BP} bp; "
                f"{slen} bp): the CDS may be a pipeline gene prediction")
    if any(k.startswith("TPA") or "Third Party" in k for k in keywords):
        res.add("annotation_source", "caveat", f"third-party annotation (TPA) record; keywords {keywords}")
    if inference:
        res.add("annotation_source", "caveat", f"/inference {inference}")
    if re.search(r"predict|putative|hypothetical", f"{note} {product or ''}", flags=re.I):
        res.add("annotation_source", "caveat", f"note/product says predicted/putative: '{note or product}'")
    if f["molecule_type"] not in ("DNA", "genomic DNA"):
        res.add("molecule", "needs_human", f"molecule type {f['molecule_type']}: not genomic DNA")

    # ---- organism, strain, taxid
    src = next((x for x in gb.features if x.type == "source"), None)
    q = src.qualifiers if src is not None else {}
    gq = gp_source.qualifiers if gp_source is not None else {}
    organism = (q.get("organism") or [None])[0]
    f["organism"] = organism
    for key in ("organism", "strain", "isolate"):
        if (q.get(key) or [None])[0] != (gq.get(key) or [None])[0]:
            res.add("source_agrees", "needs_human",
                    f"/{key} differs between nucleotide ({q.get(key)}) and GenPept ({gq.get(key)})")
    taxon_xrefs = [x.split(":", 1)[1] for x in q.get("db_xref", []) if x.startswith("taxon:")]
    if len(taxon_xrefs) != 1:
        res.add("taxid", "needs_human", f"source has {len(taxon_xrefs)} taxon db_xrefs")
        return res
    taxid = int(taxon_xrefs[0])
    strain = (q.get("strain") or [None])[0]
    isolate = (q.get("isolate") or [None])[0]
    collections = [c.replace(":", " ") for c in q.get("culture_collection", [])]
    f.update(strain_q=strain, isolate_q=isolate, culture_collection=collections)
    if strain:
        strain_name, strain_from = strain, "strain"
    elif isolate:
        strain_name, strain_from = isolate, "isolate"
        res.add("strain", "caveat", f"no /strain; strain name taken from /isolate '{isolate}'")
    else:
        strain_name, strain_from = None, None
        res.add("strain", "caveat", "no /strain or /isolate: strain unknown")
    if strain and isolate and strain != isolate:
        res.add("strain", "caveat", f"/strain '{strain}' and /isolate '{isolate}' both present; /strain used")

    tax = parse_taxonomy(fetcher.get("efetch.fcgi", db="taxonomy", id=str(taxid), retmode="xml"))
    if tax["taxid"] != taxid:
        res.add("taxid", "needs_human", f"taxid {taxid} resolves to {tax['taxid']} (merged)")
        return res
    lineage, by_rank, missing = lineage_string(tax["lineage_ex"], {"taxid": taxid, "name": tax["name"],
                                                                     "rank": tax["rank"]})
    f.update(taxid=taxid, lineage=lineage, taxon_name=tax["name"], taxon_rank=tax["rank"],
             genetic_code=tax["genetic_code"])
    if tax["name"] != organism:
        res.add("taxid", "caveat", f"taxid {taxid} is named '{tax['name']}', deposit says '{organism}'")
    if tax["rank"] != "species":
        res.add("taxid", "caveat", f"taxid {taxid} has rank '{tax['rank']}' (s__ is the species-rank ancestor)")
    if "order" in missing or "phylum" in missing:
        res.add("lineage", "needs_human", f"NCBI lineage lacks rank(s) {missing}: no order directory")
        return res
    if missing:
        res.add("lineage", "caveat", f"NCBI lineage lacks rank(s) {missing}: left empty")
    else:
        res.add("lineage", "pass", lineage)
    if tax["genetic_code"] and tax["genetic_code"] != transl_table:
        res.add("transl_table", "needs_human",
                f"CDS uses table {transl_table}, the taxon's nuclear code is {tax['genetic_code']}")

    phylum, order = by_rank["phylum"], by_rank["order"]
    lineage_taxids = {t["taxid"] for t in tax["lineage_ex"]} | {taxid}
    families = route_family(db_root, phylum, lineage_taxids, name)
    if len(families) != 1:
        res.add("family", "needs_human",
                f"{len(families)} {phylum} order.yml families cover this lineage and declare {name}: "
                f"{[n for n, _ in families]}")
        return res
    locus_name = families[0][0]
    idiomorph_values = families[0][1].get("idiomorph_values") or []
    if idiomorph not in idiomorph_values:
        res.add("family", "needs_human", f"{idiomorph} not in {locus_name} idiomorph_values {idiomorph_values}")
        return res
    res.add("family", "pass", f"{phylum}:{locus_name} (order directory {order})")

    # ---- identifiers and collisions
    if strain_name:
        slug = slugify_strain(strain_name)
        if not slug:
            res.add("record_id", "needs_human", f"strain '{strain_name}' slugs to nothing")
            return res
    else:
        n = 1
        existing = db_index["record_ids"]
        while any(r.startswith(f"{taxid}_unknown-{n}_") for r in existing):
            n += 1
        slug = f"unknown-{n}"
    record_id = f"{taxid}_{slug}_{locus_name}_{idiomorph}"
    f["record_id"] = record_id
    if not re.match(r"^[0-9]+_[a-z0-9-]+_[A-Za-z0-9-]+_[A-Za-z0-9-]+$", record_id):
        res.add("record_id", "needs_human", f"record_id {record_id} does not match the schema pattern")
        return res
    if record_id in db_index["record_ids"]:
        res.add("record_id", "needs_human", f"record_id {record_id} already exists in db/")
        return res
    bare = protein.split(".")[0]
    if bare in db_index["proteins"]:
        res.add("duplicate", "needs_human", f"protein {protein} already curated in {db_index['proteins'][bare]}")
        return res
    if nuc.split(".")[0] in db_index["nucleotides"]:
        res.add("duplicate", "needs_human",
                f"nucleotide {nuc} is already a segment of {db_index['nucleotides'][nuc.split('.')[0]]}")
        return res
    same = sorted(set(db_index["taxids"].get(taxid, [])) | set(db_index["species"].get(organism, [])))
    if same:
        res.add("same_species", "caveat", f"db/ already holds {same} for this taxid/species")
    res.add("record_id", "pass", record_id)

    # ---- the record
    tier_cited, tier_uncited = TIER_POLICIES[tier_policy]
    tier = tier_cited if pmids else tier_uncited
    cits = [{"pmid": c["pmid"], "doi": c["doi"]} for c in citations]
    deposit_desc = f"GenBank {nuc} CDS {coded_by} (protein {protein})"
    caveats = res.reasons("caveat")
    notes = (
        f"Generated {date} by scripts/deposit_to_record.py v{GENERATOR_VERSION} from NCBI E-utilities "
        f"(GenPept {protein}, nuccore {nuc}, taxonomy {taxid}"
        + (", PubMed " + ",".join(pmids) if pmids else "") + "). "
        f"Independent checks: the joined exons translate to {protein} exactly; "
        f"{'ATG start' if not partial5 else '5-partial'}, "
        f"{'stop codon included' if not partial3 else '3-partial'}; "
        f"{len(introns)} intron(s) {', '.join(d + '..' + a for _, _, d, a, _ in introns) or '-'}. "
        f"CURATOR-SUPPLIED, not derived: mating_type.system={mating_system}; evidence tier policy "
        f"{tier_policy} (tier {tier}). "
    )
    notes += f"Strain from /{strain_from}." if strain_from else "Strain unknown in the deposit."
    if caveats:
        notes += " Caveats: " + " | ".join(caveats)
    span_note = ""
    if f.get("gene_feature") and tuple(f["gene_feature"]) != (cds_lo, cds_hi):
        span_note = f" The deposit's gene feature is {f['gene_feature'][0]}..{f['gene_feature'][1]}."
    record = {
        "record_id": record_id,
        "record_version": 1,
        "taxonomy": {"taxid": taxid, "lineage": lineage, "lineage_resolved_date": date},
        "organism": {
            "species": organism,
            "strain": {"name": strain_name, "known": strain_name is not None,
                       "culture_collection_ids": collections, "differs_from_sequenced": False},
        },
        "mating_type": {"locus_name": locus_name, "idiomorphs": [idiomorph], "system": mating_system},
        "locus": {
            "assembly_accession": None,
            "coordinate_provenance": "published_explicit",
            "excluded_from_coordinate_benchmark": False,
            "core": {
                "completeness": "partial",
                "reference_orientation": f"{name}({strand})",
                "definition_note": (
                    f"Core-only record generated from {deposit_desc}; the depositor's annotation, "
                    f"{len(exons)} exon(s), {'5-partial ' if partial5 else ''}{'3-partial ' if partial3 else ''}"
                    f"{len(ncbi_protein)} aa. The segment is the CDS span; idiomorph boundaries and "
                    f"flanking genes are not recorded.{span_note}"
                    + (f" Nucleotide record: {f['nucleotide_title']}." if f.get("nucleotide_title") else "")
                    + (" Cited: " + "; ".join(c["_ref"] for c in citations) + "." if citations else "")
                ),
                "segments": [{
                    "segment_index": 0,
                    "sequence_source": {"type": "insdc_nucleotide", "accession": nuc, "seq_region": nuc},
                    "start": cds_lo, "end": cds_hi,
                    "contig_edge_distance": None, "sequence_checksum": None,
                }],
            },
            "extended_flank": [],
        },
        "genes": [{
            "gene_index": 0, "name": name, "protein_accession": f"ncbi_protein:{protein}",
            "role": "core_MAT", "present": True, "locus_tag": f["locus_tag"], "segment_index": 0,
            "start": cds_lo, "end": cds_hi, "strand": strand, "order_in_locus": 0,
            "exons": [{"start": s, "end": e} for s, e in exons],
            "codon_start": codon_start, "transl_table": transl_table,
        }],
        "evidence": {
            "locus_existence": {
                "tier": tier,
                "experimental_method": f"{name} sequenced and annotated by the depositors ({deposit_desc}).",
                "citations": cits,
            },
            "boundaries": {
                "tier": tier,
                "experimental_method": "CDS exon coordinates are the deposit's own annotation; no idiomorph "
                                       "boundary is recorded.",
                "citations": cits,
            },
            "idiomorph_assignment": {
                "tier": tier,
                "experimental_method": f"idiomorph {idiomorph} from the depositors' gene name {name} "
                                       f"({'alpha box' if name == 'MAT1-1-1' else 'HMG box'} roster gene).",
                "citations": cits,
            },
        },
        "validation": {"status": "needs_review", "rejection_reason": None},
        "curation": {
            "proposed_by": "deposit-to-record-generator",
            "proposal_dedupe_key": f"{pmids[0] if pmids else nuc}|{taxid}|{locus_name}|{idiomorph}",
            "reviewed_by": None,
            "reviewed_date": None,
            "notes": notes,
        },
        "model_provenance": None,
    }
    res.record = record
    f["order"] = order
    f["phylum"] = phylum
    return res


# --------------------------------------------------------------------------- cli

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", required=True, type=Path, help="repository/worktree root (reads db/)")
    ap.add_argument("--protein", action="append", default=[], help="protein accession (repeatable)")
    ap.add_argument("--table", type=Path, help="TSV with a `protein` column (e.g. heldout/selected.tsv)")
    ap.add_argument("--mating-system", required=True,
                    choices=["heterothallic", "homothallic", "pseudohomothallic", "unisexual_switching"],
                    help="CURATOR CHOICE written to every generated record (not derivable from a deposit)")
    ap.add_argument("--tier-policy", required=True, choices=sorted(TIER_POLICIES),
                    help="CURATOR CHOICE: evidence tier for cited / uncited deposits")
    ap.add_argument("--date", required=True, help="YYYY-MM-DD written as lineage_resolved_date and in notes")
    ap.add_argument("--cache-dir", required=True, type=Path)
    ap.add_argument("--offline", action="store_true", help="use the cache only")
    ap.add_argument("--out", type=Path, help="write <out>/<phylum>/<order>/<record_id>/metadata.yaml")
    ap.add_argument("--report", required=True, type=Path, help="report prefix (.tsv and .json are written)")
    args = ap.parse_args(argv)

    proteins = list(args.protein)
    if args.table:
        import csv
        with open(args.table) as fh:
            proteins += [r["protein"] for r in csv.DictReader(fh, delimiter="\t")]
    if not proteins:
        ap.error("give --protein or --table")
    db_root = args.repo / "db"
    db_index = load_db_index(db_root)
    fetcher = Fetcher(args.cache_dir, offline=args.offline)

    results = []
    for protein in proteins:
        try:
            res = build(protein, fetcher, db_root, db_index, date=args.date,
                        mating_system=args.mating_system, tier_policy=args.tier_policy)
        except Exception as exc:  # noqa: BLE001 -- one bad deposit never aborts the batch
            res = Result(protein=protein)
            res.add("exception", "needs_human", f"{type(exc).__name__}: {exc}")
        if res.record is not None and res.verdict != NEEDS_HUMAN and args.out:
            rid = res.record["record_id"]
            d = args.out / res.facts["phylum"] / res.facts["order"] / rid
            d.mkdir(parents=True, exist_ok=True)
            (d / "metadata.yaml").write_text(yaml.safe_dump(res.record, sort_keys=False, width=100))
            res.record_dir = str(d)
            db_index["record_ids"].add(rid)  # a second deposit of the same strain must not collide
        results.append(res)
        print(f"{protein}\t{res.verdict}\t{res.facts.get('record_id', '-')}", file=sys.stderr)

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{args.report}.tsv", "w") as fh:
        fh.write("protein\tverdict\trecord_id\torganism\tgene\tnucleotide\tstrand\texons\taa\t"
                 "partial5\tpartial3\tpmids\tneeds_human\tcaveats\n")
        for r in results:
            f = r.facts
            fh.write("\t".join(str(x) for x in [
                r.protein, r.verdict, f.get("record_id", ""), f.get("organism", ""), f.get("gene", ""),
                f.get("nucleotide", ""), f.get("strand", ""),
                ";".join(f"{s}-{e}" for s, e in f.get("exons", [])), f.get("protein_aa", ""),
                f.get("partial5", ""), f.get("partial3", ""), ",".join(f.get("pmids", [])),
                " | ".join(r.reasons("needs_human")), " | ".join(r.reasons("caveat")),
            ]) + "\n")
    with open(f"{args.report}.json", "w") as fh:
        json.dump([{"protein": r.protein, "verdict": r.verdict, "record_dir": r.record_dir,
                    "facts": r.facts, "checks": [c.__dict__ for c in r.checks]} for r in results],
                  fh, indent=1, default=str)
    print(f"{len(results)} deposits; network calls {fetcher.network_calls}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
