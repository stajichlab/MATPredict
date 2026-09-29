"""Genome-level statement when one family calls two idiomorphs.

Curator's ruling 2026-09-29: TEST FIRST. A genome with a determined Plus call
and a determined Minus call in one family is NOT a confirmed homothallic. The
causes can be homothallism, a hybrid or fusion, a duplication, a mixed culture
or heterokaryon, or an assembly artefact. So this module never asserts a
cause. It reports the arrangement, the evidence, and every possible cause with
the support each one has (possibly none).

Evidence: results/2026-09-29_mucoro_homothallism_literature/NOTE.md.
Idnurm 2011 (Syzygites, two loci, a pseudogenised flank copy) and Schulz et al.
2016 (Zygorhynchus heterogamus one locus, sexM-sexP 5.3 kb; Mycotypha africana
~150 kb apart on one scaffold; duplicated tptA/algL copies 74-76% identical).
Protoplast fusion of (+) and (-) Absidia glauca strains gives homothallic-
looking strains, so a two-idiomorph genome can also be a fusion product.

Report-level only: no call, idiomorph or confidence is changed. On for a
family through the roster field `two_idiomorphs_report` (Mucoromycota:MAT).

The rules (each is a stated threshold, not a fitted one):
* arrangement: `same_locus` when a call is `homothallic_candidate` or the
  Plus and Minus spans on one contig are within DISTANT_BP; `same_contig_distant`
  beyond it; `unlinked` on different contigs.
* duplication: more than one determined call for one idiomorph.
* mixed_culture_or_heterokaryon: the two contigs' GC content differs by more
  than GC_DIFF_PCT points (both contigs >= GC_MIN_CONTIG_BP).
* hybrid_or_fusion: a flank gene present at both calls whose two protein
  copies are >= NEAR_IDENTICAL_PCT identical (allelic copies of one species).
  At >= IDENTICAL_PCT, assembly_artefact is also supported (a duplicated
  haplotype).
* homothallism: same_locus, or a shared flank gene whose two copies are
  intact and divergent (DIVERGENT_MIN_PCT to NEAR_IDENTICAL_PCT), as in
  Z. heterogamus. A degraded flank copy is NOT detected here.
"""
from __future__ import annotations

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

CAUSES = ("homothallism", "hybrid_or_fusion", "duplication",
          "mixed_culture_or_heterokaryon", "assembly_artefact")

DISTANT_BP = 50_000
GC_DIFF_PCT = 5.0
GC_MIN_CONTIG_BP = 5_000
NEAR_IDENTICAL_PCT = 95.0
IDENTICAL_PCT = 99.5
DIVERGENT_MIN_PCT = 55.0
MIN_PROTEIN_AA = 30

_UNDETERMINED = {None, "", "undetermined"}
_CONF_RANK = {"high": 3, "medium": 2, "low": 1}


def translate_exons(contig_seq: str, exons, strand: str, genetic_code: int) -> str:
    """Protein of a model from 1-based inclusive exon coordinates."""
    parts = [contig_seq[s - 1:e] for s, e in sorted(exons)]
    dna = "".join(parts)
    if strand == "-":
        dna = str(Seq(dna).reverse_complement())
    dna = dna[: len(dna) - len(dna) % 3]
    return str(Seq(dna).translate(table=genetic_code)).rstrip("*")


def _aligner() -> PairwiseAligner:
    al = PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10
    al.extend_gap_score = -0.5
    al.end_gap_score = 0
    return al


def protein_identity(a: str, b: str) -> float | None:
    """Identical residues over aligned (non-gap) columns, in percent."""
    if len(a) < MIN_PROTEIN_AA or len(b) < MIN_PROTEIN_AA:
        return None
    aln = _aligner().align(a, b)[0]
    same = cols = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for i, j in zip(range(s1, e1), range(s2, e2)):
            cols += 1
            same += a[i] == b[j]
    return round(100.0 * same / cols, 1) if cols else None


def _gc(seq: str) -> float | None:
    acgt = sum(seq.upper().count(c) for c in "ACGT")
    if acgt < GC_MIN_CONTIG_BP:
        return None
    return 100.0 * sum(seq.upper().count(c) for c in "GC") / acgt


def _best(calls: list[dict]) -> dict:
    return max(calls, key=lambda c: (_CONF_RANK.get(c["confidence"], 0), c.get("margin") or 0))


def _flanks(call: dict) -> list[dict]:
    return [g for g in call.get("gene_evidence") or [] if str(g.get("role", "")).startswith("flanking")]


def _call_summary(call: dict) -> dict:
    return {
        "idiomorph": call["idiomorph"], "contig": call["contig"],
        "start": call["start"], "end": call["end"],
        "confidence": call["confidence"], "locus_class": call.get("locus_class"),
        "classifier_input": call.get("classifier_input"), "margin": call.get("margin"),
        "flanks": [{"gene": g["gene"], "identity": g.get("identity"), "status": g.get("status")}
                   for g in _flanks(call)],
    }


def _arrangement(plus: dict, minus: dict) -> tuple[str, int | None]:
    if "homothallic_candidate" in (plus.get("locus_class"), minus.get("locus_class")):
        return "same_locus", 0
    if plus["contig"] != minus["contig"]:
        return "unlinked", None
    gap = max(0, max(plus["start"], minus["start"]) - min(plus["end"], minus["end"]))
    return ("same_contig_distant" if gap > DISTANT_BP else "same_locus"), gap


def _shared_flanks(plus: dict, minus: dict, seqs: dict[str, str], code: int) -> list[dict]:
    out = []
    by_name = {g["gene"]: g for g in _flanks(minus)}
    for g in _flanks(plus):
        h = by_name.get(g["gene"])
        if h is None or (g["contig"], g["start"]) == (h["contig"], h["start"]):
            continue
        ident = None
        if g["contig"] in seqs and h["contig"] in seqs and g.get("exons") and h.get("exons"):
            pa = translate_exons(seqs[g["contig"]], g["exons"], g["strand"], code)
            pb = translate_exons(seqs[h["contig"]], h["exons"], h["strand"], code)
            ident = protein_identity(pa, pb)
        out.append({"gene": g["gene"], "plus_contig": g["contig"], "minus_contig": h["contig"],
                    "protein_identity": ident})
    return out


def two_idiomorph_statements(
    calls: list[dict], enabled_families: set[str], contig_seqs: dict[str, str],
    *, genetic_code: int = 1,
) -> list[dict]:
    """One statement per enabled family that has a determined Plus and Minus call.

    `calls` use the report's shape (family, contig, start, end, idiomorph,
    confidence, locus_class, classifier_input, margin, gene_evidence).
    """
    statements = []
    families = sorted({c["family"] for c in calls if c["family"] in enabled_families})
    for fam in families:
        fam_calls = [c for c in calls if c["family"] == fam and c["idiomorph"] not in _UNDETERMINED]
        by_id: dict[str, list[dict]] = {}
        for c in fam_calls:
            by_id.setdefault(c["idiomorph"], []).append(c)
        if len(by_id) < 2:
            continue
        ids = sorted(by_id)
        a, b = ids[0], ids[1]
        plus, minus = _best(by_id.get("Plus", by_id[a])), _best(by_id.get("Minus", by_id[b]))
        arrangement, gap = _arrangement(plus, minus)
        support: dict[str, list[str]] = {c: [] for c in CAUSES}
        evidence: dict = {"distance_bp": gap}

        if arrangement == "same_locus":
            support["homothallism"].append(
                "both idiomorph genes at one locus (the Z. heterogamus arrangement)")
        dups = {k: len(v) for k, v in by_id.items() if len(v) > 1}
        if dups:
            support["duplication"].append(f"more than one call for {sorted(dups)}")
        evidence["calls_per_idiomorph"] = {k: len(v) for k, v in sorted(by_id.items())}

        gc_p = _gc(contig_seqs.get(plus["contig"], "")) if plus["contig"] in contig_seqs else None
        gc_m = _gc(contig_seqs.get(minus["contig"], "")) if minus["contig"] in contig_seqs else None
        if arrangement == "unlinked" and gc_p is not None and gc_m is not None:
            diff = round(abs(gc_p - gc_m), 2)
            evidence.update(gc_plus_contig=round(gc_p, 2), gc_minus_contig=round(gc_m, 2),
                            gc_difference_pct=diff)
            if diff > GC_DIFF_PCT:
                support["mixed_culture_or_heterokaryon"].append(
                    f"contig GC differs by {diff} points (> {GC_DIFF_PCT})")
        else:
            evidence["gc_difference_pct"] = None

        shared = _shared_flanks(plus, minus, contig_seqs, genetic_code)
        evidence["shared_flanks"] = shared
        for s in shared:
            ident = s["protein_identity"]
            if ident is None:
                continue
            if ident >= NEAR_IDENTICAL_PCT:
                support["hybrid_or_fusion"].append(
                    f"{s['gene']} copies {ident}% identical (allelic-like)")
                if ident >= IDENTICAL_PCT:
                    support["assembly_artefact"].append(
                        f"{s['gene']} copies {ident}% identical (possible duplicated haplotype)")
            elif ident >= DIVERGENT_MIN_PCT:
                support["homothallism"].append(
                    f"{s['gene']} copies intact and divergent ({ident}%), as in Z. heterogamus")

        statements.append({
            "family": fam,
            "arrangement": arrangement,
            "plus_call": _call_summary(plus),
            "minus_call": _call_summary(minus),
            "evidence": evidence,
            "possible_causes": [{"cause": c, "support": support[c]} for c in CAUSES],
            "supported_causes": [c for c in CAUSES if support[c]],
            "not_assessed": ["degraded (pseudogenised) flank copies",
                             "genome-wide duplicated single-copy genes", "read depth"],
        })
    return statements


def calls_from_results(results) -> list[dict]:
    """Adapt pipeline DetectionResults to the report shape used above."""
    from MATPredict.detect.report import _family_label

    out = []
    for r in results:
        clf = r.idiomorph_classifier or {}
        out.append({
            "family": _family_label(r.family_key), "contig": r.contig,
            "start": r.start, "end": r.end, "idiomorph": r.idiomorph,
            "confidence": r.confidence, "locus_class": r.locus_class,
            "classifier_input": clf.get("classifier_input"), "margin": clf.get("margin"),
            "gene_evidence": [
                {"gene": g.gene_name, "role": g.role, "contig": g.contig, "start": g.start,
                 "end": g.end, "strand": g.strand, "identity": g.identity, "status": g.status,
                 "exons": list(g.exons) if g.exons else [(g.start, g.end)]}
                for g in r.gene_evidence
            ],
        })
    return out


def calls_from_report(doc: dict) -> list[dict]:
    """Adapt a written detection_report.yaml's `detected` list."""
    out = []
    for x in doc.get("detected") or []:
        clf = x.get("idiomorph_classifier") or {}
        out.append({
            "family": x["family"], "contig": x["contig"], "start": x["start"], "end": x["end"],
            "idiomorph": x["idiomorph"], "confidence": x["confidence"],
            "locus_class": x.get("locus_class"),
            "classifier_input": clf.get("classifier_input"), "margin": clf.get("margin"),
            "gene_evidence": [
                {**g, "exons": [(e["start"], e["end"]) for e in (g.get("exons") or [])]
                 or [(g["start"], g["end"])]}
                for g in x.get("gene_evidence") or []
            ],
        })
    return out
