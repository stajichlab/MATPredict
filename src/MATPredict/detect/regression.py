"""Pre-sign-off regression check: diff two sets of detect reports.

Curator ruling, J. Stajich, 2026-09-30: before a new record or a classifier
rebuild is signed off, a replay of existing calls must list every call whose
gene model, label, confidence or presence changed, so the curator can judge
each change. The motivating cases (results/2026-09-29_umbelopsis_merge/): a new
S. racemosum reference changed the gene model detect built for Circinella
minor (classifier 86.4 vs 111.8 bits) with no change to the classifier itself,
and new Umbelopsis flank references pushed M. pusillus out of the polish cap.

A code rule cannot tell a better model from a worse one, so this module does
not judge; it makes every change visible. Loci are matched per genome by
family + contig + span overlap. Every locus (called or withheld) on either side
gets one row; `change_types` is empty only when nothing listed below changed.
Classifier margin shifts below `score_delta` bits are recorded in the row
(`margin_delta`) but are not a change type -- rebuild-level shifts of a few
bits are expected (results/2026-09-29_aligner_default/: -8.1..+11.3 bits).
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import yaml

#: Default classifier margin shift, in bits, that counts as a change.
DEFAULT_SCORE_DELTA = 5.0

COLUMNS = [
    "genome", "family", "contig", "base_span", "cand_span", "base_state", "cand_state",
    "base_idiomorph", "cand_idiomorph", "base_confidence", "cand_confidence",
    "base_locus_class", "cand_locus_class", "base_verification", "cand_verification",
    "base_clf_input", "cand_clf_input", "base_margin", "cand_margin", "margin_delta",
    "base_best_score", "cand_best_score", "base_genes", "cand_genes", "model_changes",
    "change_types",
]


def load_reports(root: Path) -> dict[str, dict]:
    """`{genome_key: report}` for every detection_report.yaml under `root`.

    The key is the report's directory relative to `root` with `runs` parts
    dropped, so Zygo's `scaffold/runs/<org>` and `contig/runs/<org>` stay
    distinct while a plain `runs/<asmid>` layout keys by the assembly id.
    """
    root = Path(root)
    out: dict[str, dict] = {}
    for f in sorted(root.rglob("detection_report.yaml")):
        parts = [p for p in f.parent.relative_to(root).parts if p != "runs"]
        key = "/".join(parts) or f.parent.name
        rep = yaml.safe_load(f.read_text()) or {}
        rep["_capped"] = _capped_clusters(f.parent / "evidence_diagnostics.jsonl")
        out[key] = rep
    return out


def _capped_clusters(path: Path) -> list[tuple]:
    if not path.exists():
        return []
    capped = []
    for line in path.read_text().splitlines():
        try:
            row = json.loads(line)
        except ValueError:
            continue
        if row.get("kind") == "evidence" and row.get("polish_capped"):
            capped.append((row.get("family"), row.get("contig"),
                           row.get("cluster_start"), row.get("cluster_end")))
    return capped


def _loci(rep: dict) -> list[dict]:
    loci = []
    for c in rep.get("detected") or []:
        loci.append({**c, "_state": "called"})
    for s in rep.get("suppressed_loci") or []:
        reason = s.get("withheld_reason") or "unspecified"
        if any(f == s.get("family") and ct == s.get("contig")
               and _overlap(s.get("start"), s.get("end"), a, b) > 0
               for f, ct, a, b in rep.get("_capped") or []):
            reason += "+polish_capped"
        loci.append({**s, "_state": f"withheld:{reason}"})
    return loci


def _overlap(a1, a2, b1, b2) -> int:
    try:
        return max(0, min(int(a2), int(b2)) - max(int(a1), int(b1)))
    except (TypeError, ValueError):
        return 0


def _pair(base: list[dict], cand: list[dict]) -> list[tuple]:
    """Greedy best-overlap matching within family + contig; called loci first."""
    pairs, used = [], set()
    order = sorted(range(len(base)), key=lambda i: base[i]["_state"] != "called")
    for i in order:
        b, best, best_j = base[i], 0, None
        for j, c in enumerate(cand):
            if j in used or c.get("family") != b.get("family") or c.get("contig") != b.get("contig"):
                continue
            ov = _overlap(b.get("start"), b.get("end"), c.get("start"), c.get("end"))
            if ov > best:
                best, best_j = ov, j
        if best_j is None:
            pairs.append((b, None))
        else:
            used.add(best_j)
            pairs.append((b, cand[best_j]))
    pairs.extend((None, c) for j, c in enumerate(cand) if j not in used)
    return pairs


def _clf(locus: dict | None) -> dict:
    return (locus or {}).get("idiomorph_classifier") or {}


def _best(locus):
    scores = _clf(locus).get("scores") or {}
    return max(scores.values()) if scores else None


def _core_models(locus: dict | None) -> dict[str, tuple]:
    out = {}
    for g in (locus or {}).get("gene_evidence") or []:
        if g.get("role") != "core_MAT":
            continue
        exons = tuple((e.get("start"), e.get("end")) for e in g.get("exons") or [])
        length = sum(abs(int(e) - int(s)) + 1 for s, e in exons if s is not None and e is not None)
        out[g.get("gene")] = (exons, length, g.get("identity"))
    return out


def _model_changes(b, c) -> str:
    # Withheld loci carry no gene_evidence in the report, so a model cannot be
    # compared: say so instead of reporting the genes as removed or added.
    if not (b or {}).get("gene_evidence") or not (c or {}).get("gene_evidence"):
        if (b or {}).get("gene_evidence") or (c or {}).get("gene_evidence"):
            side = "candidate" if (b or {}).get("gene_evidence") else "baseline"
            return f"no_gene_evidence_in_{side}"
        return ""
    bm, cm = _core_models(b), _core_models(c)
    notes = []
    for gene in sorted(set(bm) | set(cm)):
        if gene not in bm or gene not in cm:
            notes.append(f"{gene}:{'added' if gene in cm else 'removed'}")
        elif bm[gene][0] != cm[gene][0]:
            notes.append(f"{gene}:{bm[gene][1]}bp/{bm[gene][2]}%->{cm[gene][1]}bp/{cm[gene][2]}%")
    return ";".join(notes)


def _span(locus):
    return "" if locus is None else f"{locus.get('start')}-{locus.get('end')}"


def compare_runs(base: dict[str, dict], cand: dict[str, dict],
                 score_delta: float = DEFAULT_SCORE_DELTA) -> list[dict]:
    rows = []
    for genome in sorted(set(base) | set(cand)):
        if genome not in cand:
            rows.append(_row(genome, None, None, ["genome_missing_in_candidate"]))
            continue
        if genome not in base:
            rows.append(_row(genome, None, None, ["genome_missing_in_baseline"]))
            continue
        for b, c in _pair(_loci(base[genome]), _loci(cand[genome])):
            rows.append(_row(genome, b, c, _change_types(b, c, score_delta)))
    return rows


def _change_types(b, c, score_delta) -> list[str]:
    bs = b["_state"] if b else "absent"
    cs = c["_state"] if c else "absent"
    ch = []
    if bs == "called" and cs != "called":
        ch.append("call_lost")
    elif bs != "called" and cs == "called":
        ch.append("call_gained")
    elif bs != cs:
        ch.append("withheld_reason_changed")
    if b is None or c is None:
        return ch
    for key, name in (("idiomorph", "idiomorph_changed"), ("confidence", "confidence_changed"),
                      ("locus_class", "locus_class_changed"),
                      ("verification", "verification_changed")):
        if bs == cs == "called" and b.get(key) != c.get(key):
            ch.append(name)
    if (b.get("start"), b.get("end")) != (c.get("start"), c.get("end")):
        ch.append("span_changed")
    if sorted(b.get("genes_found") or []) != sorted(c.get("genes_found") or []):
        ch.append("gene_set_changed")
    mc = _model_changes(b, c)
    if mc and not mc.startswith("no_gene_evidence"):
        ch.append("core_model_changed")
    bi, ci = _clf(b).get("classifier_input"), _clf(c).get("classifier_input")
    if bi != ci:
        ch.append("classifier_input_changed")
    bmg, cmg = _clf(b).get("margin"), _clf(c).get("margin")
    if bmg is not None and cmg is not None and abs(cmg - bmg) >= score_delta:
        ch.append("classifier_shift")
    return ch


def _row(genome, b, c, changes) -> dict:
    locus = b or c or {}
    bmg, cmg = _clf(b).get("margin"), _clf(c).get("margin")
    return {
        "genome": genome, "family": locus.get("family", ""), "contig": locus.get("contig", ""),
        "base_span": _span(b), "cand_span": _span(c),
        "base_state": b["_state"] if b else "absent", "cand_state": c["_state"] if c else "absent",
        "base_idiomorph": (b or {}).get("idiomorph", ""), "cand_idiomorph": (c or {}).get("idiomorph", ""),
        "base_confidence": (b or {}).get("confidence", ""), "cand_confidence": (c or {}).get("confidence", ""),
        "base_locus_class": (b or {}).get("locus_class", ""), "cand_locus_class": (c or {}).get("locus_class", ""),
        "base_verification": (b or {}).get("verification") or "",
        "cand_verification": (c or {}).get("verification") or "",
        "base_clf_input": _clf(b).get("classifier_input", ""), "cand_clf_input": _clf(c).get("classifier_input", ""),
        "base_margin": bmg if bmg is not None else "", "cand_margin": cmg if cmg is not None else "",
        "margin_delta": round(cmg - bmg, 1) if bmg is not None and cmg is not None else "",
        "base_best_score": _best(b) if _best(b) is not None else "",
        "cand_best_score": _best(c) if _best(c) is not None else "",
        "base_genes": ",".join(sorted((b or {}).get("genes_found") or [])),
        "cand_genes": ",".join(sorted((c or {}).get("genes_found") or [])),
        "model_changes": _model_changes(b, c) if b and c else "",
        "change_types": ",".join(changes),
    }


def write_outputs(rows: list[dict], out_dir: Path, *, title: str,
                  score_delta: float = DEFAULT_SCORE_DELTA) -> tuple[Path, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv = out_dir / "regression_diff.tsv"
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    changed = [r for r in rows if r["change_types"]]
    on_calls = [r for r in changed if "called" in (r["base_state"], r["cand_state"])]
    withheld_only = [r for r in changed if r not in on_calls]
    small = sum(1 for r in rows if not r["change_types"] and r["margin_delta"] not in ("", 0, 0.0))
    genomes = len({r["genome"] for r in rows})

    def _table(sub):
        counts = Counter(c for r in sub for c in r["change_types"].split(","))
        return (["| change | loci |", "|---|---|"]
                + ([f"| {k} | {v} |" for k, v in sorted(counts.items())] or ["| none | 0 |"]))

    def _line(r):
        return (f"- **{r['genome']}** {r['family']} {r['contig']} "
                f"[{r['change_types']}]: {r['base_state']} {r['base_idiomorph']}/{r['base_confidence']} "
                f"-> {r['cand_state']} {r['cand_idiomorph']}/{r['cand_confidence']}; "
                f"margin {r['base_margin']} -> {r['cand_margin']}; "
                f"best score {r['base_best_score']} -> {r['cand_best_score']}"
                + (f"; models {r['model_changes']}" if r["model_changes"] else "")
                + (f"; genes {r['base_genes']} -> {r['cand_genes']}"
                   if "gene_set_changed" in r["change_types"] else ""))

    lines = [f"# Regression check: {title}", "",
             f"- Genomes compared: {genomes}; loci rows: {len(rows)}; loci with a change: "
             f"{len(changed)} ({len(on_calls)} touch a call, {len(withheld_only)} are withheld on "
             f"both sides).",
             f"- Classifier margin shifts below {score_delta} bits are not listed as changes; "
             f"{small} unchanged loci carry such a shift (see `margin_delta` in the TSV).",
             "- Every withheld-only change is listed in `regression_withheld_changes.md`.",
             "", "## Change types (loci that touch a call)", ""] + _table(on_calls)
    lines += ["", "## Changes to calls", ""] + ([_line(r) for r in on_calls] or ["None."])
    md = out_dir / "regression_summary.md"
    md.write_text("\n".join(lines) + "\n")
    app = [f"# Withheld-only changes: {title}", "",
           "Loci withheld (or absent) on both sides whose span, genes or withheld reason changed.",
           "", "## Change types", ""] + _table(withheld_only)
    app += ["", "## Loci", ""] + ([_line(r) for r in withheld_only] or ["None."])
    (out_dir / "regression_withheld_changes.md").write_text("\n".join(app) + "\n")
    return tsv, md
