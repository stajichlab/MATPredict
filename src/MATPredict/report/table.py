"""`matpredict report table`: genomes.tsv (one row per genome) and loci.tsv (one row per called locus) from a directory of detect runs."""
from __future__ import annotations

import csv
from pathlib import Path

import yaml

_Loader = getattr(yaml, "CSafeLoader", yaml.SafeLoader)  # the C loader is ~10x faster on thousands of reports

REPORT_YAML = "detection_report.yaml"


def _genome_row(sample: str, doc: dict) -> dict:
    detected = doc.get("detected") or []
    idiomorphs = sorted({d["idiomorph"] for d in detected if d.get("idiomorph")})
    run = doc.get("run") or {}
    genome = run.get("genome") or {}
    return {
        "sample": sample,
        "status": "report",
        "reason": "",
        "organism": run.get("organism") or "",
        "taxid": run.get("taxid") or "",
        "phylum": run.get("phylum") or "",
        "contigs": genome.get("contigs"),
        "length_bp": genome.get("length_bp"),
        "n50": genome.get("n50"),
        "wall_seconds": run.get("wall_seconds"),
        "genetic_code": doc.get("genetic_code"),
        "families_attempted": ";".join(doc.get("families_attempted") or []),
        "routing_mode": doc.get("routing_mode"),
        "not_searched_reason": doc.get("not_searched_reason") or "",
        "n_loci": len(detected),
        "idiomorphs": ";".join(idiomorphs),
        # signals, not verdicts: a genome with one idiomorph is not shown to be heterothallic
        "two_idiomorphs": ";".join(str(t.get("arrangement")) for t in doc.get("two_idiomorphs") or []),
        "homothallic_candidate": int(any(d.get("locus_class") == "homothallic_candidate" for d in detected)),
        "zygosity": doc.get("zygosity") or "",
        "n_withheld": len(doc.get("suppressed_loci") or []),
    }


def _missing_row(run_dir: Path) -> dict:
    """A run folder with no report: failed, skipped, or nothing written."""
    for status, name in (("failed", "failed_reason.txt"), ("skipped", "skipped_reason.txt")):
        if (run_dir / name).is_file():
            return {"sample": run_dir.name, "status": status, "reason": (run_dir / name).read_text().strip()}
    return {"sample": run_dir.name, "status": "no_report", "reason": ""}


def _span(span: dict | None) -> str:
    return f"{span['start']}-{span['end']}" if isinstance(span, dict) and span.get("start") is not None else ""


def _locus_row(sample: str, d: dict) -> dict:
    supported = d.get("supported_span") if isinstance(d.get("supported_span"), dict) else {}
    return {
        "sample": sample,
        "family": d.get("family"),
        "contig": d.get("contig"),
        "start": d.get("start"),
        "end": d.get("end"),
        "idiomorph": d.get("idiomorph") or "",
        "idiomorph_class": d.get("idiomorph_class") or "",
        "confidence": d.get("confidence"),
        "locus_class": d.get("locus_class"),
        "detection_pass": d.get("detection_pass"),
        "core_span": _span(d.get("core_span")),
        "supported_span": _span(supported),
        "beyond_supported_bp": supported.get("beyond_supported_bp", ""),
        "genes_found": ";".join(map(str, d.get("genes_found") or [])),
        "genes_missing": ";".join(map(str, d.get("genes_missing") or [])),
        "polished_genes": d.get("polished_genes"),
        "merged_from": ";".join(map(str, d.get("merged_from") or [])),
        "fragmented": d.get("fragmented"),
    }


def read_samples(path: Path, key: str) -> list[dict]:
    """Rows of a CSV or TSV sample sheet; `key` names the column that matches the run folder names."""
    with open(path, newline="") as f:
        delim = "\t" if str(path).endswith((".tsv", ".txt")) else ","
        rows = list(csv.DictReader(f, delimiter=delim))
    if rows and key not in rows[0]:
        raise ValueError(f"{path} has no column {key!r} (columns: {', '.join(rows[0])})")
    return rows


def build_tables(runs: Path, samples: list[dict] | None = None, key: str = "ASMID") -> tuple[list[dict], list[dict]]:
    """(genome rows, locus rows) for every run folder under `runs`, in folder-name order.

    With a sample sheet, its columns join onto each genome row and a sheet row with no run folder gets status `not_run`.
    """
    genomes: list[dict] = []
    loci: list[dict] = []
    for run_dir in sorted(p for p in Path(runs).iterdir() if p.is_dir()):
        if not (run_dir / REPORT_YAML).is_file():
            genomes.append(_missing_row(run_dir))
            continue
        doc = yaml.load((run_dir / REPORT_YAML).read_text(), Loader=_Loader) or {}
        genomes.append(_genome_row(run_dir.name, doc))
        loci += [_locus_row(run_dir.name, d) for d in doc.get("detected") or []]
    if samples is not None:
        meta = {r[key]: {k: v for k, v in r.items() if k != key} for r in samples}
        have = {g["sample"] for g in genomes}
        genomes += [{"sample": name, "status": "not_run", "reason": "no run folder"} for name in meta if name not in have]
        genomes = [{**g, **meta.get(g["sample"], {})} for g in genomes]
    return genomes, loci


def write_tsv(path: Path, rows: list[dict]) -> None:
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
