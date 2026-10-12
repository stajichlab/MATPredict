"""`matpredict report summary`: recovery and mating-type counts by group, from genomes.tsv and loci.tsv."""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

NOT_SEARCHED = "not_searched"  # routing_mode of a genome whose phylum has no curated MAT family


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def summarize(genomes: list[dict], loci: list[dict], by: str) -> list[dict]:
    """One row per value of column `by` (plus an ALL row): how many genomes were searched, called, and what was called."""
    if genomes and by not in genomes[0]:
        raise ValueError(f"genomes.tsv has no column {by!r}; join a sample sheet that has it (columns: {', '.join(genomes[0])})")
    loci_by_sample = defaultdict(list)
    for l in loci:
        loci_by_sample[l["sample"]].append(l)
    groups: dict[str, list[dict]] = defaultdict(list)
    for g in genomes:
        groups[g.get(by) or "(none)"].append(g)
    rows = [_group_row(by, name, members, loci_by_sample) for name, members in sorted(groups.items())]
    rows.append(_group_row(by, "ALL", genomes, loci_by_sample))
    return rows


def _group_row(by: str, name: str, members: list[dict], loci_by_sample: dict) -> dict:
    status = Counter(g["status"] for g in members)
    not_searched = sum(1 for g in members if g["status"] == "report" and g.get("routing_mode") == NOT_SEARCHED)
    searched = [g for g in members if g["status"] == "report" and g.get("routing_mode") != NOT_SEARCHED]
    called = [g for g in searched if int(g.get("n_loci") or 0) > 0]
    idiomorphs = Counter(l["idiomorph"] for g in called for l in loci_by_sample[g["sample"]] if l.get("idiomorph"))
    both = [g for g in called if g.get("two_idiomorphs") or len({l["idiomorph"] for l in loci_by_sample[g["sample"]] if l.get("idiomorph")}) > 1]
    high = [g for g in called if any(l.get("confidence") == "high" for l in loci_by_sample[g["sample"]])]
    return {
        by: name,
        "genomes": len(members),
        "searched": len(searched),
        "called": len(called),
        "recovery": f"{len(called) / len(searched):.2f}" if searched else "",
        "called_high_confidence": len(high),
        "not_searched": not_searched,
        "failed": status["failed"],
        "skipped": status["skipped"],
        "no_report": status["no_report"],
        "not_run": status["not_run"],
        "both_idiomorphs": len(both),
        "homothallic_candidate": sum(1 for g in called if str(g.get("homothallic_candidate")) == "1"),
        "idiomorph_loci": ";".join(f"{k}:{v}" for k, v in idiomorphs.most_common()),
    }


def markdown(rows: list[dict], by: str) -> str:
    cols = list(rows[0]) if rows else [by]
    lines = [f"# MAT detection summary by {by}", "", "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(str(r[c]) for c in cols) + " |" for r in rows]
    lines += ["", "recovery = called / searched. A genome with one idiomorph is not shown to be heterothallic; `both_idiomorphs` and "
              "`homothallic_candidate` are signals from the detect report, not verdicts.", ""]
    return "\n".join(lines)
