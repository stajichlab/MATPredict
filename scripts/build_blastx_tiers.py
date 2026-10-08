#!/usr/bin/env python3
"""Build tiered protein references for read blastx typing from db/ records, for a target species.

usage: build_blastx_tiers.py DB_ASCOMYCOTA OUT_DIR TARGET_GENUS TARGET_ORDER TARGET_SPECIES
Tiers (each holds only that level, so accuracy can be read against taxonomic distance):
  T1 same species, T2 same genus other species, T3 same order other genus, T4 Pezizomycotina outside the order.
Header: >idiomorph|gene|species|role   (idiomorph MAT1-1 / MAT1-2 / control)
"""
import glob
import sys
from pathlib import Path

import yaml

CORE = {"MAT1-1-1", "MAT1-1-2", "MAT1-1-3", "MAT1-2-1"}
CONTROL = {"APN2", "SLA2"}


def read_fasta(path):
    out, name, seq = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(seq)))
            name, seq = dict(x.split("=", 1) for x in line[1:].split("|")[1:])["name"], []
        elif name is not None:
            seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def records(db):
    for meta in glob.glob(f"{db}/*/*/metadata.yaml"):
        d = Path(meta).parent
        y = yaml.safe_load(open(meta))
        lin = dict(x.split("__", 1) for x in y["taxonomy"]["lineage"].split(";"))
        yield y["organism"]["species"], lin, read_fasta(d / "proteins.faa")


def main(db, out, genus, order, species):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    tiers = {t: [] for t in ("T1", "T2", "T3", "T4")}
    controls = []   # APN2/SLA2 from every Pezizomycotina record: depth control, added to every tier
    for sp, lin, genes in records(db):
        if "Pezizomycotina" not in lin.get("subphylum", ""):
            continue
        if sp == species:
            t = "T1"
        elif lin.get("g") == genus:
            t = "T2"
        elif lin.get("o") == order:
            t = "T3"
        else:
            t = "T4"
        for name, seq in genes:
            if name in CORE:
                kind = "MAT1-1" if name.startswith("MAT1-1") else "MAT1-2"
            elif name in CONTROL:
                controls.append((f"control|{name}|{sp.replace(' ', '_')}", seq))
                continue
            else:
                continue
            tiers[t].append((f"{kind}|{name}|{sp.replace(' ', '_')}", seq))
    for t in tiers:
        tiers[t] += controls
    for t, items in tiers.items():
        with open(out / f"{t}.faa", "w") as fh:
            for h, s in items:
                fh.write(f">{h}\n{s}\n")
        print(t, len(items), "proteins;", sum(1 for h, _ in items if h.startswith("MAT1-1")), "MAT1-1,",
              sum(1 for h, _ in items if h.startswith("MAT1-2")), "MAT1-2,", sum(1 for h, _ in items if h.startswith("control")), "control")


if __name__ == "__main__":
    main(*sys.argv[1:6])
