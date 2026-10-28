"""Display labels for idiomorphs and MAT genes, computed from scripts/idiomorph_label_registry.json (the single source of truth).

Grammar (curator, 2026-10-28): `<Lineage>__<System>_<Cell>`, for example `Taphr__SchPom_M`, `Sacch__CanAlb_alpha`, `Peziz__AspFum_MAT1-1`, `Mucor__PhyBla_plus`.
- Lineage: a fixed five-letter code (Peziz, Sacch, Taphr, Mucor, Basid).
- System: the frozen code of the REFERENCE ORGANISM of the system (three letters of the genus, three of the species), never a clade or a genus; for the Pezizomycotina MAT family it is chosen by class
  (Sordariomycetes NeuCra, Eurotiomycetes AspFum, ...; `NoRef` when the class has no curated reference). A label says which system's panel matched, not a homology or species claim.
- Cell: `a` / `alpha`, `A` / `B`, `P` / `M`, `plus` / `minus`, `MAT1-1` / `MAT1-2`, and `undet` for an undetermined idiomorph.
Gene labels are `<System>_<gene>` (`SchPom_Mc`, `CanAlb_alpha2`, `AspFum_MAT1-1-1`); flank genes keep their names; a gene only similar to a reference gene, not in a locus the pipeline reported, gets the
suffix `_like`. A family with no approved system code (Basidiomycota) returns raw names. Labels are applied on output only; the database and the report YAML keep the raw names.
Basidiomycota will use `<Lineage>__<System>_<Locus>[-<Allele>]` once their system codes are chosen.
"""
import json
import os
import re
from functools import lru_cache

_REGISTRY_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "idiomorph_label_registry.json")


@lru_cache(maxsize=1)
def registry():
    with open(_REGISTRY_PATH) as fh:
        return json.load(fh)


LINEAGES = registry()["lineages"]


def system_of(family, class_=None):
    """System code of a family (for the Pezizomycotina MAT family, of the class); None when the family has no approved code."""
    f = registry()["families"].get(family)
    if not f: return None
    if f.get("system_by_class"):
        return registry()[f["system_by_class"]].get(class_ or "", f["default_system"])
    return f.get("system")


def label_locus(family, idiomorph, class_=None):
    """{label, lineage, system, cell, status} of a locus called by `family` with raw idiomorph name `idiomorph`. status: labelled, pending (family has no system code yet), unmapped (family or name not in the registry)."""
    f = registry()["families"].get(family)
    if not f: return dict(label=idiomorph, lineage="", system="", cell="", status="unmapped")
    sys_ = system_of(family, class_)
    if sys_ is None: return dict(label=idiomorph, lineage=f["lineage"], system="", cell="", status="pending")
    cell = f["cells"].get(idiomorph)
    if cell is None: return dict(label=idiomorph, lineage=f["lineage"], system=sys_, cell="", status="unmapped")
    return dict(label=f"{f['lineage']}__{sys_}_{cell}", lineage=f["lineage"], system=sys_, cell=cell, status="labelled")


def gene_label(gene, basis="called", family=None, class_=None):
    """Label of MAT gene `gene` of `family` (`<System>_<gene>`; `_like` when `basis` is "similarity_only"); the name itself for a flank gene or an unknown gene.
    Without `family` a name is mapped only if it is unambiguous across the registry (matPi -> SchPom_Pi; MATa1 is not: it is KomPha and OgaPol)."""
    fams = registry()["families"]
    cands = [family] if family else [k for k, f in fams.items() if gene in f.get("genes", {})]
    out = set()
    for k in cands:
        f = fams.get(k) or {}
        sys_ = system_of(k, class_)
        if sys_ is None: continue
        if gene in f.get("genes", {}): out.add(f"{sys_}_{f['genes'][gene]}")
        elif f.get("gene_pattern") and re.match(f["gene_pattern"], gene): out.add(f"{sys_}_{gene}")
    if len(out) != 1: return gene
    return out.pop() + ("_like" if basis == "similarity_only" else "")


def idiomorph_label(name, names=None):
    """Panel-level label for a panel whose two idiomorph names are `names`: only the P/M pair maps (to Taphr__SchPom_P / _M); the name itself otherwise."""
    if names is not None and tuple(sorted(names)) == ("M", "P") and name in ("P", "M"): return label_locus("Ascomycota:PM", name)["label"]
    return name


def lineage_of(names):
    """Lineage token of a panel with idiomorph names `names` (Taphr for the P/M pair); '' when unknown."""
    return "Taphr" if names is not None and tuple(sorted(names)) == ("M", "P") else ""


def idiomorph_of_gene(name):
    """P or M for a P/M-system gene name (matMc, mat1-Pi ...), '' otherwise."""
    return {"Pi": "P", "Pc": "P", "Mc": "M", "Mi": "M"}.get(registry()["families"]["Ascomycota:PM"]["genes"].get(name, ""), "")
