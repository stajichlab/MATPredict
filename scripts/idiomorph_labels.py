"""Display labels for idiomorphs and MAT genes, computed from scripts/idiomorph_label_registry.json (the single source of truth).

Grammar (curator, 2026-10-28): `<Lineage>__<System>_<Cell>`, for example `Taphr__SchPom_M`, `Sacch__CanAlb_alpha`, `Peziz__AspFum_MAT1-1`, `Mucor__PhyBla_plus`.
- Lineage: a fixed five-letter code (Peziz, Sacch, Taphr, Mucor, Basid).
- System: the frozen code of the REFERENCE ORGANISM of the system (three letters of the genus, three of the species), never a clade or a genus; for the Pezizomycotina MAT family it is chosen by class
  (Sordariomycetes NeuCra, Eurotiomycetes AspFum, ...; `NoRef` when the class has no curated reference). A label says which system's panel matched, not a homology or species claim.
- Cell (or Basidiomycota locus and allele): `a` / `alpha`, `A` / `B`, `P` / `M`, `plus` / `minus`, `MAT1-1` / `MAT1-2`, and `undet` for an undetermined idiomorph.
Gene labels are `<System>_<gene>` (`SchPom_Mc`, `CanAlb_alpha2`, `AspFum_MAT1-1-1`); flank genes keep their names; a gene only similar to a reference gene, not in a locus the pipeline reported, gets the
suffix `_like`. A family with no approved system code (Basidiomycota) returns raw names. Labels are applied on output only; the database and the report YAML keep the raw names.
Basidiomycota will use `<Lineage>__<System>_<Locus>[-<Allele>]` once their system codes are chosen.
"""
import hashlib
import json
import os
import re
import sys
from functools import lru_cache

_REGISTRY_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "idiomorph_label_registry.json")


@lru_cache(maxsize=1)
def registry():
    with open(_REGISTRY_PATH) as fh:
        return json.load(fh)


LINEAGES = registry()["lineages"]


def _canonical(reg):
    return json.dumps({k: v for k, v in reg.items() if k not in ("_comment", "registry_hash")}, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def registry_hash(reg=None):
    """First 8 hex characters of the SHA-256 of the canonical registry (without `_comment` and `registry_hash`)."""
    return hashlib.sha256(_canonical(reg or registry()).encode()).hexdigest()[:8]


def registry_version():
    """`<registry_date>+<registry_hash>` as stamped in the registry (curator, 2026-10-29, #210). The date is bumped by hand when a rule change is approved;
    `registry_problems()` (and the tests) fail when the stamped hash no longer matches the content."""
    reg = registry()
    return f"{reg['registry_date']}+{reg['registry_hash']}"


def registry_problems(reg=None):
    """Problems that stop the registry being released: a stale hash, a frozen code removed or renamed without an alias, an alias to a code that does not exist."""
    reg = reg or registry()
    out = []
    if reg.get("registry_hash") != registry_hash(reg):
        out.append(f"registry_hash is {reg.get('registry_hash')!r} but the content hashes to {registry_hash(reg)!r}: bump registry_date and run `python3 scripts/idiomorph_labels.py --stamp`")
    for kind in ("lineages", "systems"):
        aliases = reg.get("aliases", {}).get(kind, {})
        for code in reg.get("frozen_codes", {}).get(kind, []):
            if code in reg[kind]:
                continue
            if code not in aliases:
                out.append(f"frozen {kind[:-1]} code {code!r} was removed or renamed without an alias")
        for old, new in aliases.items():
            if new not in reg[kind]:
                out.append(f"alias {old!r} -> {new!r}: {new!r} is not a {kind[:-1]} code in the registry")
    return out


def canonical_code(kind, code):
    """The current code for `code` (`kind` is "lineages" or "systems"): an alias resolves to the code that replaced it."""
    return registry().get("aliases", {}).get(kind, {}).get(code, code)


# label_status (curator, 2026-10-29, #211): ONE value per locus, the first that applies in this order; `label_reasons` lists every one that applies.
STATUS_ORDER = ("family_unconfirmed", "cross_call", "no_reference", "similarity_only", "undetermined", "ok")


def locus_status(family, idiomorph, class_=None, *, family_unconfirmed=False, cross_call=False, similarity_only=False):
    """{label_status, label_reasons} of a locus. family_unconfirmed and cross_call come from the report or the table (the registry cannot know them);
    no_reference: the system is the NoRef placeholder or the family has no label mapping; undetermined: the label ends in `_undet` (a Basidiomycota locus label without an allele is complete, not undetermined)."""
    lab = label_locus(family, idiomorph, class_)
    reasons = []
    if family_unconfirmed: reasons.append("family_unconfirmed")
    if cross_call: reasons.append("cross_call")
    if lab["status"] != "labelled" or lab["system"] == "NoRef": reasons.append("no_reference")
    if similarity_only: reasons.append("similarity_only")
    if lab["cell"] == "undet": reasons.append("undetermined")      # the label carries `_undet`
    return dict(label_status=reasons[0] if reasons else "ok", label_reasons=";".join(reasons))


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
    if f.get("locus"):      # multiallelic Basidiomycota system: <Lineage>__<System>_<Locus>[-<Allele>], the allele being the raw idiomorph unless undetermined
        allele = "" if idiomorph in ("", "undetermined") else "-" + idiomorph
        return dict(label=f"{f['lineage']}__{sys_}_{f['locus']}{allele}", lineage=f["lineage"], system=sys_, cell=f["locus"] + allele, status="labelled")
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


if __name__ == "__main__":        # python3 scripts/idiomorph_labels.py [--stamp]
    if "--stamp" in sys.argv:
        reg = registry()
        reg["registry_hash"] = registry_hash(reg)
        with open(_REGISTRY_PATH, "w") as fh:
            fh.write(json.dumps(reg, indent=2, ensure_ascii=False) + "\n")
    registry.cache_clear()
    print(registry_version(), *registry_problems(), sep="\n")
