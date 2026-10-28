"""Display labels for idiomorphs and MAT genes, so that a name says which lineage and reference system it belongs to.

Grammar (curator, 2026-10-28): `<Lineage>__<System>_<Cell>`, for example `Taphr__SchPom_M`.
- Lineage: the subphylum (`Taphr` Taphrinomycotina, `Sacch` Saccharomycotina, the full word `Mucoromycotina` for that group; Pezizomycotina MAT1-1 / MAT1-2 stay as they are).
- System: a code for the REFERENCE ORGANISM of the system (three letters of the genus and three of the species), not for the species the gene was found in and not a clade: `SchPom` is the
  Schizosaccharomyces pombe P/M system, wherever its homologue sits (Pneumocystis, Taphrina, Saitoella, Protomyces, Neolecta ...), like the family name MATsc is the S. cerevisiae system.
  A label says which panel or reference system matched, not a homology or species claim.
- Cell: P / M for the P/M system (upper case, as in S. pombe), `a` / `alpha` (lower case) for a/alpha systems.
Gene labels carry the system code and the gene only (`SchPom_Mc`), the lineage being in its own column. A gene that is only similar to a reference gene, not called by the pipeline in a locus, gets the suffix `_like`
(`SchPom_Mc_like`), so a candidate is not read as an established member of the system.

Labels are applied on output (tables, FASTA headers, tree labels, the reads pages); the database and the report YAML keep the raw names until a curator ruling changes them
(docs/notes/2026-10-28_idiomorph-naming-proposal.md). Only the P/M system is mapped so far.
"""

LINEAGE_PM, SYSTEM_PM = "Taphr", "SchPom"
GENE_PM = {"matPi": "Pi", "mat1-Pi": "Pi", "matPc": "Pc", "mat1-Pc": "Pc", "matMc": "Mc", "mat1-Mc": "Mc", "matMi": "Mi", "mat1-Mi": "Mi"}


def idiomorph_label(name, names=None):
    """Label of idiomorph `name` of a panel whose two idiomorph names are `names` (a pair, any order); the name itself when no rule applies."""
    if names is not None and tuple(sorted(names)) == ("M", "P") and name in ("P", "M"): return f"{LINEAGE_PM}__{SYSTEM_PM}_{name}"
    return name


def lineage_of(names):
    """Lineage token of a panel with the idiomorph names `names`; '' when unknown."""
    return LINEAGE_PM if names is not None and tuple(sorted(names)) == ("M", "P") else ""


def gene_label(name, basis="called"):
    """Label of a MAT gene name (matMc, mat1-Mc ... -> SchPom_Mc); `_like` is added when `basis` is "similarity_only". The name itself when no rule applies."""
    g = GENE_PM.get(name)
    if g is None: return name
    return f"{SYSTEM_PM}_{g}" + ("_like" if basis == "similarity_only" else "")


def idiomorph_of_gene(name):
    """P or M for a P/M-system gene name, '' otherwise."""
    g = GENE_PM.get(name, "")
    return g[0] if g else ""
