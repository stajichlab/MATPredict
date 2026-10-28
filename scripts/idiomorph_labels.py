"""Display labels for idiomorphs and MAT genes, so that a name says which reference system it belongs to.

A bare `M` or `P` is ambiguous (Mucoromycota sexM / sexP are Minus / Plus, fission yeast P and M are cell types). The label names the REFERENCE SYSTEM the idiomorph is defined in, not the species the
gene was found in: `SchPomM` is the M idiomorph of the Schizosaccharomyces pombe P/M system, wherever its homologue sits (Pneumocystis, Taphrina, Saitoella, Protomyces, Neolecta ...), in the same way
that the family name MATsc is the S. cerevisiae system. Labels are applied on output (tables, FASTA headers, tree labels, the reads pages); the database and the report YAML keep the raw names until a curator
ruling changes them (see docs/notes/2026-10-28_idiomorph-naming-proposal.md).
"""

SCHPOM = "SchPom"
# the P/M system: families Ascomycota:PM, Ascomycota:mat1, Ascomycota:mat2 (idiomorph names P and M)
IDIOMORPH = {("P", "M"): {"P": SCHPOM + "P", "M": SCHPOM + "M"}}
GENE = {"matPi": "SchPomPi", "mat1-Pi": "SchPomPi", "matPc": "SchPomPc", "mat1-Pc": "SchPomPc",
        "matMc": "SchPomMc", "mat1-Mc": "SchPomMc", "matMi": "SchPomMi", "mat1-Mi": "SchPomMi"}


def idiomorph_label(name, names=None):
    """Label of idiomorph `name` of a panel whose two idiomorph names are `names` (a pair, any order); the name itself when no rule applies."""
    if names is not None and tuple(sorted(names)) == ("M", "P"): return IDIOMORPH[("P", "M")].get(name, name)
    return name


def gene_label(name):
    """Label of a MAT gene name (matMc, mat1-Mc ... -> SchPomMc); the name itself when no rule applies."""
    return GENE.get(name, name)
