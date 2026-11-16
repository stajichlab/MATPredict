"""Tests for scripts/flank_idiomorph_probe.py (classify_gene)."""
import importlib.util
import pathlib

SCRIPTS = pathlib.Path(__file__).resolve().parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("flank_idiomorph_probe", SCRIPTS / "flank_idiomorph_probe.py")
fip = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fip)

# genome B: contig c1 genes b0..b29 at index 0..29
IDX_B = {f"b{i}": ("c1", i) for i in range(30)}
FLANK = [f"a{i}" for i in range(1, 6)]
ORTH = {f"a{i}": f"b{10 + i}" for i in range(1, 6)}  # flank maps to b11..b15


def test_conserved_when_x_has_ortholog_in_region():
    orth = dict(ORTH, x="b12x")
    idx = dict(IDX_B, b12x=("c1", 12))
    assert fip.classify_gene("x", FLANK, orth, idx, set(), 3, 25)[0] == "conserved"


def test_absent_when_region_has_no_domain_gene():
    assert fip.classify_gene("x", FLANK, ORTH, IDX_B, set(), 3, 25)[:2] == ("absent", 5)


def test_other_domain_when_region_has_unrelated_domain_gene():
    assert fip.classify_gene("x", FLANK, ORTH, IDX_B, {"b13"}, 3, 25)[0] == "other_domain"


def test_domain_gene_far_from_region_does_not_count():
    assert fip.classify_gene("x", FLANK, ORTH, IDX_B, {"b29"}, 3, 25)[0] == "absent"


def test_no_synteny_with_too_few_flank_hits():
    orth = {"a1": "b11", "a2": "b12"}
    assert fip.classify_gene("x", FLANK, orth, IDX_B, set(), 3, 25)[0] == "no_synteny"


def test_flank_genes_split_over_contigs_are_not_one_block():
    idx = {f"b{i}": ("c1" if i < 12 else "c2", i) for i in range(30)}
    assert fip.classify_gene("x", FLANK, ORTH, idx, set(), 5, 25)[0] == "no_synteny"
