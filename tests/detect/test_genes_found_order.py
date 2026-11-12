# tests/detect/test_genes_found_order.py
"""`genes_found` lists optional slot genes in roster order, whatever the hash seed.

Two optional genes with different slots (GAP1 and BGL1, one slot each, in the Ascomycota
roster) were put in `genes_found` by iterating a set of slot names, so a genome's report
listed them in a different order from run to run (found 2026-11-12 on
GCA_028564845.1: seed 1 gave BGL1, GAP1 and seed 2 gave GAP1, BGL1). The order is output
only; no call changed. The check runs the scorer in fresh interpreters under several
`PYTHONHASHSEED` values, because the order of a string set is fixed per process.
"""
from __future__ import annotations

import os
import subprocess
import sys
import textwrap

_CODE = textwrap.dedent("""
    from MATPredict.detect.clustering import GeneCluster
    from MATPredict.detect.family_registry import Family, FamilyKey
    from MATPredict.detect.scoring import score_cluster
    from MATPredict.detect.search import SearchHit

    KEY = FamilyKey("Ascomycota", "MAT")
    GENES = [
        {"name": "MAT1-1-1", "role": "core_MAT"},
        {"name": "GAP1", "role": "flanking_conserved", "optional": True, "slot": "GAP1"},
        {"name": "BGL1", "role": "flanking_conserved", "optional": True, "slot": "BGL1"},
        {"name": "ZED1", "role": "flanking_conserved", "optional": True, "slot": "ZED1"},
    ]
    FAM = Family(KEY, "enum", ["MAT1-1"], None, GENES, [4751])
    roles = {g["name"]: g["role"] for g in GENES}
    names = ["MAT1-1-1", "BGL1", "GAP1", "ZED1"]
    hits = [SearchHit(KEY, g, roles[g], "c1", 1000 + 5000 * i, 3000 + 5000 * i, "+", 60.0, "rec1",
                      "diamond_proteome") for i, g in enumerate(names)]
    cl = GeneCluster("c1", hits[0].start, hits[-1].end, hits)
    print(",".join(score_cluster(cl, [FAM], searchable_genes={KEY: set(roles)})[0].genes_found))
""")


def _order(seed: str) -> str:
    env = dict(os.environ, PYTHONHASHSEED=seed)
    return subprocess.run([sys.executable, "-c", _CODE], env=env, check=True,
                          capture_output=True, text=True).stdout.strip()


def test_slot_genes_follow_the_roster_order_for_every_hash_seed():
    orders = {seed: _order(seed) for seed in ("0", "1", "2", "3", "4", "5")}
    assert set(orders.values()) == {"MAT1-1-1,GAP1,BGL1,ZED1"}, orders
