"""datasets/reads/sketch_store/starship_dynamics.py: Newick parsing, Fitch parsimony (against Sankoff), clone votes, thresholds, outputs."""
import csv
import importlib.util
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("starship_dynamics", ROOT / "datasets" / "reads" / "sketch_store" / "starship_dynamics.py")
sd = importlib.util.module_from_spec(spec); spec.loader.exec_module(sd)
A, P, E = sd.ABSENT, sd.PRESENT, sd.EITHER


def write(path, cols, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def tsv(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def test_parse_newick():
    ch, nm = sd.parse_newick("((a:0.1,b:-1e-05)x:0.2,(c,d),e);")
    leaves = sorted(nm[v] for v in range(len(ch)) if not ch[v])
    assert leaves == ["a", "b", "c", "d", "e"]
    assert len(ch[0]) == 3 and all(c > 0 for c in ch[0])           # trifurcating root
    assert any(nm[v] == "x" and len(ch[v]) == 2 for v in range(len(ch)))
    assert all(c > v for v in range(len(ch)) for c in ch[v])        # children have larger indices


def test_parse_deep_tree_does_not_recurse():
    n = 5000
    text = "(" * (n - 1) + "t0" + "".join(f",t{i})" for i in range(1, n)) + ";"
    ch, nm = sd.parse_newick(text)
    assert sum(1 for v in range(len(ch)) if not ch[v]) == n
    states = {f"t{i}": P if i % 2 else A for i in range(n)}
    score, rec = sd.fitch(ch, nm, states)
    assert score == sankoff(ch, nm, states) == n // 2 and sum(rec["root_absent"]) == score


def test_fitch_small():
    ch, nm = sd.parse_newick("((a,b),(c,d));")
    score, rec = sd.fitch(ch, nm, {"a": P, "b": P, "c": A, "d": A})
    assert score == 1 and rec["root_absent"] == (1, 0) and rec["root_present"] == (0, 1)
    score, _ = sd.fitch(ch, nm, {"a": P, "b": A, "c": P, "d": A})
    assert score == 2
    score, _ = sd.fitch(ch, nm, {"a": P, "c": P, "d": A})           # b is a wildcard
    assert score == 1
    assert sd.fitch(ch, nm, {})[0] == 0


def random_binary(rng, names):
    nodes = list(names)
    while len(nodes) > 1:
        i, j = rng.sample(range(len(nodes)), 2)
        a, b = nodes[i], nodes[j]
        nodes = [x for k, x in enumerate(nodes) if k not in (i, j)] + [f"({a},{b})"]
    return nodes[0] + ";"


def sankoff(ch, nm, states):
    """Minimum changes by dynamic programming (binary trees)."""
    INF = 10 ** 9
    cost = [None] * len(ch)
    for v in range(len(ch) - 1, -1, -1):
        if not ch[v]:
            m = states.get(nm[v], E)
            cost[v] = {s: (0 if m & b else INF) for s, b in ((A, A), (P, P))}
        else:
            cost[v] = {s: sum(min(cost[c][t] + (s != t) for t in (A, P)) for c in ch[v]) for s in (A, P)}
    return min(cost[0].values())


def test_fitch_matches_sankoff_and_counts_add_up():
    rng = random.Random(7)
    for _ in range(200):
        names = [f"t{i}" for i in range(rng.randint(2, 14))]
        ch, nm = sd.parse_newick(random_binary(rng, names))
        states = {n: rng.choice([A, P, E]) for n in names}
        score, rec = sd.fitch(ch, nm, states)
        assert score == sankoff(ch, nm, states)
        assert sum(rec["root_absent"]) == score and sum(rec["root_present"]) == score


def test_fitch_multifurcation_counts_add_up():
    rng = random.Random(3)
    for _ in range(100):
        names = [f"t{i}" for i in range(rng.randint(3, 12))]
        k = rng.randint(2, min(4, len(names)))
        groups = [names[i::k] for i in range(k)]
        text = "(" + ",".join("(" + ",".join(g) + ")" if len(g) > 1 else g[0] for g in groups) + ");"
        ch, nm = sd.parse_newick(text)
        score, rec = sd.fitch(ch, nm, {n: rng.choice([A, P]) for n in names})
        assert sum(rec["root_absent"]) == score and sum(rec["root_present"]) == score


def make_store(tmp_path):
    s = tmp_path / "store"; q = s / "query"
    # clones: c1 = s1+s2 (A and B 'partial-ish' 0.6), c2 = s3 (present), c3 = s4 (absent), c4 = s5 (absent), c5 = s8 (A 0.3 = missing);
    # s6 curator-excluded, s7 other species
    fr = {"s1": (1.0, 0.6), "s2": (1.0, 0.6), "s3": (1.0, 0.9), "s4": (0.0, 0.1), "s5": (0.0, 0.1), "s6": (1.0, 1.0), "s7": (1.0, 1.0), "s8": (0.3, 0.0)}
    write(q / "screen_starships_bafA" / "screen.tsv", ["sample_id", "target", "fraction"],
          [dict(sample_id=sid, target=t, fraction=f) for sid, v in fr.items() for t, f in zip(("A", "B"), v)])
    clone = dict(s1="c1", s2="c1", s3="c2", s4="c3", s5="c4", s6="c6", s7="c7", s8="c5")
    write(q / "distances" / "samples_qc.tsv", ["sample_id", "mean_ani", "clone_0.9995", "clone_0.998"],
          [dict(sample_id=sid, mean_ani="0.83" if sid == "s7" else "0.995", **{"clone_0.9995": c, "clone_0.998": "L1" if sid < "s5" else "L2"}) for sid, c in clone.items()])
    write(q / "mat" / "mat_calls.tsv", ["sample_id", "source", "study_accession", "country"],
          [dict(sample_id=sid, source="public", study_accession="P1", country="USA: CA" if sid < "s5" else "") for sid in clone])
    write(q / "qc" / "qc.tsv", ["sample_id", "flag"], [])
    write(s / "exclude.tsv", ["sample_id", "reason"], [dict(sample_id="s6", reason="low quality")])
    write(s / "manifest.tsv", ["sample_id", "strain"], [dict(sample_id=sid, strain=f"st{sid[1:]}") for sid in clone])
    (q / "distances" / "tree_fastme.nwk").write_text("((t1:0.1,t2:0.1):0.1,(t3:0.1,t4:0.1):0.1,(t5:0.1,t8:0.1):0.1);\n")
    write(q / "distances" / "tree_names.tsv", ["tip", "sample_id"], [dict(tip=f"t{i}", sample_id=f"s{i}") for i in (1, 2, 3, 4, 5, 8)])
    (q / "distances" / "tree_names.tsv").write_text("".join(f"t{i}\ts{i}\n" for i in (1, 2, 3, 4, 5, 8)))   # no header, as the build writes it
    write(tmp_path / "truth.tsv", ["sample_id", "A"], [dict(sample_id=sid, A=v) for sid, v in (("s1", 1), ("s2", 1), ("s3", 1), ("s4", 0), ("s5", 0), ("s8", 0))])
    return s


def test_end_to_end(tmp_path):
    out = tmp_path / "out"
    sd.main(["--store", str(make_store(tmp_path)), "--out", str(out), "--min-group", "1", "--truth", str(tmp_path / "truth.tsv")])
    cm = {r["clone"]: r for r in tsv(out / "clone_matrix.tsv")}
    assert set(cm) == {"c1", "c2", "c3", "c4", "c5"}                  # curator exclusion and other species dropped
    assert cm["c1"]["n_members"] == "2" and cm["c1"]["A"] == "1" and cm["c3"]["A"] == "0"
    assert cm["c5"]["A"] == "NA"                                       # 0.3 is partial -> missing, not loss
    prev = {r["element"]: r for r in tsv(out / "prevalence.tsv")}
    assert (prev["A"]["clones_typed"], prev["A"]["clones_present"]) == ("4", "2")
    assert prev["A"]["samples_present"] == "3" and prev["A"]["samples_missing"] == "1"
    par = {r["element"]: r for r in tsv(out / "parsimony.tsv")}
    assert par["A"]["fitch_score"] == "2" and par["A"]["clone_tips_present"] == "2" and par["A"]["clone_tips_typed"] == "4"
    for e in ("A", "B"):
        assert int(par[e]["gains_root_absent"]) + int(par[e]["losses_root_absent"]) == int(par[e]["fitch_score"])
    order = tsv(out / "tree_order.tsv")                                # left-to-right tree order of the clone tips
    assert [r["tip"] for r in order] == ["t1", "t3", "t4", "t5", "t8"] and [r["rank"] for r in order] == ["0", "2", "3", "4", "5"]
    co = [r for r in tsv(out / "cooccurrence.tsv") if (r["element_1"], r["element_2"]) == ("A", "B")]
    assert len(co) == 1 and co[0]["clones_typed"] == "4"
    sweep = {(r["element"], r["present_cut"]): r for r in tsv(out / "threshold_sweep.tsv")}
    r = sweep[("A", "0.5")]
    assert (r["TP"], r["FN"], r["FP"], r["TN"]) == ("3", "0", "0", "3")
    assert sweep[("A", "0.2")]["FP"] == "1"                            # s8 at 0.3 is a false positive at a 0.2 cut
    assert {r["group"] for r in tsv(out / "prevalence_by_lineage.tsv")} == {"L1", "L2"}
    assert "clone groups at ANI 0.9995: 5" in (out / "summary.txt").read_text()


def test_set_overrides_present_cut(tmp_path):
    out = tmp_path / "out"
    sd.main(["--store", str(make_store(tmp_path)), "--out", str(out), "--min-group", "1", "--set", "B=0.8"])
    cm = {r["clone"]: r for r in tsv(out / "clone_matrix.tsv")}
    assert cm["c1"]["B"] == "NA" and cm["c2"]["B"] == "1" and cm["c3"]["B"] == "0"   # 0.6 is now between the cuts
    assert cm["c1"]["A"] == "1"                                                      # other targets keep the default
    assert "B=0.8" in (out / "summary.txt").read_text()


def test_set_unknown_target_fails(tmp_path):
    try:
        sd.main(["--store", str(make_store(tmp_path)), "--out", str(tmp_path / "o"), "--set", "Nope=0.5"])
    except SystemExit as e:
        assert "Nope" in str(e)
    else:
        raise AssertionError("expected SystemExit")


def test_majority():
    assert sd.majority([True, True, False]) is True
    assert sd.majority([True, False]) is None
    assert sd.majority([None, None]) is None
    assert sd.majority([None, False]) is False
