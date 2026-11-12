"""datasets/reads/sketch_store/starship_report_md.py: markdown archive with PNG + PDF figures; local isolates are never named."""
import csv
import importlib.util
import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("starship_report_md", ROOT / "datasets" / "reads" / "sketch_store" / "starship_report_md.py")
md = importlib.util.module_from_spec(spec); spec.loader.exec_module(md)

ELS = ["bafA", "Nebuchadnezzar-h1", "Janus-h4", "Osiris-h3"]


def prev(e, k, n=100):
    return dict(element=e, clones=n, clones_typed=n, clones_present=k, frac_present=k / n, ci95_lo=max(0, k / n - 0.05), ci95_hi=k / n + 0.05,
                samples_present=k * 2, samples_missing=3, samples_absent=n - k)


def pars(e, k, score):
    return dict(element=e, clone_tips_typed=90, clone_tips_present=k, fitch_score=score, score_per_present=round(score / k, 2),
                gains_root_absent=score - 1, losses_root_absent=1, gains_root_present=score - 1, losses_root_present=1)


def family():
    gx = lambda: {str(k): dict(tp=5, fn=1, fp=0, tn=3) for k in (1, 2, 3, 4)}
    v = lambda: dict(taxa=10, matched_taxa=9, pairs=9, genes={g: gx() for g in ("bafA", "bafB", "bafC", "Pseudobaf")},
                     previous_single_target_bafA=dict(tp=4, fn=2, fp=0, tn=3), bafA_rel_abund_by_copy_number={})
    cuts = {"bafA": 2, "bafB": 1, "bafC": 2, "Pseudobaf": 2}
    return dict(screen={}, cuts=cuts, validation={"older.tab": v(), "update.tab": v()},
                population=dict(kept=110, genes={g: dict(cut=c, isolates=20 + i, clone_groups=9 + i, clone_groups_total=100) for i, (g, c) in enumerate(cuts.items())},
                                old_vs_new_bafA=dict(both=4, new_only=3, old_only=1, neither=102), bafA_all_vs_neb_h1=dict(both=4, bafA_only=3, neb_only=0, neither=103),
                                newly_called=dict(isolates=3, neb_h1_present=0),
                                classes=dict(cluster_linked=dict(isolates=5, clone_groups=4, clone_groups_total=100, neb_h1_present=4, neb_h1_absent=1, clone_groups_any_member=5, largest_clone_group=3, isolates_in_groups_over_1=4),
                                             standalone_divergent=dict(isolates=3, clone_groups=3, clone_groups_total=100, neb_h1_present=0, neb_h1_absent=3, clone_groups_any_member=3, largest_clone_group=1, isolates_in_groups_over_1=0)),
                                placement=dict(standalone_divergent=dict(genomes=2, bafA_hit=2, hrmA_or_cgnA_in_genome=0), cluster_linked=dict(genomes=4, bafA_hit=4, bafA_within_1_1kb_of_hrmA_cgnA=4)),
                                bafA_hits_histogram={"2": 5}, by_source=dict(public=dict(present=5, total=90), local=dict(present=2, total=20))))


def data(with_family=True):
    rng = random.Random(1)
    ks = dict(zip(ELS, (5, 4, 30, 25)))
    summ = ["samples in the screen: 120; kept 110; excluded {'other_species': 9, 'curator_exclude': 1}",
            "clone groups at ANI 0.9995: 100; present cut 0.5, absent cut 0.2; overrides: Osiris-h3=0.7",
            "tree tips: 100; clone groups with a tip in the tree: 90 of 100"]
    grp = lambda key, names: [dict(group=g, element=e, clones=40, clones_typed=40, clones_present=k, frac_present=k / 40, ci95_lo=0.0, ci95_hi=0.3) for g, k in zip(names, (4, 1)) for e in ELS]
    side = lambda: dict(prevalence=[prev(e, ks[e]) for e in ELS], parsimony=[pars(e, ks[e], ks[e] * 2 // 3) for e in ELS],
                        by_lineage=grp("l", ["150", "134"]), by_country=grp("c", ["USA", "China"]), summary=list(summ))
    heat = [dict(clone=f"c{i}", n=1, lineage=rng.choice(["150", "134", "78", "none"]), country="USA", in_tree=i < 90,
                 bits="".join(rng.choice("0012") for _ in ELS)) for i in range(100)]
    sweep = {e: [dict(present_cut=round(0.1 * c, 1), TP=5, FN=0, FP=max(0, 10 - c), TN=50, precision=round(5 / (5 + max(0, 10 - c)), 3), recall=1.0) for c in range(1, 11)] for e in ELS}
    cooc = [dict(element_1=a, element_2=b, clones_typed=95, both=3, only_1=2, only_2=1, neither=89, phi=0.5 if (a, b) == ("bafA", "Nebuchadnezzar-h1") else 0.1, jaccard=0.5)
            for i, a in enumerate(ELS) for b in ELS[i + 1:]]
    cand = [dict(tier="strong", sample_id="SAMN1", strain="AF1", source="public", study="PRJNA1", country="USA", bafA_frac=1.0, bafA_call="present", neb_h1_frac=0.9, neb_h1_call="present", excluded="", clone=1, lineage=150),
            dict(tier="strong", sample_id="local__SECRET-1", strain="SECRET-1", source="local", study="L", country="local", bafA_frac=1.0, bafA_call="present", neb_h1_frac=0.9, neb_h1_call="present", excluded="", clone=2, lineage=150),
            dict(tier="probable", sample_id="SAMN2", strain="AF2", source="public", study="PRJNA1", country="USA", bafA_frac=1.0, bafA_call="present", neb_h1_frac=0.6, neb_h1_call="partial", excluded="", clone=3, lineage=134)]
    tr = lambda e, fp: dict(len=100, specific=10, present=20, partial=5, absent=75, truth=dict(TP=10, FN=1, FP=fp, TN=200), dropped="")
    before = {e: tr(e, 27) for e in ELS}; after = {e: dict(tr(e, 27), partial=0) for e in ELS}
    return dict(meta=dict(code_commit="abc1234", branch="b", pr="1", rescreen_job="1", default_job="2", tuned_job="3", store="/s", results="/r", date="2026-11-11"),
                elements=ELS, pair=["bafA", "Nebuchadnezzar-h1"], default=side(), tuned=side(), heat=heat, baf_neb={"11": 4, "10": 1, "01": 0, "00": 85},
                sweep=sweep, cooccurrence=cooc, candidates=cand, baf_family=family() if with_family else None, rescreen=dict(before=before, after=after, dropped=[dict(target="bafA", hash=1, samples_alone=778)]))


def test_markdown_report(tmp_path):
    (tmp_path / "d.json").write_text(json.dumps(data()))
    out = tmp_path / "rep"
    md.main(["--data", str(tmp_path / "d.json"), "--out", str(out)])
    text = (out / "report.md").read_text()
    for name in ("fig1_prevalence", "fig2_threshold_sweep", "fig3_tree_heatmap", "fig4_gain_loss", "fig5_bafA_groups", "fig6_cooccurrence", "fig7_baf_family"):
        assert (out / "figures" / f"{name}.png").stat().st_size > 5000 and (out / "figures" / f"{name}.pdf").read_bytes().startswith(b"%PDF")
        assert f"[![" in text and f"(figures/{name}.png)](figures/{name}.pdf)" in text     # PNG embedded, linked to its PDF
    for section in ("## Summary", "## 3. How common", "## 8. bafA beyond the AF293-type allele", "## 9. Which elements", "## 10. The strains", "## 13. Traceability"):
        assert section in text
    assert "`abc1234`" in text and "bafA" in text and "2026-11-11" in text
    assert "AF293-type bafA" in text and "counting divergent alleles" in text and (out / "tables" / "baf_family_summary.json").exists()
    assert "Two classes of bafA" in text and "| standalone divergent (all-alleles screen only) | 3 | 3 of 100 (3.0%) | 0 | 3 |" in text and "not the same as** Nebuchadnezzar-h1" in text
    assert "Repeats of the same strain" in text and "the largest holds 3 isolates and 4 isolates sit" in text
    assert "have hrmA or cgnA anywhere in the genome" in text
    assert "| older.tab | bafA | 2 | 5 / 6 | 0 | 83% |" in text                 # validation table row from the family summary
    # local isolates counted, never named
    cand = (out / "tables" / "candidates_public.tsv").read_text()
    assert "SAMN1" in cand and "SAMN2" in cand and "SECRET" not in cand and "local__" not in cand
    assert "SECRET" not in text and "local__" not in text
    rows = list(csv.DictReader(open(out / "tables" / "prevalence.tsv"), delimiter="\t"))
    assert {r["element"] for r in rows} == set(ELS)
    tiers = [ln for ln in text.splitlines() if ln.startswith("| strong |")]
    cells = [c.strip() for c in tiers[0].split("|")[1:-1]]
    assert cells == ["strong", "2", "2", "1", "1"]                       # tier, all, after flagged removed, public 1, local 1 (counted only)


def test_report_without_family(tmp_path):
    (tmp_path / "d.json").write_text(json.dumps(data(with_family=False)))
    out = tmp_path / "rep"
    md.main(["--data", str(tmp_path / "d.json"), "--out", str(out)])
    text = (out / "report.md").read_text()
    assert "## 8. bafA beyond" not in text and "## 9. Which elements" in text and not (out / "figures" / "fig7_baf_family.png").exists()
    assert "counting divergent alleles" not in text and "Two classes of bafA" not in text
