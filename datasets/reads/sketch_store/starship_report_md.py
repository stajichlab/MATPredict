#!/usr/bin/env python3
"""Markdown archive of the Starship report: report.md, figures (PNG linked to PDF) and the tables behind them, for GitHub.

    starship_report_md.py --data report_data.json --out DIR [--date YYYY-MM-DD]

Reads the JSON written by starship_report_data.py (the same file that fills the HTML page), so the markdown and the page show the same
numbers. Writes DIR/report.md, DIR/figures/<name>.png and <name>.pdf (the markdown embeds the PNG and links the PDF), and DIR/tables/*.tsv.

Local (Population_Genomics) isolates count in every number but are never listed by name, as in the site report (export_report.py):
tables/candidates_public.tsv holds public isolates only; local candidates appear as counts. Figures carry no sample identifiers."""
import argparse, collections, csv, json, math, os, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

C = dict(s1="#2a78d6", s2="#eb6834", s3="#1baf7a", neg="#e34948", ink="#0b0b0b", ink2="#52514e", mute="#6f6e69", grid="#dcdcd7", soft="#eceae6")
plt.rcParams.update({"font.size": 9, "axes.edgecolor": C["grid"], "axes.labelcolor": C["ink2"], "xtick.color": C["ink2"], "ytick.color": C["ink2"],
                     "text.color": C["ink"], "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "axes.grid": False, "figure.dpi": 100})
KEY = ("bafA", "Nebuchadnezzar-h1")


def pc(x, d=1):
    return f"{100 * float(x):.{d}f}%"


def fam(e):
    return e.split("-")[0]


def order(elements):
    return sorted(elements, key=lambda e: (e != "bafA", fam(e), e))


def save(fig, out, name):
    os.makedirs(f"{out}/figures", exist_ok=True)
    fig.savefig(f"{out}/figures/{name}.png", dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{out}/figures/{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def bold_keys(ax, axis="y"):
    for t in (ax.get_yticklabels() if axis == "y" else ax.get_xticklabels()):
        if t.get_text() in KEY:
            t.set_fontweight("bold")


class R:
    """Derived values shared by the figures and the text."""
    def __init__(s, D):
        s.D, s.E = D, D["elements"]
        s.P = {r["element"]: r for r in D["default"]["prevalence"]}
        s.PS = {r["element"]: r for r in D["default"]["parsimony"]}
        s.PT = {r["element"]: r for r in D["tuned"]["prevalence"]} if D.get("tuned") else None
        s.PST = {r["element"]: r for r in D["tuned"]["parsimony"]} if D.get("tuned") else None
        s.before, s.after = D["rescreen"]["before"], D["rescreen"]["after"]
        s.m = D["meta"]

    def truth(s, e):
        t = (s.after.get(e) or {}).get("truth")
        return t

    def over(s, e):
        t = s.truth(e)
        return bool(t and t["FP"] >= 9)

    def cut(s, e):
        if not s.D.get("tuned"):
            return None
        o = re.search(re.escape(e) + r"=([\d.]+)", s.D["tuned"]["summary"][1])
        return float(o.group(1)) if o else None

    def sweep_at(s, e, c):
        return next((r for r in s.D["sweep"].get(e, []) if abs(r["present_cut"] - c) < 1e-9), None)


def fig_prevalence(r, out):
    rows = sorted(r.E, key=lambda e: -r.P[e]["frac_present"])
    fig, ax = plt.subplots(figsize=(7.6, 0.27 * len(rows) + 1.2))
    for i, e in enumerate(rows):
        p = r.P[e]; col = C["s2"] if r.over(e) else C["s1"]; y = len(rows) - 1 - i
        ax.plot([p["ci95_lo"], p["ci95_hi"]], [y, y], color=col, lw=2, solid_capstyle="butt")
        ax.plot(p["frac_present"], y, "o", color=col, ms=6, mec="white", mew=1.2)
        ax.text(p["ci95_hi"] + 0.008, y, f"{pc(p['frac_present'])} ({p['clones_present']}/{p['clones_typed']:,})", va="center", fontsize=8, color=C["mute"])
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows[::-1]); bold_keys(ax)
    ax.set_xlim(0, 0.75); ax.xaxis.set_major_formatter(lambda v, _: f"{100 * v:.0f}%"); ax.grid(axis="x", color=C["soft"]); ax.set_axisbelow(True)
    ax.set_xlabel("share of clone groups carrying the element (line = 95% Wilson interval)")
    ax.legend(handles=[plt.Line2D([], [], marker="o", color=C["s1"], ls="", label="agrees with the paper's calls (8 or fewer extra calls)"),
                       plt.Line2D([], [], marker="o", color=C["s2"], ls="", label="9 or more extra calls vs the paper: upper bound")], loc="lower right", frameon=False, fontsize=8)
    save(fig, out, "fig1_prevalence")


def fig_sweep(r, out):
    els = [e for e in ["Osiris-h3", "Janus-h4", "Janus-h2", "Janus-h1", "Janus-h3", "Janus-h5", "bafA", "Nebuchadnezzar-h1"] if e in r.D["sweep"]]
    fig, axs = plt.subplots(2, 4, figsize=(10, 4.6), sharey=True)
    for ax, e in zip(axs.flat, els):
        rows = r.D["sweep"][e]; xs = [x["present_cut"] for x in rows]
        ax.plot(xs, [x["precision"] if x["precision"] != "" else float("nan") for x in rows], color=C["s1"], lw=2, label="precision")
        ax.plot(xs, [x["recall"] if x["recall"] != "" else float("nan") for x in rows], color=C["s2"], lw=2, label="recall")
        ax.axvline(0.5, color=C["mute"], ls="--", lw=1)
        c = r.cut(e)
        if c:
            ax.axvline(c, color=C["ink"], lw=1.5); ax.text(c + 0.02, 0.06, f"tuned {c}", fontsize=7, color=C["ink2"])
        ax.set_title(e, fontsize=9, loc="left", fontweight="bold"); ax.set_ylim(0, 1.03); ax.set_xlim(0.1, 1.0); ax.set_xticks([0.1, 0.5, 0.9])
    for ax in list(axs.flat)[len(els):]:
        ax.axis("off")
    for ax in axs[-1]:
        ax.set_xlabel("present cut (hash fraction)")
    axs[0][0].legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    save(fig, out, "fig2_threshold_sweep")


def tree_bins(r, BIN=6):
    rows = [h for h in r.D["heat"] if h["in_tree"]]
    idx = {e: i for i, e in enumerate(r.E)}; ord_ = order(r.E)
    lin = collections.Counter(h["lineage"] for h in rows if h["lineage"] != "none")
    top3 = [l for l, _ in lin.most_common(3)]
    nb = math.ceil(len(rows) / BIN); mat = [[float("nan")] * nb for _ in ord_]; dom = []
    for b in range(nb):
        seg = rows[b * BIN:(b + 1) * BIN]
        dom.append(collections.Counter(h["lineage"] for h in seg).most_common(1)[0][0])
        for k, e in enumerate(ord_):
            v = [h["bits"][idx[e]] for h in seg if h["bits"][idx[e]] != "2"]
            if v:
                mat[k][b] = sum(x == "1" for x in v) / len(v)
    return rows, ord_, mat, dom, top3, nb


def fig_tree(r, out):
    import numpy as np
    rows, ord_, mat, dom, top3, nb = tree_bins(r)
    cmap = LinearSegmentedColormap.from_list("seq", ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]); cmap.set_bad(C["soft"])
    fig, (a0, a1) = plt.subplots(2, 1, figsize=(10, 6.2), gridspec_kw=dict(height_ratios=[1, 20], hspace=0.04), sharex=True)
    cols = [C["s1"], C["s2"], C["s3"]]
    strip = np.array([[matplotlib.colors.to_rgb(cols[top3.index(d)] if d in top3 else C["grid"]) for d in dom]])
    a0.imshow(strip, aspect="auto", interpolation="nearest"); a0.set_yticks([0]); a0.set_yticklabels(["lineage"], fontsize=8); a0.tick_params(length=0); [s.set_visible(False) for s in a0.spines.values()]
    a1.imshow(np.array(mat, dtype=float), aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    a1.set_yticks(range(len(ord_))); a1.set_yticklabels(ord_); bold_keys(a1)
    for k in range(1, len(ord_)):
        if fam(ord_[k - 1]) != fam(ord_[k]):
            a1.axhline(k - 0.5, color="white", lw=1.5)
    a1.set_xticks([0, nb - 1]); a1.set_xticklabels(["one end of the tree", "other end"]); a1.tick_params(length=0)
    handles = [plt.Line2D([], [], marker="s", ls="", color=cols[i], label=f"lineage {l}") for i, l in enumerate(top3)] + [plt.Line2D([], [], marker="s", ls="", color=C["grid"], label="other lineages")]
    a1.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=4, frameon=False, fontsize=8)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1)); cb = fig.colorbar(sm, ax=[a0, a1], fraction=0.015, pad=0.01); cb.set_label("share of clone groups in the bin")
    save(fig, out, "fig3_tree_heatmap")


def fig_gainloss(r, out):
    rows = r.D["default"]["parsimony"]; fig, ax = plt.subplots(figsize=(7.4, 5))
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(8, 2000); ax.set_ylim(3, 400)
    for k, ls, lab in ((1, "-", "ratio 1: every clone group its own event"), (0.25, "--", "0.25"), (0.1, "--", "0.1")):
        ax.plot([8, 2000], [8 * k, 2000 * k], color=C["mute"], ls=ls, lw=1)
    ax.text(40, 40 * 1.0 * 1.25, "ratio 1: each present clone group its own event", fontsize=7.5, color=C["mute"], rotation=29, ha="center")
    ax.text(1200, 1200 * 0.25 * 1.15, "0.25", fontsize=8, color=C["mute"]); ax.text(1200, 1200 * 0.1 * 1.15, "0.1", fontsize=8, color=C["mute"])
    off = {"Janus-h1": (-6, 6, "right"), "Nebuchadnezzar-h1": (-8, -11, "right"), "bafA": (8, 4, "left")}
    for p in rows:
        e = p["element"]; key = e in KEY; x, y = p["clone_tips_present"], max(p["fitch_score"], 3)
        ax.plot(x, y, "o", color=C["s1"], ms=9 if key else 6, alpha=1 if key else 0.8, mec="white", mew=1.2)
        if e in KEY or e in ("Lamia-h4", "Osiris-h4", "Osiris-h3", "Gnosis-h1", "Lamia-h3", "Janus-h1"):
            dx, dy, ha = off.get(e, (7, -3, "left"))
            ax.annotate(e, (x, y), xytext=(dx, dy), textcoords="offset points", ha=ha, fontsize=8, fontweight="bold" if key else "normal", color=C["ink"] if key else C["ink2"])
    ax.set_xlabel("clone groups carrying the element (log scale)"); ax.set_ylabel("minimum changes on the tree (log scale)")
    ax.grid(color=C["soft"]); ax.set_axisbelow(True)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:,.0f}"); ax.yaxis.set_major_formatter(lambda v, _: f"{v:,.0f}")
    save(fig, out, "fig4_gain_loss")


def fig_groups(r, out):
    fig, axs = plt.subplots(1, 2, figsize=(10, 5.2))
    baf = r.P["bafA"]
    for ax, key, title in ((axs[0], "by_country", "By country"), (axs[1], "by_lineage", "By lineage")):
        rows = sorted([x for x in r.D["default"][key] if x["element"] == "bafA" and x["clones_typed"] >= 10], key=lambda x: -x["frac_present"])
        for i, x in enumerate(rows):
            y = len(rows) - 1 - i
            ax.plot([max(0, x["ci95_lo"]), min(x["ci95_hi"], 0.45)], [y, y], color=C["s1"], lw=2)
            ax.plot(min(x["frac_present"], 0.45), y, "o", color=C["s1"], ms=5.5, mec="white", mew=1)
            ax.text(min(x["ci95_hi"], 0.45) + 0.008, y, f"{pc(x['frac_present'], 0)} ({x['clones_present']}/{x['clones_typed']})", va="center", fontsize=7.5, color=C["mute"])
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([str(x["group"]) for x in rows[::-1]], fontsize=8)
        ax.axvline(baf["frac_present"], color=C["mute"], ls="--", lw=1, zorder=0); ax.set_xlim(0, 0.55); ax.xaxis.set_major_formatter(lambda v, _: f"{100 * v:.0f}%")
        ax.set_title(f"{title} (10 or more clone groups)", loc="left", fontsize=9, fontweight="bold"); ax.grid(axis="x", color=C["soft"]); ax.set_axisbelow(True)
    axs[0].text(baf["frac_present"] + 0.004, len(axs[0].get_yticks()) - 0.35, f"overall {pc(baf['frac_present'])}", fontsize=7.5, color=C["mute"])
    fig.supxlabel("clone groups carrying the AF293-type bafA (95% Wilson interval; dashed = overall)", fontsize=9)
    fig.tight_layout()
    save(fig, out, "fig5_bafA_groups")


def fig_cooc(r, out):
    ord_ = order(r.E); n = len(ord_); phi = {}
    for x in r.D["cooccurrence"]:
        if x["phi"] != "":
            phi[(x["element_1"], x["element_2"])] = phi[(x["element_2"], x["element_1"])] = float(x["phi"])
    import numpy as np
    M = np.full((n, n), np.nan)
    for i, a in enumerate(ord_):
        for j, b in enumerate(ord_):
            if (a, b) in phi:
                M[i, j] = phi[(a, b)]
    cmap = LinearSegmentedColormap.from_list("div", [C["neg"], "#f0efec", C["s1"]]); cmap.set_bad(C["soft"])
    fig, ax = plt.subplots(figsize=(7.4, 6.6)); im = ax.imshow(M, cmap=cmap, vmin=-1, vmax=1)
    ax.set_xticks(range(n)); ax.set_xticklabels(ord_, rotation=60, ha="left", fontsize=8); ax.xaxis.tick_top(); ax.set_yticks(range(n)); ax.set_yticklabels(ord_, fontsize=8)
    bold_keys(ax); bold_keys(ax, "x")
    for i in range(n):
        for j in range(n):
            if not math.isnan(M[i, j]) and abs(M[i, j]) >= 0.5:
                ax.text(j, i, f"{M[i, j]:.2f}".replace("0.", "."), ha="center", va="center", fontsize=7, color="white" if abs(M[i, j]) >= 0.7 else C["ink"])
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("phi (presence correlation across clone groups)")
    for s in ax.spines.values():
        s.set_visible(False)
    save(fig, out, "fig6_cooccurrence")


def fig_baf_family(r, out):
    fam = r.D["baf_family"]; pop = fam["population"]; tname = list(fam["validation"])[-1]; v = fam["validation"][tname]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 3.9), gridspec_kw=dict(width_ratios=[1.1, 1]))
    p = r.P["bafA"]
    rows = [("bafA, AF293-type\n(Starship screen target)", p["clones_present"] / p["clones_typed"], C["s2"])]
    rows += [(lab, pop["genes"][g]["clone_groups"] / pop["genes"][g]["clone_groups_total"], C["s1"]) for g, lab in
             (("bafA", "bafA, all alleles"), ("bafB", "bafB"), ("bafC", "bafC"), ("Pseudobaf", "Pseudobaf"))]
    for i, (lab, x, col) in enumerate(rows):
        y = len(rows) - 1 - i
        a.barh(y, x, color=col, height=0.6); a.text(x + 0.01, y, pc(x), va="center", fontsize=8.5, color=C["ink2"])
    a.set_yticks(range(len(rows))); a.set_yticklabels([x[0] for x in rows][::-1], fontsize=8.5); a.set_xlim(0, 0.6)
    a.xaxis.set_major_formatter(lambda t_, _: f"{100 * t_:.0f}%"); a.set_xlabel("share of clone groups carrying the gene"); a.grid(axis="x", color=C["soft"]); a.set_axisbelow(True)
    genes = [("bafA", "bafA (bafD column)"), ("bafB", "bafB"), ("bafC", "bafC"), ("Pseudobaf", "Pseudobaf")]
    for i, (g, lab) in enumerate(genes):
        y = len(genes) - 1 - i; x = v["genes"][g][str(fam["cuts"][g])]; n = x["tp"] + x["fn"]
        b.plot(x["tp"] / n, y, "o", color=C["s1"], ms=8, mec="white", mew=1.2)
        b.text(min(x["tp"] / n, 0.9) + 0.03, y, f"{x['tp']}/{n} found, {x['fp']} extra", va="center", fontsize=8.5, color=C["ink2"])
    b.set_yticks(range(len(genes))); b.set_yticklabels([x[1] for x in genes][::-1], fontsize=8.5); b.set_xlim(0.5, 1.45); b.set_xticks([0.5, 0.75, 1.0])
    b.xaxis.set_major_formatter(lambda t_, _: f"{100 * t_:.0f}%"); b.set_xlabel(f"assembly-positive strains found by the screen ({tname.replace('.tab', '')})"); b.grid(axis="x", color=C["soft"]); b.set_axisbelow(True)
    fig.tight_layout()
    save(fig, out, "fig7_baf_family")


def md_table(head, rows, align=None):
    align = align or ["l"] + ["r"] * (len(head) - 1)
    sep = ["---:" if a == "r" else ":---" for a in align]
    return "\n".join(["| " + " | ".join(head) + " |", "| " + " | ".join(sep) + " |"] + ["| " + " | ".join(str(c) for c in row) + " |" for row in rows])


def tsv_out(path, rows, cols=None):
    cols = cols or (list(rows[0]) if rows else ["empty"])
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


def fig_md(name, alt):
    return f"[![{alt}](figures/{name}.png)](figures/{name}.pdf)\n\n*Figure: {alt}. Click the image for the vector PDF.*"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", required=True); ap.add_argument("--out", required=True); ap.add_argument("--date", default="")
    a = ap.parse_args(argv)
    D = json.load(open(a.data)); r = R(D); out = a.out; M = r.m
    os.makedirs(out, exist_ok=True)
    for f in (fig_prevalence, fig_sweep, fig_tree, fig_gainloss, fig_groups, fig_cooc) + ((fig_baf_family,) if D.get("baf_family") else ()):
        f(r, out)
    # ---- numbers
    m1 = re.search(r"samples in the screen: (\d+); kept (\d+); excluded (\{.*\})", D["default"]["summary"][0]); m2 = re.search(r"clone groups at ANI ([\d.]+): (\d+)", D["default"]["summary"][1])
    m3 = re.search(r"clone groups with a tip in the tree: (\d+) of (\d+)", " ".join(D["default"]["summary"]))
    nScreen, nKept, excl = int(m1.group(1)), int(m1.group(2)), m1.group(3).replace("'", "").strip("{}").replace("_", " ")
    ani, nClones, nTree = m2.group(1), int(m2.group(2)), int(m3.group(1)) if m3 else 0
    fp_ = (D.get("baf_family") or {}).get("population")
    baf, neb, bp, nbp = r.P["bafA"], r.P["Nebuchadnezzar-h1"], r.PS["bafA"], r.PS["Nebuchadnezzar-h1"]
    bn = D["baf_neb"] or {"11": 0, "10": 0, "01": 0, "00": 0}; cands = D["candidates"]; date = a.date or M.get("date", "")
    pub = [c for c in cands if c.get("source") != "local"]; loc = [c for c in cands if c.get("source") == "local"]
    tier = lambda rows, t, kept=False: [c for c in rows if c["tier"] == t and (not kept or not c["excluded"])]
    strong_kept = tier(cands, "strong", True)
    clones = len({c["clone"] or "solo:" + c["sample_id"] for c in strong_kept})
    # ---- tables on disk
    t = f"{out}/tables"
    tsv_out(f"{t}/prevalence.tsv", D["default"]["prevalence"]); tsv_out(f"{t}/parsimony.tsv", D["default"]["parsimony"])
    tsv_out(f"{t}/prevalence_by_country.tsv", D["default"]["by_country"]); tsv_out(f"{t}/prevalence_by_lineage.tsv", D["default"]["by_lineage"])
    tsv_out(f"{t}/cooccurrence.tsv", D["cooccurrence"])
    tsv_out(f"{t}/threshold_sweep.tsv", [dict(element=e, **x) for e, v in D["sweep"].items() for x in v])
    tsv_out(f"{t}/rescreen_before_after.tsv", [dict(element=e, specific_before=r.before[e]["specific"], specific_after=r.after[e]["specific"], present_before=r.before[e]["present"],
            present_after=r.after[e]["present"], partial_before=r.before[e]["partial"], partial_after=r.after[e]["partial"],
            **{f"{k}_{s}": ((src.get(e) or {}).get("truth") or {}).get(k, "") for s, src in (("before", r.before), ("after", r.after)) for k in ("TP", "FN", "FP")}) for e in r.E])
    if D.get("baf_family"):
        os.makedirs(t, exist_ok=True)
        json.dump(D["baf_family"], open(f"{t}/baf_family_summary.json", "w"), indent=1)
    pcols = ["tier", "sample_id", "strain", "study", "country", "bafA_frac", "bafA_call", "neb_h1_frac", "neb_h1_call", "excluded", "clone", "lineage"]
    tsv_out(f"{t}/candidates_public.tsv", pub, pcols)
    # ---- markdown
    L = []
    w = L.append
    w(f"# Starships and bafA across {nScreen:,} *Aspergillus fumigatus* isolates\n")
    w(f"Sketch-store analysis, {date}. Which isolates carry a Nebuchadnezzar-h1-type Starship with bafA, how common the other 20 defined Starship haplotypes are, "
      "and how their presence is spread over the population tree. Every figure and table names its source. An interactive version of this report exists as a private page; "
      "this file is the archived copy.\n")
    w(md_table(["Provenance", ""], [["code", f"`{M.get('code_commit', '?')}` on branch `{M.get('branch', '?')}` (PR #{M.get('pr', '?')})"], ["rescreen job", f"`{M.get('rescreen_job', '?')}`"],
      ["default run", f"`{M.get('default_job', '?')}`"], ["tuned run", f"`{M.get('tuned_job', '?')}`"], ["store (HPCC)", f"`{M.get('store', '')}`"], ["results (HPCC)", f"`{M.get('results', '')}`"]], ["l", "l"]) + "\n")
    w("## Summary\n")
    w(md_table(["", ""], [[f"{nKept:,}", f"isolates kept of {nScreen:,} sketched"], [f"{nClones:,}", f"clone groups (ANI {ani}), one vote each"],
      [pc(baf["frac_present"]), f"of clone groups carry the AF293-type bafA ({baf['clones_present']} of {baf['clones_typed']:,})" + (f"; {pc(fp_['genes']['bafA']['clone_groups'] / fp_['genes']['bafA']['clone_groups_total'])} counting divergent alleles (section 8)" if fp_ else "")],
      [f"{len(strong_kept):,}", f"isolates in the strong bafA + Nebuchadnezzar-h1 list ({clones} clone groups)"]], ["r", "l"]) + "\n")
    w(f"- **Solid.** The AF293-type bafA is present in {pc(baf['frac_present'])} of clone groups and almost always sits with Nebuchadnezzar-h1 ({bn['11']} of {bn['11'] + bn['10'] + bn['01']} clone groups with either one carry both).")
    w(f"- **Solid.** The solo-hash rule removed {r.before['bafA']['partial']:,} ambiguous partial bafA calls without changing agreement with the paper ({r.after['bafA']['truth']['TP']} true positives, {r.after['bafA']['truth']['FN']} missed).")
    if fp_:
        bfv = list(D["baf_family"]["validation"].values())[-1]["genes"]["bafA"]["2"]
        w(f"- **Solid.** A screen of the whole baf family (one target per gene, built from assembly haplotypes) matches the assembly counts: bafA {bfv['tp']} of {bfv['tp'] + bfv['fn']} assembly-positive strains found with {bfv['fp']} extra calls (the single AF293-type target found 41 of 55). Counting divergent alleles, bafA is in {fp_['genes']['bafA']['clone_groups']} of {fp_['genes']['bafA']['clone_groups_total']:,} clone groups; {fp_['newly_called']['isolates']} isolates are newly called and {fp_['newly_called']['neb_h1_present']} of them have Nebuchadnezzar-h1 present.")
    w(f"- **Caution.** AF293-type bafA presence is scattered over the tree: {bp['fitch_score']} changes for {bp['clone_tips_present']} present clone groups (ratio {bp['score_per_present']}); Nebuchadnezzar-h1 behaves the same ({nbp['score_per_present']}). That fits repeated gain or transfer, or a coarse tree; it does not by itself show either.")
    w(f"- **Caution.** Janus-h2, Janus-h4 and Osiris-h3 call more isolates than the paper does. Raising the cut for Osiris-h3 (0.5 to {r.cut('Osiris-h3') or 0.7}) removes most of the extra calls; Janus-h3, Janus-h5 and the AF293-type bafA keep theirs at every cut, with hash fraction near 1.0.")
    w("- **Open.** The cuts were picked on the same 226 truth samples (no hold-out yet). Only the 20 haplotypes of the NII pangenome study are screened; the Starbase / Gluck-Thaler set is not.\n")
    w("## 1. Data and method\n")
    w(md_table(["Step", "What", "Where"], [
      ["Read sketches", "k21 and k31, scaled 100, abundance kept (>= 2), 30x subsample per isolate", "`sketch_store/aspergillus_fumigatus/sigs/`"],
      ["Clone groups", f"average-linkage groups at ANI {ani}: {nClones:,} groups", "`query/distances/samples_qc.tsv`"],
      ["Screen", "fraction of each target's specific hashes found, 21 targets (bafA + 20 haplotypes)", "`query/screen_starships_bafA/screen.tsv`"],
      ["Solo-hash rule", f"{len(D['rescreen']['dropped'])} hashes dropped that appear alone in > 1% of isolates", f"`dropped_hashes.tsv`, job `{M.get('rescreen_job', '?')}`"],
      ["Call", "present >= 0.5, absent < 0.2, otherwise missing; clone vote = majority of its members", "`starship_dynamics.py`"],
      ["Tree", f"Fitch parsimony on the FastME distance tree, one tip per clone group ({nTree:,} of {nClones:,} have a tip)", "`query/distances/tree_fastme.nwk`"]], ["l", "l", "l"]) + "\n")
    w(f"Left out of every frequency: {excl} isolates (other species, k-mer mixtures, one curator exclusion). Partial hits (0.2 to 0.5) count as missing, never as loss.\n")
    w("## 2. Effect of the solo-hash rule\n")
    w("A hash seen alone in many isolates is shared with other sequence. \"Paper\" columns use the 226-sample truth table at a present cut of 0.5.\n")
    ar = lambda x, y: f"{x:,}" if x == y else f"{x:,} -> {y:,}"
    w(md_table(["Element", "Specific hashes", "Present", "Partial", "Paper TP", "Paper FP", "Paper FN"], [[("**%s**" % e if e in KEY else e), ar(r.before[e]["specific"], r.after[e]["specific"]), ar(r.before[e]["present"], r.after[e]["present"]),
      ar(r.before[e]["partial"], r.after[e]["partial"])] + ([ar(r.before[e]["truth"][k], r.after[e]["truth"][k]) for k in ("TP", "FP", "FN")] if r.after[e]["truth"] else ["", "", ""]) for e in r.E]) + "\n")
    w(f"The rule only helps where the extra hash was a single non-specific one: bafA partials fall from {r.before['bafA']['partial']:,} to {r.after['bafA']['partial']:,}. False positives for Janus and Osiris-h3 do not fall, because their overcall comes from hashes shared across the family.\n")
    w(f"Source: `tables/rescreen_before_after.tsv`; `query_2026-11-09_n4834/screen_starships_bafA_before_solo_rule/summary.txt` and `query/screen_starships_bafA/summary.txt`, job `{M.get('rescreen_job', '?')}`.\n")
    w("## 3. How common is each element\n")
    w(fig_md("fig1_prevalence", "Share of clone groups carrying each element") + "\n")
    w("One point is one clone group, so a study that sequenced one lineage fifty times counts once. Orange rows are upper bounds only. In this and the following sections `bafA` is the AF293-type allele (the 2-hash genomic target of the Starship screen); section 8 counts all bafA alleles.\n")
    w(md_table(["Element", "Clone groups present / typed", "Share", "95% CI", "Isolates present", "Paper TP / FP / FN"], [[("**%s**" % e if e in KEY else e), f"{r.P[e]['clones_present']} / {r.P[e]['clones_typed']:,}", pc(r.P[e]["frac_present"]),
      f"{pc(r.P[e]['ci95_lo'])} to {pc(r.P[e]['ci95_hi'])}", f"{r.P[e]['samples_present']:,}", ("%d / %d / %d" % (r.truth(e)["TP"], r.truth(e)["FP"], r.truth(e)["FN"])) if r.truth(e) else ""]
      for e in sorted(r.E, key=lambda e: -r.P[e]["frac_present"])]) + "\n")
    w("Source: `tables/prevalence.tsv` (`out_default/prevalence.tsv`); agreement from `targets/truth_store.tsv`.\n")
    w("## 4. Agreement with the paper\n")
    w(fig_md("fig2_threshold_sweep", "Precision and recall against the paper's calls by present cut") + "\n")
    tuned_rows = []
    for e in ("Osiris-h3", "Janus-h4"):
        c = r.cut(e)
        if c:
            x0, x1 = r.sweep_at(e, 0.5), r.sweep_at(e, c)
            tuned_rows.append([e, c, x0["FP"], x1["FP"], x0["TP"], x1["TP"]])
    if tuned_rows:
        w(md_table(["Element", "Cut", "False pos. at 0.5", "False pos. at cut", "True pos. at 0.5", "True pos. at cut"], tuned_rows) + "\n")
    w("Elements with only a handful of paper positives (Janus-h3, Janus-h5) cannot be tuned. For bafA and Janus-h3/h5 the extra calls have hash fraction near 1.0 and normal depth, so they look genuinely present; the paper's calls may be incomplete there (for bafA the paper table is mostly a Nebuchadnezzar-h1 proxy). The cuts for Osiris-h3 and Janus-h4 were read off this same truth table, so their gain is optimistic until tested on strains not used to choose them.\n")
    w("Source: `tables/threshold_sweep.tsv` (`out_default/threshold_sweep.tsv`, truth `targets/truth_store.tsv`, 226 samples).\n")
    w("## 5. Spread over the population tree\n")
    w(fig_md("fig3_tree_heatmap", "Share of clone groups carrying each element along the tree") + "\n")
    w(f"The {sum(1 for h in D['heat'] if h['in_tree']):,} clone groups with a tip in the tree are laid out in tree order and pooled 6 per column. A solid band means the element sits in one part of the tree; speckle means it appears here and there. The strip shows the dominant lineage per column (lineage = group at ANI 0.998). {len(D['heat']) - sum(1 for h in D['heat'] if h['in_tree']):,} clone groups have no tip in the tree and are not shown; fine order inside a lineage is not reliable.\n")
    w("Source: `out_default/clone_matrix.tsv` and `tree_order.tsv`; tree `query/distances/tree_fastme.nwk`.\n")
    w("## 6. Gains and losses on the tree\n")
    w(fig_md("fig4_gain_loss", "Minimum changes on the tree against clone groups present") + "\n")
    w(md_table(["Element", "Present clones", "Min. changes", "Changes / present", "Gains", "Losses", "Prevalence"], [[("**%s**" % p["element"] if p["element"] in KEY else p["element"]), f"{p['clone_tips_present']:,}", f"{p['fitch_score']:,}",
      p["score_per_present"], p["gains_root_absent"], p["losses_root_absent"], pc(r.P[p["element"]]["frac_present"])] for p in sorted(D["default"]["parsimony"], key=lambda p: -p["score_per_present"])]) + "\n")
    w("Fitch parsimony on the distance tree, one tip per clone group; the root state is fixed by the data, so root-present gives the same split. Read the counts as a pattern: multifurcations are resolved left to right, and overcalled elements are inflated.\n")
    w("Source: `tables/parsimony.tsv`; method in the `starship_dynamics.py` docstring, tested against a Sankoff dynamic program on 200 random trees.\n")
    w("## 7. bafA and Nebuchadnezzar-h1\n")
    w(fig_md("fig5_bafA_groups", "bafA by country and by lineage") + "\n")
    w(md_table(["", "Neb-h1 present", "Neb-h1 absent"], [["bafA present", bn["11"], bn["10"]], ["bafA absent", bn["01"], bn["00"]]]) + "\n")
    w(f"`bafA` here is the AF293-type allele; divergent alleles are in section 8. {sum(bn.values()):,} clone groups are typed for both. Country and lineage are confounded (sampling is uneven; lineage 150 contributes many bafA-positive clone groups), so differences are descriptive, not tests.\n")
    w("Source: `tables/prevalence_by_country.tsv`, `tables/prevalence_by_lineage.tsv`, `out_default/clone_matrix.tsv`.\n")
    if fp_:
        bfam = D["baf_family"]; pop = bfam["population"]; vals = bfam["validation"]
        w("## 8. bafA beyond the AF293-type allele (baf family)\n")
        w("The Starship screen's bafA target is the AF293-type genomic copy and has 2 usable hashes, so a divergent allele can be invisible to it. "
          "Assembly counts for the baf genes (`/bigdata/stajichlab/shared/projects/A_fumigatus/baf/strain_gene_counts*.tab`; the `bafD` column is bafA, 96% of its 21-mers) show this, "
          "so the family was screened with one target per gene, built from the assembly hit sequences (12 bafA, 6 bafB, 36 bafC and 67 Pseudobaf haplotypes, joined with N gaps). "
          "Gene presence = hit count at or above a cut (2 specific hashes for bafA, bafC and Pseudobaf; bafB has a single specific hash, so 1). The gene family is described in Kowalski et al. 2021, "
          "[mBio 12(1)](https://doi.org/10.1128/mBio.03579-20).\n")
        w(fig_md("fig7_baf_family", "baf family prevalence by clone group, and strains with an assembly copy that the screen finds") + "\n")
        vr = []
        for tname, v in vals.items():
            for g in ("bafA", "bafB", "bafC", "Pseudobaf"):
                x = v["genes"][g][str(bfam["cuts"][g])]
                vr.append([tname, g, bfam["cuts"][g], f"{x['tp']} / {x['tp'] + x['fn']}", x["fp"], pc(x["tp"] / max(x["tp"] + x["fn"], 1), 0)])
            x = v["previous_single_target_bafA"]
            vr.append([tname, "bafA, single AF293-type target (before)", "", f"{x['tp']} / {x['tp'] + x['fn']}", x["fp"], pc(x["tp"] / max(x["tp"] + x["fn"], 1), 0)])
        w("Agreement with the assembly counts (strains in the store that are not flagged):\n")
        w(md_table(["Counts table", "Gene", "Cut (hits)", "Assembly-positive found", "Extra calls", "Recall"], vr, ["l", "l", "r", "r", "r", "r"]) + "\n")
        w("The extra bafA calls against the older table are strains the newer table lists as carrying bafA. The bafA misses are IFM_59359 (its read sets and assembly disagree, see the earlier mixup) "
          "and one more strain that is not explained yet. Relative depth separates 2-copy from 1-copy strains only weakly, so it is not a copy-number estimate.\n")
        pg = pop["genes"]
        w(f"Across the {pop['kept']:,} kept isolates:\n")
        w(md_table(["Gene", "Cut (hits)", "Isolates", "Share", "Clone groups", "Share"],
                   [["bafA, AF293-type (Starship screen)", "", f"{r.P['bafA']['samples_present']:,}", pc(r.P['bafA']['samples_present'] / pop['kept']), f"{r.P['bafA']['clones_present']} of {r.P['bafA']['clones_typed']:,}", pc(r.P['bafA']['frac_present'])]]
                   + [[("bafA, all alleles" if g == "bafA" else g), x["cut"], f"{x['isolates']:,}", pc(x["isolates"] / pop["kept"]), f"{x['clone_groups']} of {x['clone_groups_total']:,}", pc(x["clone_groups"] / x["clone_groups_total"])] for g, x in pg.items()],
                   ["l", "r", "r", "r", "r", "r"]) + "\n")
        ov, nv = pop["old_vs_new_bafA"], pop["bafA_all_vs_neb_h1"]
        w(f"Isolates: the AF293-type target and the all-alleles call agree on {ov['both']:,} present and {ov['neither']:,} absent; {ov['new_only']:,} are called only by the all-alleles screen and {ov['old_only']:,} only by the old target (one hash in common). "
          f"Against Nebuchadnezzar-h1: {nv['both']:,} isolates have both, {nv['bafA_only']:,} have bafA without Nebuchadnezzar-h1, {nv['neb_only']} have Nebuchadnezzar-h1 without bafA. "
          f"None of the {pop['newly_called']['isolates']} newly called isolates has Nebuchadnezzar-h1 present, so these look like bafA copies that are not on that Starship. "
          "That is a pattern in sketches: the copies have not been placed in genomes. The strong candidate list in section 10 is unchanged.\n")
        w("Source: `tables/baf_family_summary.json` (`baf_family.py summarize`), screen `query/screen_baf_family/` in the store (job `31996411`, `--max-solo 0`), counts from the shared folder above. "
          "The strain-level comparison names local strains and stays in the results folder.\n")
    w("## 9. Which elements travel together\n")
    w(fig_md("fig6_cooccurrence", "Pairwise phi of presence across clone groups") + "\n")
    pairs = sorted([x for x in D["cooccurrence"] if x["phi"] != "" and fam(x["element_1"]) != fam(x["element_2"])], key=lambda x: -abs(float(x["phi"])))[:10]
    w("Strongest pairs from different families:\n")
    w(md_table(["Element 1", "Element 2", "Clone groups typed", "Both", "Only 1", "Only 2", "Neither", "phi"], [[x["element_1"], x["element_2"], f"{x['clones_typed']:,}", x["both"], x["only_1"], x["only_2"], f"{x['neither']:,}", f"{float(x['phi']):.2f}"] for x in pairs]) + "\n")
    w("Elements from the same family look correlated because they share sequence. Phi has no significance test here; treat strong pairs as leads.\n")
    w("Source: `tables/cooccurrence.tsv`.\n")
    w("## 10. The strains\n")
    cnt = [[tt, len(tier(cands, tt)), len(tier(cands, tt, True)), len(tier(pub, tt)), len(tier(loc, tt))] for tt in ("strong", "probable", "discordant")]
    w(md_table(["Tier", "All", "After flagged removed", "Public", "Local"], cnt) + "\n")
    w("**strong** = both bafA hashes found and Nebuchadnezzar-h1 present. **probable** = bafA present and Nebuchadnezzar-h1 partial. **discordant** = one of the two present without the other (check read depth before counting these). Flagged = other species, k-mer mixture or curator exclusion.\n")
    w("Local (Population_Genomics) isolates count in these numbers but are not listed by name, as in the site report. The public isolates are in [`tables/candidates_public.tsv`](tables/candidates_public.tsv).\n")
    bys = collections.Counter((c["study"], c["country"]) for c in tier(pub, "strong", True))
    w("Strong public candidates by study (top 15):\n")
    w(md_table(["Study", "Country", "Isolates"], [[s, c, n] for (s, c), n in bys.most_common(15)], ["l", "l", "r"]) + "\n")
    w("Source: `tables/candidates_public.tsv` (from `bafA_neb_h1_candidates_v2.tsv`, public rows only).\n")
    if r.PT:
        w("## 11. Sensitivity to the cuts\n")
        rows = [[e, f"{r.P[e]['clones_present']} / {r.P[e]['clones_typed']:,}", pc(r.P[e]["frac_present"]), f"{r.PT[e]['clones_present']} / {r.PT[e]['clones_typed']:,}", pc(r.PT[e]["frac_present"]),
                 r.PS[e]["score_per_present"], r.PST[e]["score_per_present"]] for e in ("Osiris-h3", "Janus-h4") if r.cut(e)]
        w(md_table(["Element", "Clones present (default)", "Share", "Clones present (tuned)", "Share", "Changes / present (default)", "Changes / present (tuned)"], rows) + "\n")
        w("The cut moves the prevalence by a few points. The scatter ratio hardly changes, so the conclusion that Osiris-h3 is scattered and Janus-h4 is clade-structured does not depend on it. bafA, Nebuchadnezzar-h1 and every other element are unaffected.\n")
    w("## 12. Limits\n")
    w("- The 20 haplotypes are those of the NII *A. fumigatus* pangenome study (Table S21). Starbase and the Gluck-Thaler set are not screened.")
    w("- Presence means most of an element's specific 21-mers are in the read sketch at 30x depth. It does not give an insertion site, an intact element, or copy number.")
    w("- The paper table used as truth covers 226 isolates, and for bafA it is mostly a Nebuchadnezzar-h1 proxy. Agreement is a check, not proof.")
    w("- Cuts for Osiris-h3 and Janus-h4 were chosen on the same truth table. They need a hold-out check.")
    w(f"- Gain/loss counts come from a sketch distance tree with {nTree:,} of {nClones:,} clone groups. They describe a pattern.\n")
    w("## 13. Traceability and regeneration\n")
    w(md_table(["What", "File or job", "Made by"], [
      ["Sketches, QC, clone groups, tree", f"`{M.get('store', '')}/{{sigs,query}}`", "store build (PR #315), queries job 31975386"],
      ["Starship screen (21 targets)", "`query/screen_starships_bafA/`", f"rescreen job `{M.get('rescreen_job', '?')}`"],
      ["Prevalence, matrix, parsimony, co-occurrence, sweep", f"`{M.get('results', '')}/out_default/`", f"`starship_dynamics.py` at `{M.get('code_commit', '?')}`, job `{M.get('default_job', '?')}`"],
      ["Tuned cuts", f"`{M.get('results', '')}/out_tuned/`", f"`--set Osiris-h3=0.7 Janus-h4=0.7`, job `{M.get('tuned_job', '?')}`"],
      ["Report data", f"`{M.get('results', '')}/report_data.json`", f"`starship_report_data.py` at `{M.get('code_commit', '?')}`"],
      *([["baf family screen", "`query/screen_baf_family/`, `tables/baf_family_summary.json`", "`baf_family.py` (targets, summarize), job `31996411`"]] if D.get("baf_family") else []),
      ["This file, figures, tables", "`report.md`, `figures/`, `tables/`", "`starship_report_md.py`"]], ["l", "l", "l"]) + "\n")
    w("Regenerate on the HPCC, from the results folder: `sbatch run_dynamics.sh; sbatch run_tuned.sh; sbatch run_report.sh` (runs both the page and this markdown). "
      "Tests: `pixi run -e test python -m pytest tests -k sketch_store`.\n")
    open(f"{out}/report.md", "w").write("\n".join(L) + "\n")
    print(f"report.md, 6 figures (png + pdf), {len(os.listdir(t))} tables, {len(pub)} public candidates, {len(loc)} local candidates counted only")


if __name__ == "__main__":
    sys.exit(main())
