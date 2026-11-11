#!/usr/bin/env python3
"""Starship (and bafA) presence/absence across a sketch store's clone groups: prevalence, co-occurrence, and gain/loss on the tree.

    starship_dynamics.py --store STORE [--query STORE/query] [--screen starships_bafA] [--tree FILE] [--tree-names FILE]
                         [--clone 0.9995] [--lineage 0.998] [--present 0.5] [--absent 0.2] [--set TARGET=FRACTION ...]
                         [--truth truth_store.tsv] [--min-group 10] --out OUTDIR

Inputs (written by the store build): QUERY/<screen>/screen.tsv (per sample and target: fraction of the target's specific hashes found),
QUERY/distances/samples_qc.tsv (mean ANI, clone groups), QUERY/mat/mat_calls.tsv (study, source, country), QUERY/qc/qc.tsv,
STORE/exclude.tsv, and the distance tree QUERY/distances/tree_fastme.nwk with tree_names.tsv (short tip id -> sample id; without
that file tip names are taken as sample ids).

State per sample and target: present if fraction >= --present (default 0.5, the screen's own cut), absent if < --absent (0.2),
else missing; partial hits are not read as loss. --set TARGET=FRACTION replaces the present cut for one target (the overlapping Janus
family and Osiris-h3 overcall at 0.5; threshold_sweep.tsv, written with --truth, shows precision and recall at each cut). Samples are
left out as in mat_frequency.py (curator exclusion, mean ANI < --min-mean-ani, k-mer QC mixture). State per clone group = majority
of its non-missing members (tie: missing).

Gain/loss: Fitch parsimony on the tree, one tip per clone group (its first member in the tree) carrying the clone state, every other
tip a wildcard. Multifurcations are resolved left to right, so scores are one resolution's. The score is the minimum number of changes.
Gains (absent -> present) and losses are counted for the two reconstructions that differ only in the root state (root absent, root
present), which gives a range. score / present-clones near 1 means scattered presence (repeated gain or transfer, or noisy calls);
near 0 means one clade. The tree is a sketch distance tree: read the counts as a pattern, not a rate.

Writes OUTDIR/clone_matrix.tsv, prevalence.tsv, prevalence_by_{lineage,country}.tsv, parsimony.tsv, cooccurrence.tsv,
threshold_sweep.tsv (with --truth) and summary.txt."""
import argparse, collections, csv, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mat_frequency import rd, wilson  # noqa: E402

ABSENT, PRESENT, EITHER = 1, 2, 3  # state sets as bitmasks for Fitch


def parse_newick(text):
    """Iterative Newick parser (the tree is ~5,000 tips deep). Returns (children, name): lists indexed by node, root = 0.
    A child always has a larger index than its parent."""
    children, name, stack, cur = [[]], [""], [], 0
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch == "(":
            children.append([]); name.append(""); children[cur].append(len(children) - 1)
            stack.append(cur); cur = len(children) - 1
        elif ch == ",":
            parent = stack[-1]
            children.append([]); name.append(""); children[parent].append(len(children) - 1); cur = len(children) - 1
        elif ch == ")":
            cur = stack.pop()
        elif ch == ";":
            break
        elif ch == ":":
            while i + 1 < n and text[i + 1] not in ",();":
                i += 1
        elif not ch.isspace():
            j = i
            while j < n and text[j] not in ",():;":
                j += 1
            name[cur] = text[i:j].strip("'"); i = j - 1
        i += 1
    return children, name


def fitch(children, name, states):
    """Fitch parsimony. states: tip name -> ABSENT | PRESENT | EITHER (tips not listed are wildcards).
    Returns (score, {'root_absent': (gains, losses), 'root_present': (gains, losses)}); gains + losses == score in both."""
    n = len(children)
    sets, chain, score = [EITHER] * n, [None] * n, 0
    for v in range(n - 1, -1, -1):
        kids = children[v]
        if not kids:
            sets[v] = states.get(name[v], EITHER)
            continue
        acc = sets[kids[0]]; seq = [acc]
        for c in kids[1:]:
            inter = acc & sets[c]
            if inter:
                acc = inter
            else:
                acc |= sets[c]; score += 1
            seq.append(acc)
        sets[v], chain[v] = acc, seq
    out = {}
    for label, pref in (("root_absent", ABSENT), ("root_present", PRESENT)):
        st, gains, losses = [0] * n, 0, 0
        st[0] = sets[0] if sets[0] in (ABSENT, PRESENT) else pref
        for v in range(n):
            kids = children[v]
            if not kids:
                continue
            x, seq = st[v], chain[v]
            for j in range(len(kids) - 1, 0, -1):
                child = x if x & sets[kids[j]] else sets[kids[j]]
                prev = x if x & seq[j - 1] else seq[j - 1]
                for b in (child, prev):
                    if b != x:
                        gains, losses = (gains + 1, losses) if x == ABSENT else (gains, losses + 1)
                st[kids[j]] = child; x = prev
            st[kids[0]] = x
        out[label] = (gains, losses)
    return score, out


def majority(vals):
    """vals: iterable of True (present) / False (absent) / None (missing) -> True / False / None (tie or nothing)."""
    vals = list(vals)
    p = sum(1 for v in vals if v is True); a = sum(1 for v in vals if v is False)
    return None if p == a else p > a


def sample_state(f, present, absent):
    return True if f >= present else False if f < min(absent, present) else None


def write_tsv(path, rows, cols=None):
    cols = cols or (list(rows[0]) if rows else ["empty"])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--store", required=True); ap.add_argument("--query"); ap.add_argument("--out", required=True)
    ap.add_argument("--screen", default="starships_bafA", help="name of the screen run (QUERY/screen_<name>/screen.tsv)")
    ap.add_argument("--tree"); ap.add_argument("--tree-names")
    ap.add_argument("--clone", default="0.9995"); ap.add_argument("--lineage", default="0.998")
    ap.add_argument("--present", type=float, default=0.5); ap.add_argument("--absent", type=float, default=0.2)
    ap.add_argument("--set", nargs="*", default=[], metavar="TARGET=FRACTION")
    ap.add_argument("--min-mean-ani", type=float, default=0.95)
    ap.add_argument("--min-group", type=int, default=10, help="groups with fewer clone groups are left out of the by-group tables")
    ap.add_argument("--truth", help="TSV with sample_id and one 0/1 column per target name (enables threshold_sweep.tsv)")
    a = ap.parse_args(argv)
    q = a.query or os.path.join(a.store, "query")
    overrides = {k: float(v) for k, v in (x.split("=") for x in a.set)}
    cut = collections.defaultdict(lambda: a.present, overrides)
    os.makedirs(a.out, exist_ok=True)
    mat, sq, qc = rd(f"{q}/mat/mat_calls.tsv"), rd(f"{q}/distances/samples_qc.tsv"), rd(f"{q}/qc/qc.tsv")
    excl, man = rd(os.path.join(a.store, "exclude.tsv")), rd(os.path.join(a.store, "manifest.tsv"))
    frac, targets = collections.defaultdict(dict), []
    for r in csv.DictReader(open(f"{q}/screen_{a.screen}/screen.tsv"), delimiter="\t"):
        frac[r["sample_id"]][r["target"]] = float(r["fraction"])
        if r["target"] not in targets:
            targets.append(r["target"])
    bad = [t for t in overrides if t not in targets]
    if bad:
        sys.exit(f"--set names not in the screen: {bad}")
    info, kept = {}, []
    for s in frac:
        d, f, m = sq.get(s, {}), qc.get(s, {}), mat.get(s, {})
        reason = ("curator_exclude" if s in excl else "other_species" if d.get("mean_ani") and float(d["mean_ani"]) < a.min_mean_ani
                  else "qc_mixture" if "mixture?" in (f.get("flag") or "") else "")
        country = (m.get("country") or "").split(":")[0].strip() or ("local" if m.get("source") == "local" else "unknown")
        info[s] = dict(source=m.get("source", ""), study=m.get("study_accession", ""), country=country, strain=man.get(s, {}).get("strain", ""),
                       excluded=reason, clone=d.get(f"clone_{a.clone}", "") or "solo:" + s, lineage=d.get(f"clone_{a.lineage}", "") or "none")
        if not reason:
            kept.append(s)
    # sample states and clone states (kept samples only)
    members = collections.defaultdict(list)
    for s in sorted(kept):
        members[info[s]["clone"]].append(s)
    cstate = {c: {t: majority(sample_state(frac[s].get(t, 0.0), cut[t], a.absent) for s in ss) for t in targets} for c, ss in members.items()}
    first = {c: info[ss[0]] for c, ss in members.items()}
    clones = sorted(members)
    write_tsv(f"{a.out}/clone_matrix.tsv", [dict(clone=c, n_members=len(members[c]), representative=members[c][0], strain=first[c]["strain"],
              study=first[c]["study"], country=first[c]["country"], lineage=first[c]["lineage"],
              **{t: {True: "1", False: "0", None: "NA"}[cstate[c][t]] for t in targets}) for c in clones])
    summary = [f"samples in the screen: {len(frac)}; kept {len(kept)}; excluded {dict(collections.Counter(i['excluded'] for i in info.values() if i['excluded']))}",
               f"clone groups at ANI {a.clone}: {len(clones)}; present cut {a.present}, absent cut {a.absent}; "
               f"overrides: {', '.join(f'{k}={v}' for k, v in overrides.items()) or 'none'}"]
    # prevalence overall and by group
    def prev_row(t, cs, **extra):
        typed = [c for c in cs if cstate[c][t] is not None]; k = sum(1 for c in typed if cstate[c][t]); lo, hi = wilson(k, len(typed))
        return dict(extra, element=t, clones=len(cs), clones_typed=len(typed), clones_present=k,
                    frac_present=f"{k / len(typed):.3f}" if typed else "", ci95_lo=f"{lo:.3f}" if typed else "", ci95_hi=f"{hi:.3f}" if typed else "")
    rows = []
    for t in targets:
        cs = collections.Counter(sample_state(frac[s].get(t, 0.0), cut[t], a.absent) for s in kept)
        rows.append(dict(prev_row(t, clones), samples_present=cs[True], samples_missing=cs[None], samples_absent=cs[False]))
    write_tsv(f"{a.out}/prevalence.tsv", rows)
    for gname, key in (("lineage", "lineage"), ("country", "country")):
        groups = collections.defaultdict(list)
        for c in clones:
            groups[first[c][key]].append(c)
        out = [prev_row(t, cs, group=g) for g, cs in sorted(groups.items(), key=lambda x: -len(x[1])) if len(cs) >= a.min_group for t in targets]
        write_tsv(f"{a.out}/prevalence_by_{gname}.tsv", out, ["group", "element", "clones", "clones_typed", "clones_present", "frac_present", "ci95_lo", "ci95_hi"])
    # co-occurrence over clones typed for both
    co = []
    for i, t1 in enumerate(targets):
        for t2 in targets[i + 1:]:
            both = [c for c in clones if cstate[c][t1] is not None and cstate[c][t2] is not None]
            n11 = sum(cstate[c][t1] and cstate[c][t2] for c in both); n10 = sum(cstate[c][t1] and not cstate[c][t2] for c in both)
            n01 = sum(not cstate[c][t1] and cstate[c][t2] for c in both); n00 = len(both) - n11 - n10 - n01
            den = math.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
            co.append(dict(element_1=t1, element_2=t2, clones_typed=len(both), both=n11, only_1=n10, only_2=n01, neither=n00,
                           phi=f"{(n11 * n00 - n10 * n01) / den:.3f}" if den else "", jaccard=f"{n11 / (n11 + n10 + n01):.3f}" if n11 + n10 + n01 else ""))
    write_tsv(f"{a.out}/cooccurrence.tsv", co)
    # parsimony on the tree
    tree_path = a.tree or f"{q}/distances/tree_fastme.nwk"
    if os.path.exists(tree_path):
        children, name = parse_newick(open(tree_path).read())
        names_path = a.tree_names or f"{q}/distances/tree_names.tsv"
        tip2sample = {}
        if os.path.exists(names_path):
            tip2sample = {r[0]: r[1] for r in csv.reader(open(names_path), delimiter="\t") if len(r) > 1}
        sample2tip = {v: k for k, v in tip2sample.items()} if tip2sample else {n: n for n in name if n}
        tips = {}
        for c in clones:
            rep = next((s for s in members[c] if s in sample2tip), None)
            if rep:
                tips[c] = sample2tip[rep]
        summary.append(f"tree tips: {sum(1 for v in range(len(children)) if not children[v])}; clone groups with a tip in the tree: {len(tips)} of {len(clones)}")
        par = []
        for t in targets:
            states = {tips[c]: (PRESENT if cstate[c][t] else ABSENT) for c in tips if cstate[c][t] is not None}
            score, rec = fitch(children, name, states)
            npres = sum(1 for v in states.values() if v == PRESENT)
            par.append(dict(element=t, clone_tips_typed=len(states), clone_tips_present=npres, fitch_score=score,
                            score_per_present=f"{score / npres:.2f}" if npres else "",
                            gains_root_absent=rec["root_absent"][0], losses_root_absent=rec["root_absent"][1],
                            gains_root_present=rec["root_present"][0], losses_root_present=rec["root_present"][1]))
        write_tsv(f"{a.out}/parsimony.tsv", par)
        summary.append("parsimony (clone tips; score = minimum changes; ratio near 1 = scattered, near 0 = one clade):")
        summary += [f"  {p['element']:20s} present {p['clone_tips_present']:5d}/{p['clone_tips_typed']:<5d} score {p['fitch_score']:4d} ratio {p['score_per_present']:>5s}  "
                    f"gains/losses root-absent {p['gains_root_absent']}/{p['losses_root_absent']}, root-present {p['gains_root_present']}/{p['losses_root_present']}" for p in par]
    else:
        summary.append(f"no tree at {tree_path}: parsimony skipped")
    # threshold sweep against the truth table
    if a.truth:
        truth = rd(a.truth); sweep = []
        for t in targets:
            if not truth or t not in next(iter(truth.values())):
                continue
            for th in [x / 10 for x in range(1, 11)]:
                tp = fn = fp = tn = 0
                for s, r in truth.items():
                    if s not in frac or r.get(t, "") not in ("0", "1"):
                        continue
                    pred, real = frac[s].get(t, 0.0) >= th, r[t] == "1"
                    tp += pred and real; fn += (not pred) and real; fp += pred and not real; tn += (not pred) and not real
                sweep.append(dict(element=t, present_cut=th, TP=tp, FN=fn, FP=fp, TN=tn,
                                  precision=f"{tp / (tp + fp):.3f}" if tp + fp else "", recall=f"{tp / (tp + fn):.3f}" if tp + fn else ""))
        write_tsv(f"{a.out}/threshold_sweep.tsv", sweep)
        summary.append(f"threshold sweep against {a.truth}: {len(sweep)} rows")
    open(f"{a.out}/summary.txt", "w").write("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    sys.exit(main())
