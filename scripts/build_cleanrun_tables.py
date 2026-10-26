#!/usr/bin/env python3
"""Per-campaign result tables for the clean full-BFD run, with the curator's post-run reporting columns.

Reads one campaign's per-genome `detection_report.yaml` files and the v0.6.0 baseline reports of the same
genomes, and writes `genomes.tsv`, `loci.tsv` and `by_<rank>.tsv` (default `by_class.tsv`). Nothing in a
report is changed: every new column is computed from the reports AFTER the run (curator rulings 2026-10-22
and 2026-10-24; issues #129, #137, #138). The columns written by the v0.6.0 `analyze_asco.py` come first and
keep their meaning, so a reader of the old tables still works.

Report layouts (found by walking the folder, so any depth works):
  Ascomycota     wave_N/small_NNN/runs/<genome>/detection_report.yaml  (clean run; v0.6.0: wave_N/runs/<genome>/)
  Basidiomycota  runs/{small_NNN,large_NNN}/runs/<genome>/...          (clean run; v0.6.0: small_N/runs/<genome>/)
  Mucoromycotina runs/<SOURCE>__<id>/...                               (--layout mucoro; SOURCE is BFD, LCG or JENA)
Either source may also be a `.tar.zst` / `.tar.gz` / `.tar` holding the same tree (for example the committed
`reports_all.tar.zst` of a v0.6.0 campaign); it is streamed, never unpacked to disk.
A genome is the name of the folder holding its report. If one genome has several folders, a non-empty report
wins over a missing one and, among non-empty ones, the lexicographically last path wins (a re-run chunk such
as `Pucciniales_long` sorts after the chunk it repeats); the count of such genomes is printed.

The listed genomes (rows of genomes.tsv) come from `lists/*.tsv` (first column; files sorted by name, rows in
file order), or for --layout mucoro from `inputs.tsv` (`source`, `id`). A listed genome with no report is
`no_report`, or `failed` when its folder holds a `stderr.log`. Reports of unlisted genomes are counted and left
out. With no list at all, every genome with a report is a row, sorted by name.

New per-locus columns (loci.tsv):
  margin_basis       `bits` when the locus has two or more idiomorph candidates (the bit-score vote, or the
                     profile-HMM classifier's scores), `identity_fallback` when `idiomorph_margin` is set with
                     fewer than two candidates (the overlap-resolution identity margin, pipeline.py, #137,
                     #145), blank when there is no margin.
  baseline_confidence  confidence of the v0.6.0 locus matched to this one (same contig, overlapping span; the
                     same detected family is preferred, then the larger overlap; each baseline locus is used
                     once, new loci taken in report order). Blank when nothing matches.
  tier_note          only when the confidence tier DROPPED against the matched v0.6.0 locus, else blank. The
                     detector's `confidence` is kept as reported. First rule that applies:
                       (blank)              the locus class also changed (the class change is the reason);
                       margin_fallback      margin_basis is identity_fallback and the margin is under 5 (the
                                            identity-point fallback crossed the tier threshold; #137);
                       vote_collapse        a bit margin under 5 (a fragment vote near a tie);
                       polished_model_lost  fewer polished gene models than in v0.6.0 (e.g. CTG-yeast MTLA1);
                       distant_fragment     a gene that is unpolished and not in `confidence_ignored_genes`
                                            now, and was not so in v0.6.0 (e.g. a distant MAT1-1-2 fragment);
                       (blank)              none of these (counted as `other` in the printed summary).
                     So a blank tier_note means "no listed mechanism", not "no drop": a drop is seen by comparing
                     `confidence` with `baseline_confidence` (curator 2026-10-24: these stay blank).
  label_unsupported  #138, `yes` or `no`: `yes` when label_unsupported_reason is not empty.
  label_unsupported_reason  `bit_margin<20` when margin_basis is `bits` and the margin is under 20, and/or
                     `genus:<G>` for every locus of Aureobasidium, Friedmanniomyces or Xylaria; joined by `;`.
                     The identity fallback is NOT marked (that question is #145). The idiomorph is unchanged.
                     Only loci that report an idiomorph are marked: an `undetermined` (or empty) idiomorph is
                     never marked, by either rule (curator 2026-10-24).
  high_basis         #129, for `high` Ascomycota:MAT loci with a called idiomorph in the scope families
                     (default Pleosporaceae; --high-basis-scope all for every family). Extra genes are the
                     genes other than the called idiomorph's core gene and the scoped flanks GAP1/BGL1 (best
                     record per gene: polished first, then coverage). An extra gene is `opposite` when it is
                     the other idiomorph's core gene or its `lost_to` is the called core gene.
                       supported   some non-opposite extra gene is polished, or reaches >=50% identity with
                                   >=30% coverage;
                       tie_break   every extra gene is opposite or is unpolished with <30% coverage (also
                                   when there is no extra gene);
                       grey        neither (e.g. an unpolished extra gene with >=30% coverage, <50% identity).
                     Blank outside the scope, below `high`, or for an `undetermined` (homothallic) locus.
                     (Curator 2026-10-24: the scope stays the BFD FAMILY label Pleosporaceae.)
  high_basis_outside_scope  the same rule for `high` Ascomycota:MAT idiomorph calls OUTSIDE the scope families,
                     reported separately so its counts never mix with the scoped column (blank inside the scope).
  defective_genes   genes of the locus whose polished model has a frameshift or an internal stop codon (issue #151,
                     `defective_genes` in the report), `|`-joined; blank when none or when the report predates the field.
                     Report-only: the tier and `polished_genes` do not read it. It says "defective in this assembly";
                     a real pseudogene shows in many isolates of a species, an assembly indel in few.
New per-genome columns (genomes.tsv): genus, baseline_status, change_vs_baseline (the phase-5 judgement's
category, computed on the first locus of each sorted locus list, so it can differ from the per-locus view),
n_tier_drop (loci whose tier dropped), n_tier_note, n_label_unsupported, high_basis (distinct values over the
genome's loci), n_defective_loci (loci with at least one defective gene). by_class.tsv ends with `defective_loci`.

Taxonomy: --samples (BFD samples.csv: ASMID, SPECIES, CLASS, ORDER, FAMILY, GENUS). For --layout mucoro, the
names and orders come from --baseline-calls (the v0.6.0 `calls.tsv`) and then `inputs.tsv`.

Deterministic: the same inputs give byte-identical outputs (no dates, sorted where order is not given).
Needs PyYAML (and `zstd` for .tar.zst). Usage, for example:
  python3 scripts/build_cleanrun_tables.py --campaign results/2026-10-23_ascomycota_cleanrun \\
      --baseline results/2026-10-03_ascomycota_v060/reports_all.tar.zst --samples samples.csv \\
      --out results/2026-10-23_ascomycota_cleanrun
"""
import argparse
import csv
import io
import os
import re
import statistics
import subprocess
import sys
import tarfile
from collections import Counter, defaultdict

import yaml

LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)

POLISHED = frozenset({"polished_agree", "polished_single", "polished_disagree"})
CONF_RANK = {"high": 3, "medium": 2, "low": 1}
TIER_MARGIN = 5.0          # DEFAULT_MIN_IDIOMORPH_MARGIN (detect/family_registry.py): the margin tier cap
LABEL_MARGIN_BITS = 20.0   # #138: a bit margin under 20 does not support the idiomorph label
LABEL_GENERA = ("Aureobasidium", "Friedmanniomyces", "Xylaria")   # #138 genus rule (curator 2026-10-22/24)
HIGH_BASIS_FAMILY_CALLED = "Ascomycota:MAT"
HIGH_BASIS_SCOPE = ("Pleosporaceae",)
CORE_GENE = {"MAT1-1": "MAT1-1-1", "MAT1-2": "MAT1-2-1"}
SCOPED_FLANKS = frozenset({"GAP1", "BGL1"})
SPAN_SHIFT_BP = 1000       # a coordinate shift up to this is not a change (judgement, as in the pilot)

TIER_NOTES = ("margin_fallback", "vote_collapse", "polished_model_lost", "distant_fragment")

GENOME_COLS = ["genome", "species", "class_", "order", "family", "status", "routing", "n_loci", "families_called", "wall_s",
               "genus", "not_searched_reason", "baseline_status", "change_vs_baseline", "n_tier_drop", "n_tier_note", "n_label_unsupported",
               "high_basis", "n_defective_loci"]
LOCUS_COLS = ["genome", "class_", "order", "family_called", "contig", "start", "end", "idiomorph", "confidence", "locus_class",
              "detection_pass", "genes_found", "verification",
              "genus", "family", "polished_genes", "idiomorph_margin", "margin_basis", "n_candidates", "baseline_confidence",
              "tier_note", "label_unsupported", "label_unsupported_reason", "high_basis", "high_basis_outside_scope",
              "defective_genes"]


# ----------------------------------------------------------------------------------------------- reading reports

def _fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def summarise_locus(x):
    """The fields of one `detected` entry that the tables and the rules need."""
    v = x.get("verification")
    note = x.get("homothallic_note") if isinstance(x.get("homothallic_note"), dict) else {}
    return {
        "family": x.get("family"), "contig": x.get("contig"), "start": x.get("start"), "end": x.get("end"),
        "idiomorph": x.get("idiomorph"), "confidence": x.get("confidence"), "locus_class": x.get("locus_class"),
        "detection_pass": x.get("detection_pass"), "genes_found": list(x.get("genes_found") or []),
        "verification": v.get("status") if isinstance(v, dict) else (v or ""),
        "polished_genes": x.get("polished_genes"), "idiomorph_margin": _fnum(x.get("idiomorph_margin")),
        "n_candidates": len(x.get("idiomorph_candidates") or []),
        "vote_idiomorph": note.get("vote_idiomorph"),
        "ignored": sorted(set(x.get("confidence_ignored_genes") or [])),
        "defective": sorted(set(x.get("defective_genes") or [])),   # genes whose polished model has a frameshift or an internal stop
        "genes": [{"gene": e.get("gene"), "status": e.get("status"), "identity": _fnum(e.get("identity")),
                   "coverage": _fnum(e.get("coverage")), "lost_to": e.get("lost_to")}
                  for e in (x.get("gene_evidence") or [])],
    }


def summarise_report(text):
    """Compact summary of one detection_report.yaml (text or bytes); None when it is empty or unreadable."""
    if not text:
        return None
    doc = yaml.load(text, Loader=LOADER) or {}
    r = doc.get("routing")
    return {
        "routing": str(r.get("mode")) if isinstance(r, dict) else str(doc.get("routing_mode", "")),
        "not_searched_reason": doc.get("not_searched_reason") or "",
        "loci": [summarise_locus(x) for x in (doc.get("detected") or [])],
    }


def _member_kind(name):
    base = os.path.basename(name)
    if base in ("detection_report.yaml", "wall_seconds", "stderr.log"):
        return os.path.dirname(name), base
    return None, None


def _open_tar(path):
    if path.endswith(".zst"):
        proc = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE)
        return proc, tarfile.open(fileobj=proc.stdout, mode="r|")
    return None, tarfile.open(path, mode="r|*")


def read_campaign(src):
    """({run folder: run}, {list file name: text}) for a campaign folder or a tarball of one.

    A run is {"path": report path or None, "summary": parsed summary (tarballs) or None, "has_report": bool,
    "wall": str, "stderr": bool}. Folder reports are parsed later (`parse_runs`); tarball reports are parsed
    while the tarball streams. Run folder paths are relative to `src` and normalised (no leading ./)."""
    runs = defaultdict(lambda: {"path": None, "summary": None, "has_report": False, "wall": "", "stderr": False})
    lists = {}
    if os.path.isdir(src):
        for root, dirs, files in os.walk(src):
            dirs.sort()
            rel_root = os.path.normpath(os.path.relpath(root, src))
            for f in sorted(files):
                if f.endswith(".tsv") and (rel_root == "lists" or (rel_root == "." and f == "inputs.tsv")):
                    with open(os.path.join(root, f)) as fh:
                        lists[f if rel_root == "." else f"lists/{f}"] = fh.read()
            if "detection_report.yaml" in files:
                p = os.path.join(root, "detection_report.yaml")
                runs[rel_root].update(path=p, has_report=os.path.getsize(p) > 0)
            if "wall_seconds" in files:
                with open(os.path.join(root, "wall_seconds")) as fh:
                    runs[rel_root]["wall"] = fh.read().strip()
            if "stderr.log" in files:
                runs[rel_root]["stderr"] = True
        return dict(runs), lists
    proc, tf = _open_tar(src)
    with tf:
        for m in tf:
            if not m.isfile():
                continue
            name = os.path.normpath(m.name)
            if re.fullmatch(r"lists/[^/]+\.tsv|inputs\.tsv", name):
                lists[name] = tf.extractfile(m).read().decode()
                continue
            rel, kind = _member_kind(name)
            if kind == "detection_report.yaml":
                data = tf.extractfile(m).read()
                runs[rel].update(summary=summarise_report(data), has_report=bool(data))
            elif kind == "wall_seconds":
                runs[rel]["wall"] = tf.extractfile(m).read().decode().strip()
            elif kind == "stderr.log":
                runs[rel]["stderr"] = True
    if proc is not None and proc.wait():
        raise RuntimeError(f"zstd failed on {src}")
    return dict(runs), lists


def choose_run(candidates):
    """The run folder to use for a genome with several: `candidates` is [(path, has_non_empty_report)]; a non-empty
    report wins, then the lexicographically last path. Shared with build_locus_shards.py, so the tables and the
    locus maps always describe the same run."""
    return max(candidates, key=lambda c: (bool(c[1]), c[0]))[0]


def pick_runs(runs):
    """{genome: run} choosing one folder per genome (see the module docstring); also returns the duplicate count."""
    by_genome = defaultdict(list)
    for rel in sorted(runs):
        by_genome[os.path.basename(rel)].append(rel)
    out, dup = {}, 0
    for g, rels in by_genome.items():
        if len(rels) > 1:
            dup += 1
        out[g] = runs[choose_run([(r, runs[r]["has_report"]) for r in rels])]
    return out, dup


def _summarise_path(path):
    if not path or not os.path.getsize(path):
        return None
    with open(path, "rb") as fh:
        return summarise_report(fh.read())


def parse_runs(runs, jobs=1):
    """Fill in the summaries of folder reports, in sorted genome order (the result does not depend on `jobs`)."""
    todo = sorted(g for g, r in runs.items() if r["path"] and r["summary"] is None)
    paths = [runs[g]["path"] for g in todo]
    if jobs > 1 and len(paths) > 1:
        from multiprocessing import Pool
        with Pool(jobs) as pool:
            sums = list(pool.imap(_summarise_path, paths, chunksize=64))
    else:
        sums = [_summarise_path(p) for p in paths]
    for g, s in zip(todo, sums):
        runs[g]["summary"] = s
    return runs


def listed_genomes(lists, layout):
    """Genome ids in list order. Mucoromycotina ids are `<SOURCE>__<id>`, as the run folders are named."""
    out = []
    if layout == "mucoro":
        text = lists.get("inputs.tsv")
        if text:
            for r in csv.DictReader(io.StringIO(text), delimiter="\t"):
                if r.get("source") and r.get("source") != "source":
                    out.append(f"{r['source']}__{r['id']}")
        return out
    for name in sorted(k for k in lists if k.startswith("lists/")):
        for line in lists[name].splitlines():
            if line.strip():
                out.append(line.split("\t")[0])
    return out


# ----------------------------------------------------------------------------------------------- taxonomy

def load_samples(path):
    if not path:
        return {}
    with open(path, newline="") as fh:
        return {r["ASMID"]: r for r in csv.DictReader(fh)}


def genus_from_species(species):
    tok = (species or "").replace("[", "").replace("]", "").split()
    return tok[0] if tok and re.fullmatch(r"[A-Z][a-z-]+", tok[0]) else ""


def taxonomy(genome, samples, mucoro_names, inputs_rows, mucoro=False):
    """(species, class_, order, family, genus) for one genome id (`SOURCE__id` when `mucoro`)."""
    if mucoro:
        name, order = mucoro_names.get(genome, ("", ""))
        inp = inputs_rows.get(genome, {})
        name = name or inp.get("species", "") or genome.partition("__")[2].replace("_", " ")
        order = order or inp.get("order", "")
        return name, "Mucoromycotina", order, "", genus_from_species(name)
    s = samples.get(genome, {})
    species = s.get("SPECIES", "")
    return species, s.get("CLASS", "") or "NA", s.get("ORDER", ""), s.get("FAMILY", ""), s.get("GENUS", "") or genus_from_species(species)


def load_mucoro_names(calls_path):
    """{SOURCE__genome: (name, order)} from a v0.6.0 Mucoromycotina calls.tsv."""
    out = {}
    if calls_path:
        with open(calls_path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("source") and r["source"] != "source":
                    out.setdefault(f"{r['source']}__{r['genome']}", (r.get("name", ""), r.get("bfd_order", "")))
    return out


# ----------------------------------------------------------------------------------------------- the rules

def margin_basis(locus):
    if locus["idiomorph_margin"] is None:
        return ""
    return "bits" if locus["n_candidates"] >= 2 else "identity_fallback"


def label_unsupported(locus, genus):
    """#138 reasons for one locus ('' when the label is not marked). A locus that reports no idiomorph
    (`undetermined`, e.g. a homothallic candidate) has no label to support and is never marked."""
    if locus["idiomorph"] in (None, "", "undetermined"):
        return ""
    why = []
    if margin_basis(locus) == "bits" and locus["idiomorph_margin"] < LABEL_MARGIN_BITS:
        why.append(f"bit_margin<{LABEL_MARGIN_BITS:g}")
    if genus in LABEL_GENERA:
        why.append(f"genus:{genus}")
    return ";".join(why)


def _unpolished_counted(locus):
    ign = set(locus["ignored"])
    return {e["gene"] for e in locus["genes"] if e["status"] not in POLISHED and e["gene"] not in ign}


def tier_note(new, old):
    """The reason a locus' tier dropped against its v0.6.0 match, '' when it did not drop or no rule applies.

    Returns (note, dropped)."""
    if old is None or CONF_RANK.get(new["confidence"], 0) >= CONF_RANK.get(old["confidence"], 0):
        return "", False
    if new["locus_class"] != old["locus_class"]:
        return "", True
    m = new["idiomorph_margin"]
    if m is not None and m < TIER_MARGIN:
        return ("margin_fallback" if new["n_candidates"] <= 1 else "vote_collapse"), True
    if (new["polished_genes"] or 0) < (old["polished_genes"] or 0):
        return "polished_model_lost", True
    if _unpolished_counted(new) - _unpolished_counted(old):
        return "distant_fragment", True
    return "", True


def high_basis(locus, tax_family, scope):
    """#129 split of a `high` Ascomycota:MAT idiomorph call; '' when it does not apply."""
    if locus["family"] != HIGH_BASIS_FAMILY_CALLED or locus["confidence"] != "high":
        return ""
    if scope is not None and tax_family not in scope:
        return ""
    core = CORE_GENE.get(locus["idiomorph"])
    if core is None:      # undetermined (homothallic candidate): the split does not apply
        return ""
    opposite = {g for k, g in CORE_GENE.items() if k != locus["idiomorph"]}
    best = {}
    for e in locus["genes"]:
        g = e["gene"]
        if g == core or g in SCOPED_FLANKS:
            continue
        rec = (e["status"] in POLISHED, e["identity"] or 0.0, e["coverage"] or 0.0,
               g in opposite or (e["lost_to"] is not None and e["lost_to"] == core))
        if g not in best or (rec[0], rec[2]) > (best[g][0], best[g][2]):
            best[g] = rec
    if any(not opp and (pol or (idt >= 50 and cov >= 30)) for pol, idt, cov, opp in best.values()):
        return "supported"
    if all(opp or (not pol and cov < 30) for pol, idt, cov, opp in best.values()):
        return "tie_break"
    return "grey"


def match_baseline(new_loci, old_loci):
    """For each new locus, the index of its v0.6.0 locus or None (see `baseline_confidence`)."""
    used, out = set(), []
    for x in new_loci:
        best, key = None, None
        if x["start"] is None or x["end"] is None:   # no span: nothing to overlap (a null span never aborts the table)
            out.append(None)
            continue
        for j, o in enumerate(old_loci):
            if j in used or o["contig"] != x["contig"] or o["start"] is None or o["end"] is None:
                continue
            ov = min(int(x["end"]), int(o["end"])) - max(int(x["start"]), int(o["start"])) + 1
            if ov <= 0:
                continue
            k = (o["family"] == x["family"], ov, -j)
            if key is None or k > key:
                best, key = j, k
        if best is not None:
            used.add(best)
        out.append(best)
    return out


def _tup(x):
    return (x["family"] or "", x["idiomorph"] or "", x["confidence"] or "", x["locus_class"] or "", x["contig"] or "",
            int(x["start"] or 0), int(x["end"] or 0), "|".join(x["genes_found"]))


def change_category(new, old):
    """Genome-level category of the phase-5 judgement (cmp.py, then minor.py's split).

    A genome with no baseline report counts as uncalled in the baseline, as in the judgement."""
    if new is None:
        return "no new report"
    nl, ol = sorted(_tup(x) for x in new["loci"]), sorted(_tup(x) for x in (old or {}).get("loci", []))
    if not nl and not ol:
        return "uncalled both"
    if nl and not ol:
        return "NEW call"
    if ol and not nl:
        return "LOST call"
    if len(nl) != len(ol):
        return "locus count changed"
    o1, n1 = ol[0], nl[0]
    if o1[1] != n1[1]:
        return "idiomorph -> undetermined" if n1[1] == "undetermined" else "idiomorph flip"
    if CONF_RANK.get(n1[2], 0) > CONF_RANK.get(o1[2], 0):
        return "confidence up"
    if CONF_RANK.get(n1[2], 0) < CONF_RANK.get(o1[2], 0):
        return "confidence DOWN"
    if o1[3] != n1[3]:
        return "class changed"
    if abs(o1[5] - n1[5]) > SPAN_SHIFT_BP or abs(o1[6] - n1[6]) > SPAN_SHIFT_BP:
        return "span changed >1 kb"
    if sorted(x[7] for x in ol) != sorted(x[7] for x in nl):
        return "genes_found changed only"
    for p, q in zip(ol[1:], nl[1:]):
        if any(p[i] != q[i] for i in (0, 1, 2, 3, 4, 7)):
            return "secondary locus changed"
    return "identical call"


# ----------------------------------------------------------------------------------------------- tables

def build(new, base, listed, samples, mucoro_names, inputs_rows, scope, mucoro=False):
    """Returns (genome rows, locus rows, per-category tallies) as lists of dicts. An empty `base` (no
    --baseline) leaves the baseline columns blank."""
    base_given = bool(base)
    G, L = [], []
    tally = Counter()
    for g in listed:
        run = new.get(g)
        species, cls, order, family, genus = taxonomy(g, samples, mucoro_names, inputs_rows, mucoro)
        row = dict(genome=g, species=species, class_=cls, order=order, family=family, status="no_report", routing="",
                   n_loci=0, families_called="", wall_s="", genus=genus, not_searched_reason="")
        b = base.get(g)
        bsum = b["summary"] if b else None
        row["baseline_status"] = (("called" if bsum["loci"] else "uncalled") if bsum else "no_report") if base_given else ""
        summ = run["summary"] if run else None
        if run:
            row["wall_s"] = run["wall"]
        if summ is None:
            if run and run["stderr"]:
                row["status"] = "failed"
            row.update(change_vs_baseline=change_category(None, bsum) if base_given else "", n_tier_drop=0,
                       n_tier_note=0, n_label_unsupported=0, high_basis="", n_defective_loci=0)
            G.append(row)
            continue
        loci = summ["loci"]
        row.update(status="called" if loci else "uncalled", n_loci=len(loci), routing=summ["routing"],
                   not_searched_reason=summ["not_searched_reason"],
                   families_called="|".join(sorted({x["family"] for x in loci})))
        row["change_vs_baseline"] = change_category(summ, bsum) if base_given else ""
        old_loci = bsum["loci"] if bsum else []
        match = match_baseline(loci, old_loci)
        n_note = n_lab = n_drop = 0
        bases = set()
        for x, j in zip(loci, match):
            o = old_loci[j] if j is not None else None
            note, dropped = tier_note(x, o)
            lab = label_unsupported(x, genus)
            hb = high_basis(x, family, scope)
            hb_out = high_basis(x, family, None) if scope is not None and family not in scope else ""
            n_note += bool(note)
            n_lab += bool(lab)
            n_drop += dropped
            if hb:
                bases.add(hb)
            if dropped:
                tally["tier_drop:" + (note or "other")] += 1
            L.append(dict(genome=g, class_=cls, order=order, family_called=x["family"], contig=x["contig"], start=x["start"],
                          end=x["end"], idiomorph=x["idiomorph"], confidence=x["confidence"], locus_class=x["locus_class"],
                          detection_pass=x["detection_pass"], genes_found="|".join(x["genes_found"]),
                          verification=x["verification"], genus=genus, family=family,
                          polished_genes="" if x["polished_genes"] is None else x["polished_genes"],
                          idiomorph_margin="" if x["idiomorph_margin"] is None else x["idiomorph_margin"],
                          margin_basis=margin_basis(x), n_candidates=x["n_candidates"],
                          baseline_confidence=o["confidence"] if o else "", tier_note=note,
                          label_unsupported="yes" if lab else "no", label_unsupported_reason=lab, high_basis=hb,
                          high_basis_outside_scope=hb_out, defective_genes="|".join(x["defective"])))
        row.update(n_tier_drop=n_drop, n_tier_note=n_note, n_label_unsupported=n_lab, high_basis="|".join(sorted(bases)),
                   n_defective_loci=sum(1 for x in loci if x["defective"]))
        G.append(row)
    return G, L, tally


def by_rank(G, L, rank):
    """Per-class (or per-order) table: the v0.6.0 columns, then locus counts for the new columns."""
    col = "class_" if rank == "class" else rank
    by = defaultdict(list)
    for r in G:
        by[r[col]].append(r)
    loc = defaultdict(list)
    for x in L:
        loc[x[col]].append(x)

    def q(v, p):
        v = sorted(v)
        return v[min(len(v) - 1, int(round((len(v) - 1) * p)))]

    head = [rank, "genomes", "reports", "called", "pct_called", "routing", "wall_median_s", "wall_p90_s", "wall_max_s",
            "loci", "n_label_unsupported", "tier_dropped", "tier_dropped_noted"] + list(TIER_NOTES) + ["high_supported", "high_tie_break", "high_grey",
                                    "outside_scope_high_supported", "outside_scope_high_tie_break", "outside_scope_high_grey",
                                    "defective_loci"]
    rows = []
    for c, rs in sorted(by.items(), key=lambda x: (-len(x[1]), x[0])):
        rep = [r for r in rs if r["status"] in ("called", "uncalled")]
        wall = [int(r["wall_s"]) for r in rep if str(r["wall_s"]).isdigit()]
        called = sum(r["status"] == "called" for r in rs)
        ls = loc.get(c, [])
        dropped = sum(1 for x in ls if x["baseline_confidence"] and
                      CONF_RANK.get(x["confidence"], 0) < CONF_RANK.get(x["baseline_confidence"], 0))
        notes = Counter(x["tier_note"] for x in ls)
        hb = Counter(x["high_basis"] for x in ls)
        ho = Counter(x["high_basis_outside_scope"] for x in ls)
        rows.append([c, len(rs), len(rep), called, f"{100 * called / max(1, len(rep)):.1f}",
                     "; ".join(f"{k}:{n}" for k, n in sorted(Counter(r["routing"] for r in rep).items(), key=lambda kv: (-kv[1], kv[0]))),
                     int(statistics.median(wall)) if wall else "", q(wall, .9) if wall else "", max(wall) if wall else "",
                     len(ls), sum(x["label_unsupported"] == "yes" for x in ls), dropped, sum(notes[n] for n in TIER_NOTES)]
                    + [notes[n] for n in TIER_NOTES] + [hb["supported"], hb["tie_break"], hb["grey"],
                                                                     ho["supported"], ho["tie_break"], ho["grey"],
                                                                     sum(1 for x in ls if x["defective_genes"])])
    return head, rows


def write_tsv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "") for c in cols] if isinstance(r, dict) else r)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--campaign", required=True, help="campaign folder, or a tarball of its report tree")
    ap.add_argument("--baseline", help="v0.6.0 report folder or reports_all.tar.zst (blank baseline columns without it)")
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--layout", choices=("bfd", "mucoro"), default="bfd",
                    help="bfd: Ascomycota/Basidiomycota (lists/*.tsv); mucoro: runs/<SOURCE>__<id> and inputs.tsv")
    ap.add_argument("--lists", nargs="*", help="genome list files to use instead of the campaign's lists/ (or inputs.tsv)")
    ap.add_argument("--samples", help="BFD samples.csv for taxonomy")
    ap.add_argument("--baseline-calls", help="mucoro: v0.6.0 calls.tsv, for genome names and orders")
    ap.add_argument("--rank", choices=("class", "order"), default="class", help="grouping of by_<rank>.tsv")
    ap.add_argument("--high-basis-scope", default=",".join(HIGH_BASIS_SCOPE),
                    help="comma-separated taxonomic families for high_basis, or 'all' (default %(default)s)")
    ap.add_argument("--jobs", type=int, default=1, help="processes for YAML parsing (output does not depend on it)")
    a = ap.parse_args(argv)

    runs, lists = read_campaign(a.campaign)
    if a.lists:
        lists = {}
        for p in a.lists:
            with open(p) as fh:
                lists["inputs.tsv" if a.layout == "mucoro" else f"lists/{os.path.basename(p)}"] = fh.read()
    picked, dup = pick_runs(runs)
    new = parse_runs(picked, a.jobs)
    base = {}
    if a.baseline:
        bruns, _ = read_campaign(a.baseline)
        bpicked, _ = pick_runs(bruns)
        base = parse_runs(bpicked, a.jobs)
    listed = listed_genomes(lists, a.layout)
    seen = set()
    listed = [g for g in listed if not (g in seen or seen.add(g))]
    if not listed:
        listed = sorted(g for g, r in new.items() if r["summary"] is not None)
    unlisted = sorted(set(g for g, r in new.items() if r["summary"] is not None) - set(listed))
    inputs_rows = {}
    if a.layout == "mucoro" and "inputs.tsv" in lists:
        for r in csv.DictReader(io.StringIO(lists["inputs.tsv"]), delimiter="\t"):
            if r.get("source") and r["source"] != "source":
                inputs_rows[f"{r['source']}__{r['id']}"] = r
    scope = None if a.high_basis_scope == "all" else tuple(s for s in a.high_basis_scope.split(",") if s)
    G, L, tally = build(new, base, listed, load_samples(a.samples), load_mucoro_names(a.baseline_calls), inputs_rows, scope,
                     a.layout == "mucoro")

    os.makedirs(a.out, exist_ok=True)
    genome_cols, locus_cols = list(GENOME_COLS), list(LOCUS_COLS)
    if a.layout == "mucoro":   # the run folder is SOURCE__id; split it back out for the site
        for r in G:
            r["source"], _, r["genome"] = r["genome"].partition("__")
        for x in L:
            x["source"], _, x["genome"] = x["genome"].partition("__")
        genome_cols.insert(0, "source")
        locus_cols.insert(0, "source")
    write_tsv(os.path.join(a.out, "genomes.tsv"), genome_cols, G)
    write_tsv(os.path.join(a.out, "loci.tsv"), locus_cols, L)
    head, rows = by_rank(G, L, a.rank)
    write_tsv(os.path.join(a.out, f"by_{a.rank}.tsv"), head, rows)

    print(f"{len(G)} genomes listed: {dict(sorted(Counter(r['status'] for r in G).items()))}; {len(L)} loci", file=sys.stderr)
    if dup:
        print(f"{dup} genomes had more than one run folder (rule: non-empty report, then the last path)", file=sys.stderr)
    if unlisted:
        print(f"{len(unlisted)} reports of unlisted genomes left out (e.g. {unlisted[:3]})", file=sys.stderr)
    print("tier drops by note:", dict(sorted(tally.items())), file=sys.stderr)
    print("label_unsupported loci:", sum(x["label_unsupported"] == "yes" for x in L), file=sys.stderr)
    print("high_basis:", dict(sorted(Counter(x["high_basis"] for x in L if x["high_basis"]).items())), file=sys.stderr)
    print("high_basis outside the scope:", dict(sorted(Counter(x["high_basis_outside_scope"] for x in L
                                                              if x["high_basis_outside_scope"]).items())), file=sys.stderr)


if __name__ == "__main__":
    main()
