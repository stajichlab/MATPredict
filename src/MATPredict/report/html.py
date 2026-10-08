"""Render one `detect` run (`detection_report.yaml`) as a self-contained HTML page.

No network, no external CSS, fonts or scripts: the page opens from disk, an
e-mail attachment or a Galaxy history, and prints to PDF with the print
stylesheet below (Chrome and WeasyPrint; only CSS both support). Every string
from the report is escaped: contig names and record ids come from user files.
Design review of the first version: analysis/2026-10-07_report-design-review.md.
"""
from __future__ import annotations

from datetime import datetime, timezone
from html import escape

from MATPredict import __version__
from MATPredict.report.svg import figure_keys, locus_svg

_CONF_RANK = {"high": 3, "medium": 2, "low": 1}

LOCUS_CLASS = {
    "mat_locus": "Complete locus: MAT gene(s) with their flanking genes.",
    "idiomorph_gene_only": "A MAT gene with no flanking gene beside it; not a full locus.",
    "homothallic_candidate": "Genes of both idiomorphs at one locus.",
    "partial_locus": "Provisional: the call only tied the admission bar, or came from the relaxed pass.",
}

ROUTING = {
    "lineage": ("By taxonomy", "The genome's taxonomy falls inside a curated family's scope; only those families "
                               "were searched."),
    "phylum_fallback": ("Whole phylum", "No curated family covers this lineage, so every family of the genome's "
                                        "phylum was searched. Family labels are less reliable here."),
    "explicit_phylum": ("Phylum given", "The phylum was given when the run was started; taxonomy was not used."),
    "exhaustive": ("All families", "No curated family covers this genome, so every family was searched on request. "
                                   "Calls from another phylum's family are unconfirmed."),
    "not_searched": ("Not searched", "No curated family covers this genome."),
}

STATUS = {
    "polished_agree": "gene model",
    "polished_disagree": "models differ",
    "polished_single": "gene model",
    "unpolished": "hit only",
    "not_polish_candidate": "hit only",
}

METHOD = {
    "exonerate_refine": "exonerate",
    "miniprot_refine": "miniprot",
    "tblastn_genome": "tblastn",
    "diamond_proteome": "DIAMOND, proteome",
    "caax_scan": "motif scan (CAAX)",
}

BASIS = {"first_pass_identity": "higher identity"}
CLF_INPUT = {"polished_model": "gene model", "hsp_fragment": "hit fragment"}

WITHHELD = {
    "modelled_gene_bar": "fewer than 2 genes could be modelled",
    "below_fraction_floor": "too few of the family's genes found",
    "mat_gene_gate": "no MAT gene passed the MAT-gene check",
    "paralog_class": "core protein closer to a known non-MAT paralog",
    "flank_carried_core_outside_flank_span": "MAT hit outside the flanking genes' span",
}

ROLE = {"core_MAT": "core MAT", "flanking_conserved": "conserved flank", "flanking_variable": "variable flank"}

ARRAY_REASONS = {
    "array_size>=2": "array of 2 or more receptor loci",
    "single_locus": "a single receptor locus",
    "no_precursor_homology": "no pheromone-precursor homology hit",
    "precursor_homology": "a pheromone-precursor homology hit",
}


def array_reason(r) -> str:
    r = str(r)
    if r in ARRAY_REASONS:
        return ARRAY_REASONS[r]
    if r.startswith("strict_caax_orfs="):
        n = r.split("=", 1)[1].split("(")[0]
        return f"{n} strict-CAAX pheromone ORF{'s' if n != '1' else ''}"
    return r.replace("_", " ")


CAUSES = {
    "homothallism": "homothallism (a self-fertile species)",
    "hybrid_or_fusion": "hybrid or fused nuclei",
    "duplication": "a duplication or paralog",
    "mixed_culture_or_heterokaryon": "a mixed culture or heterokaryon",
    "assembly_artefact": "an assembly artefact",
}

ARRANGEMENT = {
    "same_locus": "at the same locus",
    "same_contig_distant": "on the same contig, more than 50 kb apart",
    "unlinked": "on different contigs",
}


# ---------------------------------------------------------------- formatting

def e(value) -> str:
    """Escaped text; None as an en dash."""
    return "–" if value is None or value == "" else escape(str(value))


def wbr(value) -> str:
    """Escaped identifier with break opportunities after `_` and `:`, so long ids wrap between tokens."""
    if value is None or value == "":
        return "–"
    return escape(str(value)).replace("_", "_<wbr>").replace(":", ":<wbr>")


def plain(value, table: dict) -> str:
    """A pipeline key in plain words: the mapped label, else the key with spaces."""
    if value is None:
        return "–"
    return escape(table.get(value, str(value).replace("_", " ")))


def num(n) -> str:
    return "–" if n is None else f"{int(n):,}"


def one(v) -> str:
    return "–" if v is None else f"{float(v):,.1f}"


def evalue(v) -> str:
    if v is None:
        return "–"
    v = float(v)
    return "&lt; 1e-180" if v == 0 else f"{v:.1e}"


def size(bp) -> str:
    if bp is None:
        return "–"
    if bp >= 1_000_000:
        return f"{bp / 1e6:,.1f} Mb"
    return f"{bp / 1000:,.1f} kb" if bp >= 1000 else f"{bp:,} bp"


def plural(n: int, word: str, many: str | None = None) -> str:
    return f"{n:,} {word if n == 1 else (many or word + 's')}"


def family_parts(fam: str | None) -> tuple[str, str]:
    phylum, _, name = (fam or "").partition(":")
    return (name or phylum, phylum if name else "")


def badge(conf: str | None) -> str:
    c = conf or "none"
    return f'<span class="badge conf-{e(c)}">{e(c)} confidence</span>'


def location(contig, start, end) -> str:
    return f'<span class="mono">{wbr(contig)}:<wbr>{num(start)}–{num(end)}</span>'


def _has_flank(call: dict) -> bool:
    return any(str(g.get("role", "")).startswith("flanking") for g in call.get("gene_evidence") or [])


def _locus_class_text(call: dict) -> str:
    cls = call.get("locus_class")
    if cls == "mat_locus" and not _has_flank(call):
        return "MAT genes found; no flanking gene identified."
    return LOCUS_CLASS.get(cls, (cls or "–").replace("_", " "))


def _ignored(call: dict) -> set[str]:
    return set(call.get("confidence_ignored_genes") or [])


# ---------------------------------------------------------------- result

def _gene_order(call: dict) -> str:
    """Genes on the call's main contig, left to right; core genes in bold, ignored cross-hits left out."""
    ign = _ignored(call)
    genes = sorted((g for g in call.get("gene_evidence") or []
                    if g.get("contig") == call.get("contig") and g.get("gene") not in ign),
                   key=lambda g: g.get("start") or 0)
    seen: list[str] = []
    for g in genes:
        name = f'<span class="nowrap">{e(g.get("gene"))}</span>'
        name = f"<b>{name}</b>" if g.get("role") == "core_MAT" else name
        if name not in seen:
            seen.append(name)
    return " · ".join(seen)


def _call_sentence(call: dict, show_family: bool) -> str:
    name, _ = family_parts(call.get("family"))
    cls = call.get("locus_class")
    if cls == "mat_locus":
        what = "complete locus" if _has_flank(call) else "MAT genes without flanking genes"
    else:
        what = {"idiomorph_gene_only": "MAT gene without flanking genes",
                "homothallic_candidate": "locus with genes of both idiomorphs"}.get(cls, "locus")
    span = (call.get("end") or 0) - (call.get("start") or 0) + 1
    fam = f" ({e(name)} family)" if show_family else ""
    order = _gene_order(call)
    return (f"<b>{e(call.get('idiomorph'))}</b>{fam}: {e(what)}, {e(call.get('confidence'))} confidence, "
            f"{size(span)} on <span class='mono'>{wbr(call.get('contig'))}</span>"
            + (f". Gene order: {order}." if order else "."))


def _two_idiomorph_block(t: dict) -> str:
    plus, minus = t.get("plus_call") or {}, t.get("minus_call") or {}
    supported = [c for c in t.get("supported_causes") or []]
    rest = [c for c in CAUSES if c not in supported]
    ev = t.get("evidence") or {}
    extra = []
    if ev.get("gc_difference_pct") is not None:
        extra.append(f"The two contigs differ in GC content by {ev['gc_difference_pct']} percentage points, which "
                     "suggests they may come from different genomes")
    if ev.get("weak_calls"):
        extra.append("The " + " and ".join(ev["weak_calls"]) + " call is the weaker one, so it may be a paralog "
                     "rather than a second idiomorph")
    return (
        f'<p class="verdict review">Two idiomorphs found ({e(plus.get("idiomorph"))} and {e(minus.get("idiomorph"))}): '
        'needs review</p>'
        f'<p>{e(plus.get("idiomorph"))} ({e(plus.get("confidence"))} confidence) and {e(minus.get("idiomorph"))} '
        f'({e(minus.get("confidence"))} confidence) were found {e(ARRANGEMENT.get(t.get("arrangement"), t.get("arrangement")))}. '
        + "".join(f"{e(x)}. " for x in extra)
        + "This is a flag, not a verdict on the strain's mating system. "
        + (f'Supported by the evidence: {e(", ".join(CAUSES.get(c, c) for c in supported))}. ' if supported
           else "No cause is favoured by the evidence. ")
        + (f'Not ruled out: {e(", ".join(CAUSES[c] for c in rest))}.' if rest else "")
        + "</p>")


def _result(doc: dict, phylum: str | None) -> str:
    calls = sorted(doc.get("detected") or [], key=lambda c: -_CONF_RANK.get(c.get("confidence"), 0))
    families = {c.get("family") for c in calls}
    show_family = len(families) > 1
    mode = doc.get("routing_mode")
    two = doc.get("two_idiomorphs") or []
    notes: list[tuple[str, str]] = []
    if mode == "not_searched":
        where = e(phylum) if phylum else "this lineage"
        body = (f'<p class="verdict none">Not searched</p>'
                f'<p>No curated MAT family covers {where}, so this genome was not searched. To search anyway, re-run '
                'with every family (<code>--exhaustive</code>) or with a chosen phylum (<code>--phylum</code>); '
                'calls made that way are unconfirmed.</p>')
    elif not calls:
        body = ('<p class="verdict none">No MAT locus called</p>'
                '<p>No family reached the evidence bar. That is weak evidence that the locus is absent: a fragmented '
                'assembly, or a lineage far from every reference, often explains it.</p>'
                + ('<p><b>Next:</b> check the withheld candidate loci below against a known locus, or type the '
                   'strain from its reads.</p>' if doc.get("suppressed_loci") else
                   '<p><b>Next:</b> type the strain from its reads, or check the assembly around the expected '
                   'flanking genes.</p>'))
    else:
        chips = "".join(
            f'<a class="chip" href="#locus-{i}"><span class="chip-id">{e(c.get("idiomorph"))}</span>'
            + (f'<span class="chip-fam">{e(family_parts(c.get("family"))[0])}</span>' if show_family else "")
            + f'{badge(c.get("confidence"))}</a>'
            for i, c in enumerate(calls))
        if two:
            lead = "".join(_two_idiomorph_block(t) for t in two)
            body = lead + f'<div class="chips small-chips">{chips}</div>'
        else:
            label = "Mating type" if len(calls) == 1 else "Mating-type loci"
            body = f'<p class="lead-label">{label}</p><div class="chips">{chips}</div>'
        body += "<ul class='sentences'>" + "".join(f"<li>{_call_sentence(c, show_family)}</li>" for c in calls) + "</ul>"
    if doc.get("routing_error"):
        notes.append(("warn", "The taxonomy lookup failed, so more families were searched than this genome needs. "
                              f"Re-run when the lookup works. ({doc['routing_error']})"))
    z = doc.get("zygosity")
    if z and z.get("status") in ("unknown", "unchecked"):
        notes.append(("info", f"Zygosity {z.get('status')}: {z.get('reason')}"))
    for g in doc.get("assembly_gap_at_locus") or []:
        notes.append(("info", f"Assembly gap where the {family_parts(g.get('family'))[0]} locus is expected: "
                              f"{g.get('contig')}:{g.get('start'):,}–{g.get('end'):,}, {g.get('n_bases'):,} N bases "
                              f"between {', '.join(g.get('anchors') or [])}. The assembly never resolved this "
                              "region, so a missing call here is not evidence of absence."))
    note_html = "".join(f'<li class="{lvl}">{e(t)}</li>' for lvl, t in notes)
    return f'''
<section class="result" aria-labelledby="result-h">
  <h2 id="result-h">Result</h2>
  {body}
  {f'<ul class="flags">{note_html}</ul>' if note_html else ""}
</section>'''


# ---------------------------------------------------------------- locus card

def _flags(call: dict) -> list[tuple[str, str]]:
    """(level, text) notes for one call; 'warn' items print with a "Check:" prefix."""
    out: list[tuple[str, str]] = []
    ver = call.get("verification")
    if ver:
        reason = str(ver.get("reason") or ver.get("status") or "")
        if "CAAX" in reason:
            out.append(("warn", "The pheromone-precursor gene was found only by a sequence-motif scan (CAAX), with "
                                "no reference match. Treat this call as unconfirmed."))
        else:
            out.append(("warn", f"Unconfirmed call: {reason.rstrip('.')}."))
    if call.get("idiomorph") == "undetermined":
        out.append(("warn", "Idiomorph undetermined: the evidence for two idiomorphs is tied or too close to call."))
    if call.get("fragmented"):
        out.append(("warn", "The locus is split across contigs or interrupted; confidence was lowered one step."))
    if call.get("idiomorph_unmodelled"):
        out.append(("warn", "The idiomorph rests on MAT-gene hits that no tool could model."))
    if call.get("span_exceeds_plausible_bound"):
        out.append(("info", "The locus is longer than this family's usual span. Large loci are real; check the figure."))
    if call.get("ambiguous_with"):
        out.append(("info", "The same genes also match: " + ", ".join(call["ambiguous_with"]) + "."))
    sl = call.get("split_locus")
    if sl:
        out.append(("info", "Locus split into parts: " + str(sl if isinstance(sl, str) else sl.get("status", "see the YAML report"))))
    if call.get("receptor_array_id"):
        out.append(("info", f"Part of pheromone-receptor array {call['receptor_array_id']} "
                            f"({plural(int(call.get('receptor_array_size') or 0), 'locus', 'loci')}, "
                            f"{call.get('receptor_array_support')}). Arrays hold mating and non-mating receptors, so "
                            "membership is not evidence of mating function."))
    for g in sorted(_ignored(call)):
        out.append(("info", f"{g} is a weak hit to the other idiomorph's gene and was not counted against the call."))
    return out


def _gene_rows(call: dict) -> str:
    ign = _ignored(call)
    rows = []
    for g in sorted(call.get("gene_evidence") or [], key=lambda g: (g.get("contig") or "", g.get("start") or 0)):
        status = g.get("status") or ""
        subs = [plain(g.get("method"), METHOD)]
        alt = g.get("alternate_model")
        if alt:
            subs.append(f"{plain(alt.get('method'), METHOD)}: {num(alt.get('start'))}–{num(alt.get('end'))}, "
                        f"{one(alt.get('identity'))}%")
        if g.get("frameshifts"):
            subs.append(f'<span class="warn-text">{plural(int(g["frameshifts"]), "frameshift")}</span>')
        if g.get("caax_motif"):
            subs.append(f"motif {e(g.get('caax_motif'))}, {num(g.get('orf_length_aa'))} aa ORF")
        role = g.get("role") or ""
        role_txt = e(ROLE.get(role, role.replace("_", " ") or "–"))
        cls = ""
        if g.get("gene") in ign:
            role_txt += '<span class="sub">other idiomorph, not counted</span>'
            cls = ' class="ignored"'
        cov = f'<span class="sub">cov {one(g.get("coverage"))}%</span>' if g.get("coverage") is not None else ""
        bits = (f'<span class="sub nowrap">{float(g["bitscore"]):,.0f} bits</span>'
                if g.get("bitscore") is not None else "")
        exons = g.get("exons")
        rows.append(
            f'<tr{cls}><th scope="row"><span class="swatch {_role_css(role)}" aria-hidden="true"></span>{e(g.get("gene"))}</th>'
            f'<td class="nowrap">{role_txt}</td>'
            f'<td><span class="status st-{e(status)}">{e(STATUS.get(status, status.replace("_", " ") or "–"))}</span>'
            + "".join(f'<span class="sub">{s}</span>' for s in subs) + "</td>"
            f'<td class="num">{one(g.get("identity"))}{cov}</td>'
            f'<td class="num">{evalue(g.get("evalue"))}{bits}</td>'
            f'<td class="mono num">{num(g.get("start"))}–{num(g.get("end"))} <span class="strand">{e(g.get("strand"))}</span></td>'
            f'<td class="num">{len(exons) if exons else "–"}</td>'
            f'<td class="mono small rec">{wbr(g.get("reference_record"))}</td></tr>'
        )
    for name in call.get("genes_missing") or []:
        rows.append(f'<tr class="missing"><th scope="row">{e(name)}</th><td colspan="7">expected, not found</td></tr>')
    for name in call.get("genes_not_searchable") or []:
        rows.append(f'<tr class="missing"><th scope="row">{e(name)}</th>'
                    '<td colspan="7">not searchable: no reference protein in the database</td></tr>')
    return "".join(rows)


def _role_css(role: str) -> str:
    return {"core_MAT": "r-core", "flanking_conserved": "r-flank", "flanking_variable": "r-var"}.get(role, "r-other")


def _idiomorph_evidence(call: dict) -> str:
    bits = []
    cands = call.get("idiomorph_candidates") or []
    clf = call.get("idiomorph_classifier")
    if cands:
        items = "".join(f'<li><b>{e(c.get("idiomorph"))}</b> <span class="num">{one(c.get("score"))}</span></li>'
                        for c in cands)
        gloss = ("HMM bitscore of the MAT protein against each idiomorph's profile" if clf
                 else "bitscore of the best core-gene hit for each idiomorph")
        bits.append(f'<div><h5>Idiomorph scores</h5><ul class="scores">{items}</ul>{_margin_text(call, cands)}'
                    f'<p class="small muted">Score: {gloss}.</p></div>')
    if clf:
        scores = "".join(f'<li><b>{e(k)}</b> <span class="num">{one(v)}</span></li>'
                         for k, v in (clf.get("scores") or {}).items())
        bits.append(
            f'<div><h5>Profile-HMM classifier</h5><ul class="scores">{scores}</ul>'
            f'<p class="small">Margin <b class="num">{one(clf.get("margin"))}</b> (a call needs {one(clf.get("min_margin"))} '
            f'or more). Scored: {e(", ".join(clf.get("genes_scored") or []))} ({plain(clf.get("classifier_input"), CLF_INPUT)}).</p></div>')
    res = call.get("idiomorph_resolutions") or []
    if res:
        rows = "".join(
            f'<tr><td class="nowrap">{e(r.get("winner"))}</td><td class="num">{one(r.get("winner_identity"))}</td>'
            f'<td class="nowrap">{e(r.get("loser"))}</td><td class="num">{one(r.get("loser_identity"))}</td>'
            f'<td class="num">{one((r.get("overlap_fraction") or 0) * 100)}</td><td>{plain(r.get("basis"), BASIS)}</td></tr>'
            for r in res)
        bits.append(
            '<div class="wide"><h5>Overlapping hits resolved</h5>'
            '<p class="small">Where genes of two idiomorphs hit the same place, the stronger hit was kept.</p>'
            '<div class="table-wrap"><table class="compact"><thead><tr><th>kept</th><th class="num">identity %</th>'
            '<th>set aside</th><th class="num">identity %</th><th class="num">overlap %</th><th>basis</th></tr></thead>'
            f'<tbody>{rows}</tbody></table></div></div>')
    if not bits:
        return ""
    return f'<div class="evidence">{"".join(bits)}</div>'


def _margin_text(call: dict, cands: list[dict]) -> str:
    """What `idiomorph_margin` measures depends on how the idiomorph was decided (pipeline: the vote's score gap;
    with one candidate, the narrowest overlap resolution in identity points; with a classifier, its HMM margin)."""
    margin = call.get("idiomorph_margin")
    res = call.get("idiomorph_resolutions") or []
    if len(cands) > 1 and margin is not None:
        return f'<p class="small">Margin over the next idiomorph: <b class="num">{one(margin)}</b></p>'
    if margin is not None and res:
        narrow = min(res, key=lambda r: (r.get("winner_identity") or 0) - (r.get("loser_identity") or 0))
        return (f'<p class="small">Margin over the set-aside <span class="nowrap">{e(narrow.get("loser"))}</span> hit: '
                f'<b class="num">{one(margin)}</b> identity points (the narrowest overlap resolution below).</p>')
    if margin is not None:
        return f'<p class="small">Margin: <b class="num">{one(margin)}</b></p>'
    return '<p class="small">No other idiomorph scored.</p>'


def _subloci(call: dict) -> str:
    subs = call.get("subloci") or []
    if not subs:
        return ""
    rows = "".join(
        f'<tr><td>{e(s.get("sublocus"))}</td><td>{e(s.get("idiomorph"))}</td>'
        f'<td>{location(s.get("contig"), s.get("start"), s.get("end"))}</td>'
        f'<td>{e(", ".join(s.get("genes") or []))}</td><td>{e(s.get("completeness"))}'
        + (f'<span class="sub">missing: {e(", ".join(s.get("genes_missing") or []))}</span>' if s.get("genes_missing") else "")
        + "</td></tr>" for s in subs)
    return ('<h4>Subloci</h4><div class="table-wrap"><table class="compact"><thead><tr><th>sublocus</th><th>allele</th>'
            f'<th>location</th><th>genes</th><th>completeness</th></tr></thead><tbody>{rows}</tbody></table></div>')


def _legend(keys: set[str]) -> str:
    items = [
        ("r-core", '<span class="swatch r-core"></span>core MAT gene'),
        ("r-flank", '<span class="swatch r-flank"></span>conserved flank'),
        ("r-var", '<span class="swatch r-var"></span>variable flank'),
        ("r-other", '<span class="swatch r-other"></span>other gene'),
        ("s-model", '<span class="key key-model"></span>gene model'),
        ("s-disagree", '<span class="key key-disagree"></span>tools disagree'),
        ("alt", '<span class="key key-alt"></span>alternate model'),
        ("s-hit", '<span class="key key-hit"></span>search hit only'),
    ]
    if not keys & {"s-disagree", "s-hit"}:
        keys = keys - {"s-model"}
    lis = "".join(f"<li>{html}</li>" for k, html in items if k in keys).replace(
        '"></span>', '" aria-hidden="true"></span>')
    return f'<ul class="legend" aria-label="Figure legend">{lis}</ul>'


def _locus_card(i: int, call: dict, show_family: bool) -> str:
    name, phylum = family_parts(call.get("family"))
    span = (call.get("end") or 0) - (call.get("start") or 0) + 1
    contigs: list[str] = []
    for g in call.get("gene_evidence") or []:
        if g.get("contig") and g["contig"] not in contigs:
            contigs.append(g["contig"])
    if not contigs:
        contigs = [call.get("contig")]
    label = f"{call.get('idiomorph')} locus"
    figs = []
    for j, contig in enumerate(contigs):
        part = f" (part {j + 1} of {len(contigs)})" if len(contigs) > 1 else ""
        figs.append(f'<figure class="fig"><div class="fig-scroll">{locus_svg(call, contig, f"fig{i}-{j}", label)}</div>'
                    f'<figcaption>Gene order on <span class="mono">{wbr(contig)}</span>{part}.'
                    '<span class="scroll-hint"> Scroll sideways to see the whole figure.</span></figcaption></figure>')
    edge = [s.get("contig_edge_distance") for s in call.get("segments") or [] if s.get("contig_edge_distance") is not None]
    flags = "".join(f'<li class="{lvl}">{e(text)}</li>' for lvl, text in _flags(call))
    n_hits = len(call.get("gene_evidence") or [])
    region = f'{call.get("contig")}:{call.get("start")}-{call.get("end")}'
    facts = [
        ("Location", location(call.get("contig"), call.get("start"), call.get("end")) + f" · {size(span)}"),
        ("Locus", e(_locus_class_text(call))),
        ("Genes modelled", f'{e(call.get("polished_genes"))} of {plural(n_hits, "gene hit")}'),
        ("Nearest contig end", size(min(edge)) if edge else "–"),
        ("Region (samtools)", f'<code class="region">{e(region)}</code>'),
    ]
    core = call.get("core_span")
    if core:
        beyond = core.get("beyond_core_bp") or 0
        note = (f' · <span class="warn-text">the span runs {size(beyond)} beyond these genes</span>' if beyond > 0 else "")
        facts.insert(1, ("Own genes", location(call.get("contig"), core.get("start"), core.get("end"))
                         + f" · {size((core.get('end') or 0) - (core.get('start') or 0) + 1)}{note}"))
    if call.get("detection_pass") and call.get("detection_pass") != "strict":
        facts.append(("Detection pass", e(call.get("detection_pass"))))
    if call.get("idiomorph_class") and call.get("idiomorph_class") != call.get("idiomorph"):
        facts.insert(1, ("Cross-lineage class", e(call.get("idiomorph_class"))))
    dl = "".join(f"<div><dt>{k}</dt><dd>{v}</dd></div>" for k, v in facts)
    keys = figure_keys(call.get("gene_evidence") or [])
    return f'''
<section class="card locus" id="locus-{i}" aria-labelledby="locus-{i}-h">
  <div class="locus-top">
  <header class="card-head">
    <h3 id="locus-{i}-h"><span class="idiomorph">{e(call.get("idiomorph"))}</span>
      <span class="fam">{e(name)} family <span class="muted">· {e(phylum)}</span></span></h3>
    {badge(call.get("confidence"))}
  </header>
  {"".join(figs)}
  {_legend(keys)}
  {f'<ul class="flags">{flags}</ul>' if flags else ""}
  <dl class="facts">{dl}</dl>
  </div>
  {_idiomorph_evidence(call)}
  {_subloci(call)}
  <h4>Gene evidence</h4>
  <div class="table-wrap"><table class="genes">
    <thead><tr><th scope="col">Gene</th><th scope="col">Role</th><th scope="col">Model</th>
    <th scope="col" class="num">Identity %</th><th scope="col" class="num">E-value</th>
    <th scope="col" class="num">Position (strand)</th><th scope="col" class="num">Exons</th>
    <th scope="col">Reference record</th></tr></thead>
    <tbody>{_gene_rows(call)}</tbody>
  </table></div>
</section>'''


# ---------------------------------------------------------------- other sections

def _search_summary(doc: dict) -> str:
    mode = doc.get("routing_mode")
    label, text = ROUTING.get(mode, (mode or "–", ""))
    nd = doc.get("not_detected") or []
    nd_rows = "".join(
        f"<tr><td>{e(n.get('family'))}</td><td>{e(n.get('reason'))}</td>"
        f"<td class='num'>{one((n.get('best_fraction_found') or 0) * 100) if n.get('best_fraction_found') is not None else '–'}</td>"
        f"<td>{e(', '.join(n.get('genes_found') or []))}</td></tr>" for n in nd)
    gc = doc.get("genetic_code")
    return f'''
<section class="card" id="search" aria-labelledby="search-h">
  <h2 id="search-h">What was searched</h2>
  <dl class="facts">
    <div class="wide-fact"><dt>Routing</dt><dd><b>{e(label)}.</b> {e(text)}</dd></div>
    <div><dt>Families searched</dt><dd>{e(", ".join(doc.get("families_attempted") or []) or "none")}</dd></div>
    <div><dt>Genetic code</dt><dd>{e(gc)}{" · " + e(doc.get("genetic_code_error")) if doc.get("genetic_code_error") else ""}</dd></div>
    <div><dt>Taxonomy</dt><dd>{e(doc.get("taxonomy_source"))}</dd></div>
  </dl>
  {f"""<h4>Families searched but not called</h4><div class="table-wrap"><table class="compact"><thead><tr>
  <th>family</th><th>reason</th><th class="num">best % of genes found</th><th>genes found</th></tr></thead>
  <tbody>{nd_rows}</tbody></table></div>""" if nd else ""}
  {_receptor_arrays(doc)}
</section>'''


def _receptor_arrays(doc: dict) -> str:
    arrays = doc.get("receptor_arrays") or []
    if not arrays:
        return ""
    rows = "".join(
        f"<tr><td>{e(a.get('receptor_array_id'))}</td><td>{location(a.get('contig'), a.get('start'), a.get('end'))}</td>"
        f"<td class='num'>{e(a.get('receptor_array_size'))}</td><td class='num'>{e(a.get('calls'))}</td>"
        f"<td>{e(a.get('receptor_array_support'))}</td>"
        f"<td class='small'>{e('; '.join(array_reason(r) for r in a.get('receptor_array_support_reasons') or []))}</td></tr>"
        for a in arrays)
    return (f'<h4>Pheromone-receptor arrays</h4><p class="small">{e(doc.get("receptor_arrays_note"))}</p>'
            '<div class="table-wrap"><table class="compact"><thead><tr><th>array</th><th>location</th>'
            '<th class="num">loci</th><th class="num">calls</th><th>support</th><th>reasons</th></tr></thead>'
            f'<tbody>{rows}</tbody></table></div>')


def _withheld(doc: dict) -> str:
    loci = doc.get("suppressed_loci") or []
    if not loci:
        return ""
    counts: dict[str, int] = {}
    for s in loci:
        counts[s.get("withheld_reason")] = counts.get(s.get("withheld_reason"), 0) + 1
    summary = "; ".join(f"{WITHHELD.get(r, str(r).replace('_', ' '))} ({n})"
                        for r, n in sorted(counts.items(), key=lambda kv: -kv[1]))
    rows = "".join(
        f"<tr><td>{e(family_parts(s.get('family'))[0])}</td><td>{location(s.get('contig'), s.get('start'), s.get('end'))}</td>"
        f"<td>{e(s.get('idiomorph'))}</td><td>{e(', '.join(s.get('genes_found') or []))}</td>"
        f"<td class='num'>{e(s.get('polished_genes'))}</td><td>{plain(s.get('withheld_reason'), WITHHELD)}</td></tr>"
        for s in loci)
    n = len(loci)
    return f'''
<section class="card" id="withheld" aria-labelledby="withheld-h">
  <h2 id="withheld-h">Withheld candidate loci <span class="count">{n}</span></h2>
  <p>Places with some MAT-gene evidence that did not pass a quality gate. Most genomes have several. They are listed so
  a known locus can be checked against them, not because they are likely calls.</p>
  <p class="small">Reasons: {e(summary)}.</p>
  <details><summary>Show the {plural(n, "withheld locus", "withheld loci")}</summary>
  <h4 class="print-only">Withheld loci ({n})</h4>
  <div class="table-wrap"><table class="compact"><thead><tr><th>family</th><th>location</th><th>idiomorph</th>
  <th>genes</th><th class="num">modelled</th><th>reason</th></tr></thead><tbody>{rows}</tbody></table></div>
  </details>
</section>'''


def _glossary(has_calls: bool) -> str:
    items = []
    if has_calls:
        items += [
            ("High confidence", "Every expected core MAT gene was found and modelled, with a conserved flanking gene "
                                "where the family has one, on one contig."),
            ("Medium confidence", "The core genes were found, but a gene could not be modelled, no conserved flank "
                                  "was found, or the locus is split (one step down)."),
            ("Low confidence", "A single isolated gene hit, or no core MAT gene."),
            ("Model status", "Gene model: the reference protein was aligned to the genome to build the gene "
                             "(exonerate, miniprot). Models differ: the two tools disagree on the structure; the "
                             "alternate model is shown. Hit only: a search hit with no model, weaker evidence."),
            ("Identity", "Percent amino-acid identity between the gene (or hit) and the closest reference protein."),
            ("Reference record", "The curated database record whose protein matched best. Each record's metadata "
                                 "cites the publication or deposit it comes from."),
        ]
    items.insert(0, ("Idiomorph", "The form of the mating-type locus a strain carries (for example MAT1-1 or MAT1-2, "
                                  "Plus or Minus). Basidiomycota HD and PR alleles are named after the reference they "
                                  "match."))
    dl = "".join(f"<div><dt>{k}</dt><dd>{v}</dd></div>" for k, v in items)
    return f'''
<section class="card glossary-card" id="glossary" aria-labelledby="glossary-h">
  <h2 id="glossary-h">How to read this report</h2>
  <dl class="glossary">{dl}</dl>
  <h4>Limits</h4>
  <ul class="plain">
    <li>MATPredict reads an assembly; it is not a mating assay. Two idiomorphs in one assembly can come from a
    homothallic species, a mixed culture or a merged diploid assembly.</li>
    <li>Fragmented assemblies can split a locus or lose its flanking genes; an uncalled genome is weak evidence of absence.</li>
    <li>Lineages with no curated record are searched with the whole phylum's references; family labels there are less reliable.</li>
  </ul>
</section>'''


def _provenance(doc: dict, files: list[str]) -> str:
    run = doc.get("run")
    files_txt = e(", ".join(files)) if files else "–"
    if not run:
        return f'''
<section class="card" id="provenance" aria-labelledby="prov-h">
  <h2 id="prov-h">Provenance</h2>
  <p>This report predates provenance recording: the run did not record its input file, MATPredict version, database
  or parameters. Re-run <code>matpredict detect</code> to record them.</p>
  <dl class="facts"><div><dt>Taxonomy</dt><dd>{e(doc.get("taxonomy_source"))}</dd></div>
  <div><dt>Files</dt><dd>{files_txt}</dd></div></dl>
</section>'''
    g = run.get("genome") or {}
    db = run.get("database") or {}
    params = run.get("parameters") or {}
    rows = [
        ("Genome file", wbr(g.get("file"))),
        ("Assembly", f'{plural(int(g["contigs"]), "contig")} · {size(g.get("length_bp"))} · N50 {size(g.get("n50"))}'
                     if g.get("length_bp") else e(g.get("error"))),
        ("Genome SHA-256", f'<span class="mono small break">{e(g.get("sha256"))}</span>' if g.get("sha256") else None),
        ("MATPredict", e(run.get("matpredict_version"))),
        ("Database", f'<span class="mono small">{e((db.get("content_sha256") or "")[:16])}</span> · '
                     f'{plural(int(db.get("records") or 0), "record")}' if db else None),
        ("Taxonomy", e(run.get("taxonomy_source") or doc.get("taxonomy_source"))),
        ("Parameters", ", <wbr>".join(f'<span class="nowrap">{wbr(k)}={e(v)}</span>'
                                      for k, v in params.items() if v not in (None, False, "")) or None),
        ("Run started", e(run.get("started"))),
        ("Run time", f'{run["wall_seconds"]:,} s' if run.get("wall_seconds") is not None else None),
        ("Files", files_txt),
    ]
    dl = "".join(f"<div><dt>{k}</dt><dd>{v}</dd></div>" for k, v in rows if v and v != "–")
    return f'''
<section class="card" id="provenance" aria-labelledby="prov-h">
  <h2 id="prov-h">Provenance</h2>
  <dl class="facts">{dl}</dl>
  <p class="small">Cite MATPredict as given in CITATION.cff at github.com/stajichlab/MATPredict. The reference records
  named in the gene tables are curated database entries; each record's metadata cites its source publication.</p>
</section>'''


# ---------------------------------------------------------------- page

def render_genome_report(doc: dict, *, sample: str | None = None, sample_from_folder: bool = False,
                         files: list[str] | None = None, print_mode: bool = False,
                         generated: str | None = None) -> str:
    """The full HTML page for one `detection_report.yaml` document.

    `print_mode` opens every collapsed section (PDF engines print a closed
    <details> closed). `generated` fixes the timestamp (tests)."""
    run = doc.get("run") or {}
    sample = sample or run.get("sample") or "unnamed sample"
    organism, taxid, phylum = run.get("organism"), run.get("taxid"), run.get("phylum")
    generated = generated or datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    calls = sorted(doc.get("detected") or [], key=lambda c: -_CONF_RANK.get(c.get("confidence"), 0))
    show_family = len({c.get("family") for c in calls}) > 1
    loci_html = "".join(_locus_card(i, c, show_family) for i, c in enumerate(calls))
    sub = " · ".join(x for x in [
        f"<i>{e(organism)}</i>" if organism else None,
        f"taxid {e(taxid)}" if taxid else None,
        e(phylum) if phylum else None,
    ] if x) or '<span class="muted">organism not recorded</span>'
    run_v = run.get("matpredict_version")
    meta = [f"Run by MATPredict {e(run_v)}" if run_v else "Run version not recorded"]
    if run.get("database"):
        meta.append(f'database <span class="mono">{e((run["database"].get("content_sha256") or "")[:12])}</span>')
    if run.get("taxonomy_source"):
        meta.append(e(run["taxonomy_source"]))
    meta.append(f"rendered {e(generated)} by MATPredict {e(__version__)}")
    folder_note = ' <span class="muted small">(name from folder)</span>' if sample_from_folder else ""
    body = f'''
<header class="masthead"><div class="mast-inner">
  <div class="mast-top"><div class="brand">MATPredict · mating-type report</div>
  <button type="button" class="print-btn" data-print>Save as PDF</button></div>
  <h1>{e(sample)}{folder_note}</h1>
  <p class="subtitle">{sub}</p>
  <div class="meta">{" · ".join(meta)}</div>
</div></header>
<main>
  {_result(doc, phylum)}
  {f'<h2 class="section-h">Called loci <span class="count">{len(calls)}</span></h2>' if calls else ""}
  {loci_html}
  {_search_summary(doc)}
  {_withheld(doc)}
  {_provenance(doc, files or [])}
  {_glossary(bool(calls))}
</main>
<footer class="foot">{e(sample)} · MATPredict mating-type report</footer>'''
    if print_mode:
        body = body.replace("<details>", "<details open>")
    title = f"{sample} – MATPredict report"
    return (
        "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">"
        f"<meta name=\"generator\" content=\"MATPredict {e(__version__)}\">"
        f"<title>{e(title)}</title><style>{_css(sample)}</style></head>"
        f"<body>{body}<script>{_JS}</script></body></html>\n"
    )


_JS = """
(function(){
  var b=document.querySelector('[data-print]');
  document.querySelectorAll('.fig-scroll').forEach(function(box){
    var svg=box.querySelector('svg'); var cx=svg&&parseFloat(svg.getAttribute('data-core-x'));
    if(box.scrollWidth>box.clientWidth&&cx){box.scrollLeft=Math.max(0,cx*svg.clientWidth/960-40);}
  });
  if(b){b.addEventListener('click',function(){window.print();});}
  var opened=[];
  window.addEventListener('beforeprint',function(){
    document.querySelectorAll('details:not([open])').forEach(function(d){d.open=true;opened.push(d);});
  });
  window.addEventListener('afterprint',function(){opened.forEach(function(d){d.open=false;});opened=[];});
})();
"""


def _css_string(s: str) -> str:
    """A CSS string literal; `<` and `>` as CSS escapes so a name cannot close the <style> element."""
    s = s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ")
    return '"' + s.replace("<", "\\3C ").replace(">", "\\3E ") + '"'


_LIGHT = """
  --bg:#ffffff; --surface:#ffffff; --surface-2:#f4f5f7; --ink:#1b1f24; --ink-2:#424a54; --muted:#5f6873;
  --line:#dde1e6; --line-2:#c4cad1; --accent:#1f5f99;
  --core:#b8471a; --flank:#2a6fb0; --other:#6b7480;
  --high:#1d7a46; --high-bg:#e3f3ea; --medium:#7a4f00; --medium-bg:#fbf0d9; --low:#4f5864; --low-bg:#eceef1;
  --warn:#8a4b00; --warn-bg:#fff6e5; --info-bg:#eef4fa;"""
_DARK = """
  --bg:#111418; --surface:#171b20; --surface-2:#1e232a; --ink:#e6e9ed; --ink-2:#c3c9d0; --muted:#9aa3ad;
  --line:#2c333b; --line-2:#3a434d; --accent:#7fb6ea;
  --core:#f08a4b; --flank:#5ea3e6; --other:#9aa3ad;
  --high:#6fd39b; --high-bg:#173527; --medium:#f0c060; --medium-bg:#3a2e12; --low:#b4bcc5; --low-bg:#262c33;
  --warn:#f0b866; --warn-bg:#33270f; --info-bg:#172634;"""


def _css(sample: str) -> str:
    return (":root{" + _LIGHT + """
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, "Liberation Mono", monospace;
  --sans: system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
}
@media (prefers-color-scheme: dark){ :root:not([data-theme="light"]){""" + _DARK + """} }
:root[data-theme="dark"]{""" + _DARK + """}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 var(--sans);font-variant-numeric:tabular-nums}
.masthead,main,.foot{max-width:1040px;margin:0 auto;padding:0 16px}
.mast-inner{padding:28px 0 10px;border-bottom:1px solid var(--line)}
.brand{font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);font-weight:600}
h1{font-size:28px;line-height:1.2;margin:6px 0 2px;overflow-wrap:anywhere}
.subtitle{margin:0;color:var(--ink-2);font-size:16px}
.mast-top{display:flex;justify-content:space-between;align-items:center;gap:12px}
.meta{margin:8px 0 2px;color:var(--muted);font-size:13px}
.print-btn{font:inherit;font-size:13px;padding:5px 12px;border:1px solid var(--line-2);border-radius:6px;
  background:var(--surface);color:var(--ink);cursor:pointer}
.print-btn:hover{border-color:var(--accent);color:var(--accent)}
:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
h2{font-size:20px;margin:0 0 12px}
h3{font-size:18px;margin:0}
h4{font-size:13px;margin:20px 0 8px;color:var(--ink-2);text-transform:uppercase;letter-spacing:.05em}
h5{font-size:13px;margin:0 0 6px;color:var(--ink-2)}
p{max-width:75ch}
code{font-family:var(--mono);font-size:.92em}
.section-h{margin:28px 0 12px}
.count{display:inline-block;min-width:1.6em;padding:0 6px;border-radius:10px;background:var(--surface-2);
  color:var(--ink-2);font-size:13px;text-align:center;vertical-align:middle}
.result{margin:20px 0 8px;padding:20px;border:1px solid var(--line);border-left:4px solid var(--accent);
  border-radius:8px;background:var(--surface)}
.verdict{font-size:24px;font-weight:700;margin:0 0 8px;line-height:1.25}
.verdict.none{color:var(--ink-2)}
.verdict.review{color:var(--warn)}
.lead-label{margin:0 0 6px;font-size:13px;text-transform:uppercase;letter-spacing:.05em;color:var(--muted)}
.chips{display:flex;flex-wrap:wrap;gap:10px;margin:0 0 10px}
.chip{display:inline-flex;align-items:center;gap:10px;padding:8px 14px;border:1px solid var(--line-2);border-radius:8px;
  text-decoration:none;color:var(--ink);background:var(--surface-2)}
.chip:hover{border-color:var(--accent)}
.chip-id{font-size:26px;font-weight:700;line-height:1.1;white-space:nowrap}
.chip-fam{color:var(--ink-2);font-size:14px;white-space:nowrap}
.small-chips .chip{padding:4px 10px}
.small-chips .chip-id{font-size:16px}
.sentences{margin:6px 0 0;padding-left:20px}
.sentences li{margin:3px 0}
.badge{display:inline-block;padding:1px 9px;border-radius:999px;font-size:12px;font-weight:600;
  border:1px solid currentColor;white-space:nowrap}
.conf-high{color:var(--high);background:var(--high-bg)}
.conf-medium{color:var(--medium);background:var(--medium-bg)}
.conf-low,.conf-none{color:var(--low);background:var(--low-bg)}
.flags{list-style:none;padding:0;margin:12px 0}
.flags li{padding:8px 12px;border-radius:6px;margin:6px 0;font-size:14px}
.flags li.warn{background:var(--warn-bg);color:var(--ink);border-left:3px solid var(--warn)}
.flags li.warn::before{content:"Check: ";font-weight:700;color:var(--warn)}
.flags li.info{background:var(--info-bg);border-left:3px solid var(--accent)}
.card{margin:16px 0;padding:20px;border:1px solid var(--line);border-radius:8px;background:var(--surface)}
.card-head{display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin-bottom:10px}
.idiomorph{font-weight:700;margin-right:6px}
.fam{font-weight:500;color:var(--ink-2)}
.muted{color:var(--muted);font-weight:400}
.fig{margin:0 0 4px}
.fig-scroll{overflow-x:auto}
.locus-fig{width:100%;height:auto;display:block}
figcaption{font-size:12px;color:var(--muted);margin-top:2px}
.scroll-hint{display:none}
.locus-fig .axis,.locus-fig .tick{stroke:var(--line-2)}
.locus-fig .tick-label,.locus-fig .end-label{fill:var(--muted)}
.locus-fig .contig-end{stroke:var(--ink)}
.locus-fig .bracket{stroke:var(--ink-2)} .locus-fig .bracket-label{fill:var(--ink-2)}
.locus-fig .intron{stroke:var(--ink-2)}
.locus-fig .r-core{--c:var(--core)} .locus-fig .r-flank,.locus-fig .r-var{--c:var(--flank)} .locus-fig .r-other{--c:var(--other)}
.locus-fig .exon{stroke:var(--c);fill:var(--c)}
.locus-fig .r-var:not(.s-hit) .exon{fill:var(--surface)}
.locus-fig .s-disagree .exon{stroke:var(--ink)}
.locus-fig .gene-label{fill:var(--ink)} .locus-fig .s-hit .gene-label{fill:var(--ink-2)}
.locus-fig .alt{fill:var(--ink-2)}
.legend{list-style:none;display:flex;flex-wrap:wrap;gap:6px 18px;padding:0;margin:8px 0 4px;font-size:12px;color:var(--ink-2)}
.legend li{display:flex;align-items:center;gap:6px;white-space:nowrap}
.swatch{display:inline-block;width:12px;height:12px;border-radius:2px;vertical-align:-1px;margin-right:6px;flex:none}
.legend .swatch{margin-right:0}
.r-core{background:var(--core)} .r-flank{background:var(--flank)} .r-other{background:var(--other)}
.swatch.r-var{background:var(--surface);border:2px solid var(--flank)}
.key{display:inline-block;width:22px;height:11px;border-radius:2px;flex:none}
.key-model{background:var(--surface);border:1.5px solid var(--ink-2)}
.key-disagree{background:var(--other);border:1.5px dashed var(--ink)}
.key-alt{height:3px;background:var(--ink-2)}
.key-hit{border:1.5px dotted var(--ink-2);background:var(--line)}
.facts{display:flex;flex-wrap:wrap;gap:10px 28px;margin:14px 0}
.facts>div{flex:1 1 200px;min-width:0}
.facts>.wide-fact{flex-basis:100%}
dt{font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
dd{margin:2px 0 0;overflow-wrap:anywhere}
.region{user-select:all;background:var(--surface-2);padding:1px 5px;border-radius:4px;white-space:nowrap;
  display:inline-block;max-width:100%;overflow-x:auto;vertical-align:bottom}
.evidence{display:flex;flex-wrap:wrap;gap:16px 40px;margin:12px 0;padding:12px 14px;background:var(--surface-2);border-radius:6px}
.evidence>div{min-width:0}
.evidence .wide{flex-basis:100%}
.evidence .muted{color:var(--ink-2)}
.scores{list-style:none;padding:0;margin:0;display:flex;gap:16px;flex-wrap:wrap}
.table-wrap{overflow-x:auto;-webkit-overflow-scrolling:touch;border:1px solid var(--line);border-radius:6px}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{padding:6px 10px;text-align:left;vertical-align:top;border-bottom:1px solid var(--line)}
thead th{background:var(--surface-2);font-weight:600;color:var(--ink-2);white-space:nowrap;font-size:12px}
tbody tr:last-child th,tbody tr:last-child td{border-bottom:0}
tbody th{font-weight:600;white-space:nowrap}
.num{text-align:right;white-space:nowrap} .nowrap{white-space:nowrap}
.mono{font-family:var(--mono);font-size:.92em}
.small{font-size:12px} .break{word-break:break-all}
.sub{display:block;font-size:11px;color:var(--muted);font-weight:400;font-style:normal;white-space:normal}
.num .sub{text-align:right}
.strand{color:var(--muted);font-family:var(--sans)}
td.rec{min-width:18ch;overflow-wrap:anywhere}
.status{font-weight:600;white-space:nowrap}
.st-polished_disagree,.warn-text{color:var(--warn)}
.st-unpolished,.st-not_polish_candidate{color:var(--muted);font-style:italic}
tr.missing th,tr.missing td,tr.ignored th,tr.ignored td{color:var(--muted)}
tr.missing td{font-style:italic}
details summary{cursor:pointer;color:var(--accent);margin:8px 0;font-weight:500}
.print-only{display:none}
.glossary{display:flex;flex-wrap:wrap;gap:12px 28px;margin:0}
.glossary>div{flex:1 1 280px}
.glossary dt{text-transform:none;letter-spacing:0;font-size:14px;font-weight:600;color:var(--ink)}
.glossary dd{color:var(--ink-2);font-size:14px}
.plain{padding-left:20px;color:var(--ink-2);font-size:14px}
.foot{color:var(--muted);font-size:12px;padding-top:12px;padding-bottom:32px}
a{color:var(--accent)}
@media (max-width:600px){
  body{font-size:14px}
  h1{font-size:22px}
  .card,.result{padding:14px}
  .chip-id{font-size:22px}
  .locus-fig{min-width:960px}
  .scroll-hint{display:inline}
}
@page{margin:16mm 14mm 18mm;
  @top-left{content:""" + _css_string(sample) + """;font:8pt sans-serif;color:#5f6873}
  @top-center{content:string(locus, start);font:8pt sans-serif;color:#5f6873}
  @top-right{content:"MATPredict mating-type report";font:8pt sans-serif;color:#5f6873}
  @bottom-right{content:"Page " counter(page) " of " counter(pages);font:8pt sans-serif;color:#5f6873}}
@media print{
  :root,:root[data-theme="dark"]{""" + _LIGHT + """}
  body{font-size:9.5pt;-webkit-print-color-adjust:exact;print-color-adjust:exact}
  .masthead,main{max-width:none;padding:0}
  .mast-inner{padding-top:0}
  .print-btn,.foot,.scroll-hint,details summary{display:none}
  .print-only{display:block}
  .locus-fig{min-width:0}
  .region{white-space:normal;overflow:visible;overflow-wrap:anywhere;display:inline}
  .card,.result{padding:10px 12px;margin:8px 0}
  .card-head h3{string-set:locus content()}
  #search h2,#withheld h2,#provenance h2,#glossary h2,.section-h,#result-h{string-set:locus ""}
  .locus-top,.facts>div,.glossary>div,.legend,.flags li,tr,.chips,.result,.scores{break-inside:avoid}
  h2,h3,h4,.card-head{break-after:avoid}
  .table-wrap{overflow:visible;border:0}
  table{font-size:8pt}
  table .small,table .mono,table .sub{font-size:inherit}
  table .sub{font-size:7pt}
  th,td{padding:3px 5px}
  /* A table must not run past its card on paper: let headers, numbers and record names wrap. */
  thead th,tbody th,.num,.status{white-space:normal}
  th,td{overflow-wrap:anywhere}
  td.rec{min-width:0}
  thead{display:table-header-group}
  .glossary-card{break-inside:avoid}
  .glossary>div{flex-basis:30%}
  .glossary dt,.glossary dd,.plain{font-size:8pt}
  a{color:inherit;text-decoration:none}
}
""")
