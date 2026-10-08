"""Gene-order figure of one called locus, as inline SVG.

Written here rather than with pyGenomeViz so the figure is vector in the PDF,
a few kB, and needs no plotting dependency in `detect`'s environment.
Colours are written as SVG presentation attributes (light theme), so engines
that ignore page CSS inside SVG (WeasyPrint) draw them; the page's CSS class
rules win over attributes and re-theme the figure on screen (dark mode).

Encoding (the page's legend shows only the keys a figure uses):
- fill = gene role; the roles differ in lightness as well as hue, so they
  survive greyscale printing: core MAT dark orange and drawn taller,
  conserved flank mid blue, variable flank white with a blue outline, other
  grey;
- outline = model status: solid (a gene model), dashed dark outline (the two
  modelling tools disagree; the alternate model is a thin bar under the
  gene), pale fill with a dashed outline (a search hit only, no model);
- arrowhead = strand; blocks joined by a line = exons and introns;
- a bracket under the axis = the reported locus; a cap or an arrow at the
  axis end = where the contig ends.
Genes and their labels are packed into lanes by their drawn width, so no
label is drawn over another.
"""
from __future__ import annotations

from html import escape

WIDTH = 960
MARGIN_X = 28
LANE_H = 36
GENE_H = 14
CORE_H = 18
LABEL_H = 16
AXIS_Y = 22
TRACK_TOP = 52
CHAR_PX = 7.0

ROLE_FILL = {"r-core": "#b8471a", "r-flank": "#2a6fb0", "r-var": "#ffffff", "r-other": "#6b7480"}
ROLE_STROKE = {"r-core": "#b8471a", "r-flank": "#2a6fb0", "r-var": "#2a6fb0", "r-other": "#6b7480"}
INK, INK_2, MUTED, LINE_2 = "#1b1f24", "#4a525c", "#5f6873", "#9aa3ad"
FONT = 'font-family="Helvetica, Arial, sans-serif"'
_ROLE_WORDS = {"core_MAT": "core MAT gene", "flanking_conserved": "conserved flank", "flanking_variable": "variable flank"}
_STATUS_WORDS = {"s-model": "gene model", "s-disagree": "the two modelling tools disagree", "s-hit": "search hit only"}


def role_class(role: str | None) -> str:
    return {"core_MAT": "r-core", "flanking_conserved": "r-flank",
            "flanking_variable": "r-var"}.get(role or "", "r-other")


def status_class(gene: dict) -> str:
    status = gene.get("status") or ""
    if status == "polished_disagree":
        return "s-disagree"
    if status.startswith("polished"):
        return "s-model"
    return "s-hit"


def _exons(gene: dict) -> list[tuple[int, int]]:
    out = [(int(e["start"]), int(e["end"])) for e in (gene.get("exons") or []) if e]
    return sorted(out) or [(int(gene["start"]), int(gene["end"]))]


def assign_lanes(extents: list[tuple[float, float]], modelled: list[bool], gap_px: float = 6.0) -> list[int]:
    """Greedy packing of pixel extents (gene body plus label); modelled genes first, so they keep lane 0."""
    order = sorted(range(len(extents)), key=lambda i: (not modelled[i], extents[i][0]))
    lane_ends: list[float] = []
    lanes = [0] * len(extents)
    for i in order:
        lo, hi = extents[i]
        for lane, end in enumerate(lane_ends):
            if lo > end + gap_px:
                lanes[i] = lane
                lane_ends[lane] = hi
                break
        else:
            lanes[i] = len(lane_ends)
            lane_ends.append(hi)
    return lanes


def _nice_step(span: int) -> int:
    for step in (100, 200, 500, 1000, 2000, 5000, 10_000, 20_000, 50_000, 100_000, 200_000, 500_000):
        if span / step <= 7:
            return step
    return 1_000_000


def _fmt_kb(bp: int, step: int) -> str:
    return f"{bp / 1000:,.0f} kb" if step >= 1000 else f"{bp / 1000:,.1f} kb"


def contig_end(segment: dict) -> tuple[str, int] | None:
    """('left'|'right', position of the contig end) from a segment's `contig_edge_distance`.

    The distance is to the nearer contig end; when it equals start - 1 the
    nearer end is the contig's first base, otherwise it lies `distance` bases
    past the segment end."""
    d = segment.get("contig_edge_distance")
    if d is None:
        return None
    if int(segment["start"]) - 1 == int(d):
        return "left", 1
    return "right", int(segment["end"]) + int(d)


def figure_keys(genes: list[dict]) -> set[str]:
    """Legend keys this figure uses: role classes, status classes, 'alt'."""
    keys: set[str] = set()
    for g in genes:
        keys.add(role_class(g.get("role")))
        sc = status_class(g)
        keys.add(sc)
        if sc == "s-disagree" and g.get("alternate_model"):
            keys.add("alt")
    return keys


def locus_svg(call: dict, contig: str, fig_id: str, label: str) -> str:
    """SVG of the genes of `call` on `contig` (one figure per contig of a split locus)."""
    genes = [g for g in (call.get("gene_evidence") or []) if g.get("contig") == contig
             and g.get("start") is not None and g.get("end") is not None]
    segs = [s for s in (call.get("segments") or []) if s.get("contig") == contig]
    lo_bp = min([g["start"] for g in genes] + [s["start"] for s in segs] or [call["start"]])
    hi_bp = max([g["end"] for g in genes] + [s["end"] for s in segs] or [call["end"]])
    pad = max(300, int((hi_bp - lo_bp) * 0.05))
    lo, hi = max(1, lo_bp - pad), hi_bp + pad
    ends = [c for c in (contig_end(s) for s in segs) if c]
    for side, pos in ends:  # show a nearby contig end inside the window
        if side == "left" and lo - pos < pad * 3:
            lo = pos
        if side == "right" and pos - hi < pad * 3:
            hi = pos
    span = max(hi - lo, 1)
    plot_w = WIDTH - 2 * MARGIN_X

    def x(bp: float) -> float:
        return MARGIN_X + (bp - lo) / span * plot_w

    extents, modelled = [], []
    for g in genes:
        x0, x1 = x(g["start"]), x(g["end"])
        half = (len(str(g["gene"])) * CHAR_PX + 4) / 2
        c = (x0 + x1) / 2
        extents.append((min(x0, c - half), max(x1, c + half)))
        modelled.append(status_class(g) != "s-hit")
    lanes = assign_lanes(extents, modelled)
    n_lanes = max(lanes, default=0) + 1
    height = TRACK_TOP + n_lanes * (LANE_H + LABEL_H) - LABEL_H + 10

    ordered = sorted(genes, key=lambda g: g["start"])
    desc = (f"{label}: genes from left to right on {contig}, {lo_bp:,} to {hi_bp:,}: "
            + "; ".join(f"{g['gene']} ({_ROLE_WORDS.get(g.get('role'), 'other gene')}, "
                        f"{'minus' if g.get('strand') == '-' else 'plus'} strand)" for g in ordered) + ".")
    core = [x(g["start"]) for g in ordered if role_class(g.get("role")) == "r-core"]
    core_attr = f' data-core-x="{core[0]:.0f}"' if core else ""
    parts = [
        f'<svg class="locus-fig"{core_attr} viewBox="0 0 {WIDTH} {height}" role="img" '
        f'aria-labelledby="{fig_id}-t {fig_id}-d" preserveAspectRatio="xMinYMin meet">'
        f'<title id="{fig_id}-t">Gene order of the {escape(label)} on {escape(contig)}</title>'
        f'<desc id="{fig_id}-d">{escape(desc)}</desc>'
    ]
    # Axis with ticks.
    step = _nice_step(span)
    parts.append(f'<line class="axis" x1="{x(lo):.1f}" y1="{AXIS_Y}" x2="{x(hi):.1f}" y2="{AXIS_Y}" '
                 f'stroke="{LINE_2}" stroke-width="1.5"/>')
    t = (lo // step + 1) * step
    while t < hi:
        tx = x(t)
        parts.append(f'<line class="tick" x1="{tx:.1f}" y1="{AXIS_Y - 4}" x2="{tx:.1f}" y2="{AXIS_Y + 4}" stroke="{LINE_2}"/>')
        parts.append(f'<text class="tick-label" x="{tx:.1f}" y="{AXIS_Y - 8}" text-anchor="middle" '
                     f'fill="{MUTED}" font-size="11" {FONT}>{_fmt_kb(t, step)}</text>')
        t += step
    # Contig ends: a cap when inside the window, else an arrow with the distance.
    for side, pos in ends:
        if lo <= pos <= hi:
            cx = x(pos)
            parts.append(f'<line class="contig-end" x1="{cx:.1f}" y1="{AXIS_Y - 9}" x2="{cx:.1f}" y2="{AXIS_Y + 9}" '
                         f'stroke="{INK}" stroke-width="2.5"/>')
            anchor, dx = ("start", 5) if side == "left" else ("end", -5)
            parts.append(f'<text class="end-label" x="{cx + dx:.1f}" y="{AXIS_Y + 14}" text-anchor="{anchor}" '
                         f'fill="{INK_2}" font-size="11" {FONT}>contig end</text>')
        else:
            dist = (lo - pos) if side == "left" else (pos - hi)
            dist_txt = f"{dist / 1000:,.1f} kb" if dist >= 1000 else f"{dist:,} bp"
            if side == "left":
                parts.append(f'<text class="end-label" x="{x(lo):.1f}" y="{AXIS_Y + 14}" text-anchor="start" '
                             f'fill="{MUTED}" font-size="11" {FONT}>← {dist_txt} more to contig end</text>')
            else:
                parts.append(f'<text class="end-label" x="{x(hi):.1f}" y="{AXIS_Y + 14}" text-anchor="end" '
                             f'fill="{MUTED}" font-size="11" {FONT}>{dist_txt} more to contig end →</text>')
    # Reported locus: a neutral bracket.
    for s in segs or [{"start": call["start"], "end": call["end"]}]:
        bx0, bx1, by = x(s["start"]), x(s["end"]), AXIS_Y + 26
        parts.append(f'<path class="bracket" d="M{bx0:.1f},{by - 4} V{by} H{bx1:.1f} V{by - 4}" fill="none" '
                     f'stroke="{INK_2}" stroke-width="1"/>')
        parts.append(f'<text class="bracket-label" x="{(bx0 + bx1) / 2:.1f}" y="{by - 3}" text-anchor="middle" '
                     f'fill="{INK_2}" font-size="10" {FONT}>called locus</text>')

    for g, lane in zip(genes, lanes):
        rc, sc = role_class(g.get("role")), status_class(g)
        h = CORE_H if rc == "r-core" else GENE_H
        mid = TRACK_TOP + lane * (LANE_H + LABEL_H) + LABEL_H + CORE_H / 2
        top = mid - h / 2
        x0, x1 = x(g["start"]), x(g["end"])
        minus = g.get("strand") == "-"
        tip = min(8.0, max(3.0, (x1 - x0) * 0.3))
        fill, stroke = ROLE_FILL[rc], ROLE_STROKE[rc]
        if sc == "s-hit":
            paint = f'fill="{stroke}" fill-opacity="0.3" stroke="{stroke}" stroke-width="1.5" stroke-dasharray="3 2"'
        elif sc == "s-disagree":
            paint = f'fill="{fill}" stroke="{INK}" stroke-width="1.5" stroke-dasharray="4 2"'
        else:
            paint = f'fill="{fill}" stroke="{stroke}" stroke-width="{2 if rc == "r-var" else 1.2}"'
        tooltip = (f"{g['gene']}: {_ROLE_WORDS.get(g.get('role'), 'other gene')}, {g['start']:,}–{g['end']:,} "
                   f"({g.get('strand') or '?'}), identity {g.get('identity')}%, {_STATUS_WORDS.get(sc)}")
        parts.append(f'<g class="gene {rc} {sc}"><title>{escape(tooltip)}</title>')
        parts.append(f'<line class="intron" x1="{x0:.1f}" y1="{mid:.1f}" x2="{x1:.1f}" y2="{mid:.1f}" '
                     f'stroke="{INK_2}" stroke-width="1"/>')
        exons = _exons(g)
        for i, (es, ee) in enumerate(exons):
            ex0, ex1 = x(es), x(ee)
            is_tip = (i == 0 and minus) or (i == len(exons) - 1 and not minus)
            if is_tip and ex1 - ex0 > 1.5:
                tl = min(tip, ex1 - ex0)
                if minus:
                    pts = [(ex0, mid), (ex0 + tl, top), (ex1, top), (ex1, top + h), (ex0 + tl, top + h)]
                else:
                    pts = [(ex0, top), (ex1 - tl, top), (ex1, mid), (ex1 - tl, top + h), (ex0, top + h)]
                parts.append('<polygon class="exon" points="' + " ".join(f"{a:.1f},{b:.1f}" for a, b in pts)
                             + f'" {paint}/>')
            else:
                parts.append(f'<rect class="exon" x="{ex0:.1f}" y="{top:.1f}" width="{max(1.5, ex1 - ex0):.1f}" '
                             f'height="{h}" {paint}/>')
        lx, anchor = (x0 + x1) / 2, "middle"
        if lx < MARGIN_X + 30:
            lx, anchor = max(x0, 2.0), "start"
        elif lx > WIDTH - MARGIN_X - 30:
            lx, anchor = min(x1, WIDTH - 2.0), "end"
        style = (f'fill="{INK_2}" font-style="italic"' if sc == "s-hit"
                 else f'fill="{INK}" font-weight="{700 if rc == "r-core" else 600}"')
        parts.append(f'<text class="gene-label" x="{lx:.1f}" y="{top - 4:.1f}" text-anchor="{anchor}" '
                     f'font-size="12" {FONT} {style}>{escape(str(g["gene"]))}</text>')
        alt = g.get("alternate_model")
        if sc == "s-disagree" and alt and alt.get("contig") == contig:
            ay = top + h + 3
            for e in (alt.get("exons") or [{"start": alt["start"], "end": alt["end"]}]):
                parts.append(f'<rect class="alt" x="{x(e["start"]):.1f}" y="{ay:.1f}" '
                             f'width="{max(1.0, x(e["end"]) - x(e["start"])):.1f}" height="3" fill="{INK_2}" opacity="0.75"/>')
        parts.append("</g>")
    parts.append("</svg>")
    return "".join(parts)
