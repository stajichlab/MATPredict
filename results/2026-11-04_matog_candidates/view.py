"""Local view of the candidate MATog records and the leave-one-out table: python3 view.py -> view.html"""
import csv, glob, html, yaml
COL = {"core_MAT": "#c0392b", "flanking_conserved": "#2e6da4"}
def svg(rec):
    out = []
    for seg in rec["locus"]["core"]["segments"]:
        gs = [g for g in rec["genes"] if g["segment_index"] == seg["segment_index"]]
        lo, hi = seg["start"], seg["end"]; W = 900; sc = W / max(1, hi - lo)
        s = f'<div class="seg">{seg["sequence_source"]["accession"]}:{lo:,}-{hi:,} ({hi-lo+1:,} bp)</div><svg width="{W+20}" height="64">'
        s += f'<line x1="10" x2="{W+10}" y1="32" y2="32" stroke="#999"/>'
        for g in gs:
            x1 = 10 + (g["start"] - lo) * sc; x2 = 10 + (g["end"] - lo) * sc; up = g["strand"] == "+"
            pts = f"{x1},20 {x2-8 if up else x2},20 {x2 if up else x2},32 {x2-8 if up else x2},44 {x1},44" if up else f"{x2},20 {x1+8},20 {x1},32 {x1+8},44 {x2},44"
            s += f'<polygon points="{pts}" fill="{COL[g["role"]]}" opacity=".85"/><text x="{x1}" y="14" font-size="11">{g["name"]} {g["end"]-g["start"]+1}bp</text>'
        out.append(s + "</svg>")
    return "".join(out)
body = ""
for f in sorted(glob.glob("records/*.yaml")):
    r = yaml.safe_load(open(f)); hl = " flag" if r["record_id"].startswith("51923") else ""
    body += f'<section class="{hl}"><h3>{r["record_id"]}</h3><p>{r["organism"]["species"]} {r["organism"]["strain"]["name"]} - {r["locus"]["assembly_accession"]}</p>{svg(r)}</section>'
rows = list(csv.DictReader(open("compare.tsv"), delimiter="\t"))
tab = "<table><tr>" + "".join(f"<th>{h}</th>" for h in ("species", "base cassettes", "change", "base calls", "new calls")) + "</tr>"
for r in rows:
    if r["change"] in ("calls_differ", "gained_idiomorph", "LOST_idiomorph"):
        fl = ' class="flag"' if "MSSE01000034" in r["new_calls"] else ""
        tab += f'<tr{fl}><td>{r["species"]}</td><td>{r["base_cassettes"]}</td><td>{r["change"]}</td><td>{html.escape(r["base_calls"])}</td><td>{html.escape(r["new_calls"])}</td></tr>'
open("view.html", "w").write(f"""<!doctype html><meta charset=utf-8><title>MATog candidates</title><style>body{{font:14px system-ui;margin:2em;max-width:1100px}}section{{margin:1em 0;padding:.5em 1em;border:1px solid #ccc}}.flag{{background:#fff4d6}}table{{border-collapse:collapse;font-size:12px}}td,th{{border:1px solid #ccc;padding:2px 6px;vertical-align:top}}.seg{{font-size:12px;color:#555}}</style>
<h1>Candidate MATog records (red = MAT gene, blue = flank)</h1>{body}<h2>Leave-one-out: genomes whose calls changed (yellow = extra boidinii call)</h2>{tab}</table>""")
