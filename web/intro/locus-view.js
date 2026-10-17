// Locus maps for the explorer: one card per called genome, one small gene map per locus, contig accession linked to NCBI.
const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const ACC = /^(?:[A-Z]{1,2}_)?[A-Z]{1,6}\d{5,}(?:\.\d+)?$/;   // GenBank / RefSeq nucleotide accession; lab scaffold names do not match
const GCA = /^(GC[AF]_\d+\.\d+)/;
const fmt = n => n.toLocaleString("en-US");

export async function loadShard(path) {
  const buf = await (await fetch(path)).arrayBuffer();
  const text = await new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream("gzip"))).text();
  return JSON.parse(text).loci;
}

export const contigLink = (c, s, e) => ACC.test(c)
  ? `<a href="https://www.ncbi.nlm.nih.gov/nuccore/${encodeURIComponent(c)}?from=${s}&to=${e}" target="_blank" rel="noopener">${esc(c)}</a>` : `<span>${esc(c)}</span>`;

const ROLE = { core_MAT: "var(--c1)", flanking_conserved: "var(--c3)", flanking_variable: "var(--c5)" };
const SOLID = new Set(["polished_agree", "polished_single", "polished_disagree"]);

// SVG gene map. Genes are placed in lanes so neither arrows nor labels overlap.
export function locusSvg(locus) {
  const [contig, start, end, , , , genes] = locus;
  const W = 720, pad = 14, top = 8;
  const lo0 = Math.min(start, ...genes.map(g => g[1])), hi0 = Math.max(end, ...genes.map(g => g[2]));
  const margin = Math.max(200, (hi0 - lo0) * 0.04), lo = lo0 - margin, hi = hi0 + margin, span = hi - lo;
  const x = p => pad + (W - 2 * pad) * (p - lo) / span;
  const lanes = [];
  const placed = [...genes].sort((a, b) => a[1] - b[1]).map(g => {
    const x0 = x(g[1]), x1 = x(g[2]), labelEnd = x0 + g[0].length * 6.4 + 10, need = Math.max(x1, labelEnd);
    let lane = lanes.findIndex(end => end < x0); if (lane < 0) { lane = lanes.length; lanes.push(0); }
    lanes[lane] = need + 4; return { g, x0, x1, lane };
  });
  const laneH = 30, axisY = top + Math.max(1, lanes.length) * laneH + 6, H = axisY + 22;
  let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Gene map of ${esc(contig)}:${start}-${end}">`;
  s += `<rect x="${x(start)}" y="${top - 4}" width="${Math.max(2, x(end) - x(start))}" height="${axisY - top + 4}" fill="var(--locus)" opacity=".85"><title>Called locus ${fmt(start)}–${fmt(end)}</title></rect>`;
  for (const { g, x0, x1, lane } of placed) {
    const [name, gs, ge, strand, role, status] = g, y = top + lane * laneH + 13, h = 12, w = Math.max(3, x1 - x0), hd = Math.min(7, w);
    const fwd = strand !== "-", pts = fwd ? `${x0},${y} ${x0 + w - hd},${y} ${x0 + w},${y + h / 2} ${x0 + w - hd},${y + h} ${x0},${y + h}`
                                         : `${x0 + hd},${y} ${x0 + w},${y} ${x0 + w},${y + h} ${x0 + hd},${y + h} ${x0},${y + h / 2}`;
    const col = ROLE[role] || "var(--other)", solid = SOLID.has(status);
    s += `<polygon points="${pts}" fill="${solid ? col : "none"}" stroke="${col}" stroke-width="${solid ? 0 : 1.5}" ${solid ? "" : 'stroke-dasharray="3 2"'}><title>${esc(name)} ${fmt(gs)}–${fmt(ge)} (${esc(strand)}), ${esc(role)}, ${esc(status)}</title></polygon>`;
    const right = x0 + name.length * 6.4 > W - pad;
    s += `<text x="${right ? Math.max(x1, x0) : x0}" y="${y - 3}" font-size="11"${right ? ' text-anchor="end"' : ""}>${esc(name)}</text>`;
  }
  s += `<line x1="${pad}" x2="${W - pad}" y1="${axisY}" y2="${axisY}" stroke="var(--muted)"/>`;
  s += `<text x="${pad}" y="${axisY + 14}" font-size="11" class="m">${fmt(Math.round(lo))}</text><text x="${W - pad}" y="${axisY + 14}" font-size="11" text-anchor="end" class="m">${fmt(Math.round(hi))}</text>`;
  const tick = [100, 500, 1000, 2000, 5000, 10000, 20000, 50000].find(t => (W - 2 * pad) * t / span > 40 && (W - 2 * pad) * t / span < 220) || 1000;
  s += `<line x1="${W / 2 - (W - 2 * pad) * tick / span / 2}" x2="${W / 2 + (W - 2 * pad) * tick / span / 2}" y1="${axisY + 10}" y2="${axisY + 10}" stroke="var(--fg)" stroke-width="2"/><text x="${W / 2}" y="${axisY + 20}" font-size="10" text-anchor="middle" class="m">${tick >= 1000 ? tick / 1000 + " kb" : tick + " bp"}</text>`;
  return s + "</svg>";
}

export const LEGEND = `<div class="legend"><span><i style="background:var(--c1)"></i>core MAT gene</span><span><i style="background:var(--c3)"></i>conserved flank gene</span>` +
  `<span><i style="background:var(--c5)"></i>variable flank gene</span><span><i style="border:1.5px dashed var(--muted)"></i>unpolished (search hit only)</span><span><i style="background:var(--locus)"></i>called locus span</span></div>`;

// genomes: [{dataset, genome, species, status, call}]; shards: {dataset: {genome: [locus...]}}
export function renderCards(el, genomes, shards, { limit = 40 } = {}) {
  const called = genomes.filter(g => (shards[g.dataset] || {})[g.genome]);
  const uncalled = genomes.filter(g => !(shards[g.dataset] || {})[g.genome]);
  called.sort((a, b) => ((a.species_recorded || a.species) + a.genome).localeCompare((b.species_recorded || b.species) + b.genome));
  const card = g => {
    const asm = GCA.exec(g.genome), loci = shards[g.dataset][g.genome];
    return `<div class="card locus-card"><div><b><i>${esc(g.species_recorded || g.species || "unnamed")}</i></b> · ` +
      (asm ? `<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/${asm[1]}/" target="_blank" rel="noopener">${esc(g.genome)}</a>` : esc(g.genome)) +
      ` <span class="tag">${esc(g.call)}</span></div>` +
      loci.map(l => `<div class="muted" style="margin-top:8px">${contigLink(l[0], l[1], l[2])}:${fmt(l[1])}–${fmt(l[2])} · ${esc(l[3])} · ${esc(l[4])} · ${esc(l[5].replace(/_/g, " "))}</div>${locusSvg(l)}`).join("") + `</div>`;
  };
  let shown = 0;
  const more = () => { const part = called.slice(shown, shown + limit); shown += part.length; el.querySelector(".cards-list").insertAdjacentHTML("beforeend", part.map(card).join("")); el.querySelector(".more").hidden = shown >= called.length; };
  el.innerHTML = `<p><b>${fmt(called.length)}</b> of ${fmt(genomes.length)} genomes in this view have a call.</p>${called.length ? LEGEND : ""}<div class="cards-list"></div>` +
    `<p><button class="more">Show more</button></p>` +
    (uncalled.length ? `<details><summary>${fmt(uncalled.length)} without a call</summary><ul>${uncalled.map(g => `<li><i>${esc(g.species_recorded || g.species || "unnamed")}</i> · ${esc(g.genome)} <span class="muted">(${esc(g.call)})</span></li>`).join("")}</ul></details>` : "");
  el.querySelector(".more").onclick = more; more();
}
