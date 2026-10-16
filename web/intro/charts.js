// Tiny SVG chart helpers (no dependencies). Bars start at zero.
export const COLORS = ["var(--c1)","var(--c2)","var(--c3)","var(--c4)","var(--c5)","var(--c6)","var(--c7)","var(--c8)"];
export const fmt = n => n.toLocaleString("en-US");
export const pct = (a, b, d = 1) => b ? (100 * a / b).toFixed(d) + "%" : "n/a";
const NEUTRAL = { "other allele": "var(--c8)", "other single allele": "var(--c8)", "two alleles": "#8c8c91", "undetermined": "var(--c5)", "uncalled": "var(--line)" };
export const colorFor = (c, j) => NEUTRAL[c] || COLORS[j % 6];
const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

export async function loadData() {
  const [taxa, summary] = await Promise.all([fetch("data/taxa.json").then(r => r.json()), fetch("data/summary.json").then(r => r.json())]);
  const S = taxa.strings, F = taxa.fields;
  const rows = taxa.rows.map(r => {
    const o = {};
    F.forEach((f, i) => { if (i < 8) o[f] = S[r[i]]; else if (i < 12) o[f] = r[i]; });
    o.alleles = Object.fromEntries(Object.entries(r[12]).map(([k, v]) => [S[k], v]));
    o.types = Object.fromEntries(Object.entries(r[13]).map(([k, v]) => [S[k], v]));
    return o;
  });
  return { rows, summary };
}

// Horizontal bars. items: [{label, value, max, note}]; value shown as fraction of max (default 1).
export function hbar(el, items, { fmtv = v => v, max = null, labelw = 190 } = {}) {
  const W = 760, rh = 26, m = max ?? Math.max(...items.map(i => i.value), 1e-9), bw = W - labelw - 150;
  let s = `<svg viewBox="0 0 ${W} ${items.length * rh + 4}" role="img">`;
  items.forEach((it, i) => {
    const y = i * rh + 2, w = Math.max(0, bw * it.value / m);
    s += `<text x="${labelw - 8}" y="${y + 16}" text-anchor="end">${esc(it.label)}</text><rect x="${labelw}" y="${y + 3}" width="${bw}" height="16" rx="3" fill="var(--card)"/>` +
      `<rect x="${labelw}" y="${y + 3}" width="${w}" height="16" rx="3" fill="var(--c1)"><title>${esc(it.label)}: ${esc(it.note || fmtv(it.value))}</title></rect>` +
      `<text x="${labelw + bw + 8}" y="${y + 16}" class="m">${esc(it.note || fmtv(it.value))}</text>`;
  });
  el.innerHTML = s + "</svg>";
}

// 100% stacked bars. items: [{label, parts: {cat: n}, total?}], cats: ordered list of category names.
export function stacked(el, items, cats, { labelw = 190, rh = 30 } = {}) {
  const W = 760, bw = W - labelw - 70;
  let s = `<svg viewBox="0 0 ${W} ${items.length * rh + 4}" role="img">`;
  items.forEach((it, i) => {
    const tot = it.total ?? Object.values(it.parts).reduce((a, b) => a + b, 0), y = i * rh + 2;
    s += `<text x="${labelw - 8}" y="${y + 17}" text-anchor="end">${esc(it.label)}</text>`;
    let x = labelw;
    cats.forEach((c, j) => {
      const v = it.parts[c] || 0; if (!v || !tot) return;
      const w = bw * v / tot;
      s += `<rect x="${x}" y="${y + 3}" width="${w}" height="20" fill="${colorFor(c, j)}"><title>${esc(it.label)} — ${esc(c)}: ${fmt(v)} (${pct(v, tot)})</title></rect>`;
      if (w > 34) s += `<text x="${x + w / 2}" y="${y + 17}" text-anchor="middle" style="fill:#fff;font-size:11px">${(100 * v / tot).toFixed(0)}%</text>`;
      x += w;
    });
    s += `<text x="${labelw + bw + 8}" y="${y + 17}" class="m">${fmt(tot)}</text>`;
  });
  el.innerHTML = s + "</svg>";
}

export function legend(el, cats) {
  el.innerHTML = cats.map((c, j) => `<span><i style="background:${colorFor(c, j)}"></i>${esc(c)}</span>`).join("");
}
