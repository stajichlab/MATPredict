// What to do with a GmG record: which of the page's phylum options to select, or why not. The same rules as
// src/MATPredict/detect/gmg_guess.py (`guess_from_record`, `decide`), so the browser and `detect --phylum auto` agree;
// gmg_route.test.mjs holds the shared cases. A wrong phylum is worse than none (the right family is never searched), so a
// call GmG marks weak, or an input it says is not a fungus, is not used unless the visitor accepts weak calls.

// The GmG record's routing_hint (GmG schema 9), or the same rule applied to an older record.
export function guessFromRecord(rec) {
  const calls = (rec.guess && rec.guess.calls) || {};
  const flags = rec.qc_flags || [];
  const call = (rank) => calls[rank] || {};
  let hint = rec.routing_hint;
  if (!hint) {
    const phylum = call("phylum").taxon || null;
    const reason = rec.status === "error" ? "error"
      : !phylum ? "no_phylum_call"
      : flags.includes("not_fungal") ? "not_fungal"
      : flags.includes("weak_phylum") ? "weak_phylum" : null;
    hint = { usable: reason === null, reason, phylum };
  }
  return { status: rec.status || null, usable: !!hint.usable, reason: hint.reason || null, phylum: hint.phylum || null, flags };
}

// -> { phylum: string|null, weak: bool, note: string }. `available` is the phylum names the page offers.
export function decide(guess, available, acceptWeak = false) {
  const ph = guess.phylum;
  if (guess.reason === "error") return { phylum: null, weak: false, note: "GmG could not read this genome." };
  if (!ph || guess.reason === "no_phylum_call") return { phylum: null, weak: false, note: "GmG made no phylum call for this genome." };
  if (guess.reason === "not_fungal") return { phylum: null, weak: false, note: `GmG says this is probably not a fungus (nearest phylum ${ph}).` };
  const weak = guess.reason === "weak_phylum";
  if (weak && !acceptWeak) {
    return { phylum: null, weak: true, note: `GmG's guess is ${ph}, but it is a weak call (best-hit ANI below 0.85). Tick “Accept weak guesses” to use it.` };
  }
  if (guess.reason && !weak) return { phylum: null, weak: false, note: `GmG marks its call ${ph} as not usable (${guess.reason}).` };
  if (!available.includes(ph)) return { phylum: null, weak: false, note: `GmG calls phylum ${ph}, which has no curated MAT locus family here.` };
  return { phylum: ph, weak, note: `GmG calls phylum ${ph}${weak ? " (a weak call, accepted)" : ""}.` };
}
