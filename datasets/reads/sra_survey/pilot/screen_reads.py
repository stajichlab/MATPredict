#!/usr/bin/env python3
"""Screen a read set for the wrong organism or a mixture, from a GmG result. Standard library only.

The raw-reads typing finds the MAT k-mers of the expected species in whatever the reads hold, so a read set that is another
species, or a mix of two, is still typed. GmG (https://github.com/stajichlab/GmG) classifies the same reads against a fungal
reference and reports a mixture block. This turns its record into flags that travel with the typing result; it never
stops a run from being typed.

  screen_reads.py trigger  --results results.jsonl --run RUN              prints "deeper" or "ok"
  screen_reads.py finalize --results results.jsonl --run RUN --group SLUG --spots N [--deep-results results.jsonl --deep-spots M]
                           --out RUN.screen.json

The rule (calibrated in results/2026-11-02_contam_spikein/NOTE.md):
  mixed             GmG's own `mixture_suspected` flag, at 2M spots or at the deeper screen. A contaminant is flagged at about
                    2% of the pairs at 2M spots; 0.5% needs about 6M.
  genus_mismatch    the genus of the species the run was selected for is neither GmG's genus call nor in the lineage of its top
                    hits: the whole read set is another organism. It can also be a renamed genus, so it flags and does not exclude.
  borderline        not flagged, but a second genus with at least 5 hashes and a concentration of at least 0.7 (most of its
                    hashes in one reference genome; chance second genera have a concentration near 0.2). That triggers a deeper
                    screen; a set still not flagged after it, or that cannot be read deeper, keeps `borderline`.
  no_call           GmG made no genus call.
  ok                none of the above.
`status` is the first that applies, in the order genus_mismatch, mixed, borderline, no_call, ok; `flags` lists all that apply.
"""
import argparse
import json
import sys

MIN_HASHES, MIN_CONCENTRATION = 5, 0.7
PRECEDENCE = ("genus_mismatch", "mixed", "borderline", "no_call")


def expected_genus(group):
    """The genus of the species a run was selected for, from its list label ("aspergillus_fumigatus" gives Aspergillus), or ""."""
    first = (group or "").replace(" ", "_").split("_")[0]
    return first.capitalize()


def load_record(path, run):
    """The GmG record (a dict) for `run` in a results.jsonl, or None."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                rec = json.loads(line)
                if rec.get("id") == run:
                    return rec
    return None


def facts(rec):
    """The fields of one GmG record that the rule uses; every missing piece is None."""
    if not rec:
        return {"present": False}
    qc = rec.get("qc") or {}
    mix = qc.get("mixture") or {}
    second = mix.get("second") or {}
    calls = (rec.get("guess") or {}).get("calls") or {}
    hits = (rec.get("guess") or {}).get("hits") or []
    species = calls.get("species") or {}
    return {
        "present": True, "status": rec.get("status"), "gmg_flags": list(rec.get("qc_flags") or []),
        "genus": (calls.get("genus") or {}).get("taxon") or None, "species": species.get("taxon") or species.get("vote") or None,
        "hit_lineages": [h.get("lineage") or [] for h in hits[:5]],
        "mixed": bool(mix.get("suspected")), "primary": mix.get("primary"), "explained": mix.get("explained_fraction"),
        "second": second.get("taxon") or None, "second_hashes": second.get("hashes"), "second_explained": second.get("explained_fraction"),
        "second_concentration": second.get("concentration"), "coverage_est": qc.get("coverage_est"),
    }


def is_borderline(f):
    """A second genus that looks like a real organism and not a chance match, in a set GmG did not flag."""
    return bool(f.get("present") and not f["mixed"] and f["second"] and (f["second_hashes"] or 0) >= MIN_HASHES
                and (f["second_concentration"] or 0) >= MIN_CONCENTRATION)


def trigger(rec):
    """"deeper" when a deeper screen should be run, else "ok"."""
    return "deeper" if is_borderline(facts(rec)) else "ok"


def genus_found(f, genus):
    """True when `genus` is GmG's genus call or appears among the genera of its top hits."""
    if not genus:
        return True
    return genus == f.get("genus") or any(genus in lin for lin in f.get("hit_lineages", []))


def finalize(rec, group, spots, deep_rec=None, deep_spots=None, deep_possible=True):
    """The screen result of one run: status, flags and the numbers behind them. `rec` is the GmG record at `spots`
    spots; `deep_rec` the record of the deeper screen when one was run; `deep_possible` is False when the run has no
    more reads to give."""
    f = facts(rec)
    exp = expected_genus(group)
    flags = []
    out = {"expected_genus": exp or None, "screened_spots": spots, "rescreened": deep_rec is not None, "rescreen_spots": deep_spots if deep_rec else None}
    if not f["present"]:
        out.update({"status": "no_result", "flags": ["no_result"]})
        return out
    deeper = facts(deep_rec) if deep_rec else None
    used = deeper if (deeper and deeper["present"]) else f
    if not genus_found(f, exp):
        flags.append("genus_mismatch")
    if f["mixed"] or (deeper and deeper["present"] and deeper["mixed"]):
        flags.append("mixed")
    elif is_borderline(f) and not (deeper and deeper["present"] and not deeper["mixed"] and not is_borderline(deeper)):
        flags.append("borderline")
    if not f["genus"] and not (deeper and deeper["genus"]):
        flags.append("no_call")
    status = next((p for p in PRECEDENCE if p in flags), "ok")
    out.update({
        "status": status, "flags": flags, "genus": used["genus"], "species": used["species"], "gmg_flags": used["gmg_flags"],
        "second_genus": used["second"], "second_hashes": used["second_hashes"], "second_explained": used["second_explained"],
        "second_concentration": used["second_concentration"], "coverage_est": used["coverage_est"],
        "deeper_unavailable": bool(is_borderline(f) and deep_rec is None and not deep_possible) or None,
    })
    if deeper and deeper["present"]:
        out["at_screened_spots"] = {k: f[k] for k in ("second", "second_hashes", "second_explained", "second_concentration", "mixed")}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("trigger")
    t.add_argument("--results", required=True)
    t.add_argument("--run", required=True)
    z = sub.add_parser("finalize")
    z.add_argument("--results", required=True)
    z.add_argument("--run", required=True)
    z.add_argument("--group", default="")
    z.add_argument("--spots", type=int, required=True)
    z.add_argument("--deep-results")
    z.add_argument("--deep-spots", type=int)
    z.add_argument("--no-deeper", action="store_true", help="the run has no more reads than the first screen used")
    z.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "trigger":
        print(trigger(load_record(a.results, a.run)))
        return 0
    deep = load_record(a.deep_results, a.run) if a.deep_results else None
    res = finalize(load_record(a.results, a.run), a.group, a.spots, deep, a.deep_spots, deep_possible=not a.no_deeper)
    with open(a.out, "w") as fh:
        json.dump(res, fh)
    print(json.dumps({k: res.get(k) for k in ("status", "flags", "genus", "second_genus", "second_hashes", "second_explained")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
