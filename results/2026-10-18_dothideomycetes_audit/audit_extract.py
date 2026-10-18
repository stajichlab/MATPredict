"""Per-locus audit table for the Dothideomycetes run: report fields, core-gene evidence, flank genes, alpha-box/HMG-box
profile scores (own-genus-out profile where a reference of the genome's genus exists; order-out as a second view) and a
Pfam domain check (GA cutoffs) of the core and optional MAT proteins.

    audit_extract.py --runs RUNS --genomes genome_list.tsv --lib LIB --db DB --pfam-subset pfam_subset.hmm --out OUTDIR
Reads only; changes no code or database.  Run with PYTHONPATH=<worktree>/src and the pixi `test` env on PATH.
"""
import argparse, csv, gzip, json, os, sys, tempfile
from collections import defaultdict
from pathlib import Path

import yaml
from Bio import SeqIO
from Bio.Seq import Seq

L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
CORE = ("MAT1-1-1", "MAT1-2-1")
OWN_CORE = {"MAT1-1": "MAT1-1-1", "MAT1-2": "MAT1-2-1"}
MARGIN = 25.0  # provisional profile margin used by scripts/score_core_profiles.py
MODELLED = ("polished_agree", "polished_single", "polished_disagree")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def optional_and_idiomorph(db: Path):
    d = yaml.safe_load((db / "Ascomycota" / "order.yml").read_text())
    locus = next(l for l in d["loci"] if l["locus_name"] == "MAT")
    opt, idio, role = {}, {}, {}
    for g in locus["genes"]:
        opt[g["name"]] = bool(g.get("optional"))
        idio[g["name"]] = (g.get("present_in_idiomorphs") or [None])[0]
        role[g["name"]] = g.get("role")
    return opt, idio, role


def cds(entry, seq):
    if entry.get("exons"):
        ex = sorted((int(x["start"]), int(x["end"])) for x in entry["exons"])
    else:
        ex = [(int(entry["start"]), int(entry["end"]))]
    dna = "".join(seq[a - 1:b] for a, b in ex)
    return dna.translate(COMP)[::-1] if entry["strand"] == "-" else dna


def proteins(entry, seq):
    dna = cds(entry, seq)
    if entry.get("exons"):
        return [str(Seq(dna[: len(dna) // 3 * 3]).translate()).replace("*", "")]
    return [str(Seq(dna[f:][: (len(dna) - f) // 3 * 3]).translate()).replace("*", "X") for f in range(3)]


def contigs(lib, genome, wanted):
    """Only the wanted contigs, read with a plain line parser (Biopython over every contig was the slow step)."""
    p = Path(lib) / (genome + ".fa.gz")
    out, cur, buf = {}, None, []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf).upper()
                name = line[1:].split()[0]
                cur = name if name in wanted else None
                buf = []
                if len(out) == len(wanted) and cur is None and out:
                    pass
            elif cur is not None:
                buf.append(line.strip())
    if cur is not None:
        out[cur] = "".join(buf).upper()
    return out


def build_profiles(db, drop_ids, work):
    from MATPredict.detect import classifier_build as cb
    from MATPredict.detect.family_registry import FamilyKey, load_all_families
    family = next(f for f in load_all_families(db) if f.key == FamilyKey("Ascomycota", "MAT"))
    rows = [r for r in cb.training_set(db, family, work) if r["gene"] in CORE
            and r["id"].split("|")[1] not in drop_ids]
    out, n = {}, {}
    for gene in CORE:
        seqs = {r["id"]: r["sequence"] for r in rows if r["gene"] == gene}
        n[gene] = len(seqs)
        if not seqs:
            return None, n
        out[gene], _ = cb.build_hmm(seqs, gene, work)
    return out, n


def hmm_score(hmm, prots):
    from MATPredict.detect.classifier_build import score
    sc = score(hmm, {f"p{i}": p for i, p in enumerate(prots) if p})
    return max(sc.values(), default=0.0)


def prof_scores(hmms, prots):
    return {g: hmm_score(h, prots) for g, h in hmms.items()}


def pcall(a, h):
    if a - h >= MARGIN:
        return "MAT1-1"
    if h - a >= MARGIN:
        return "MAT1-2"
    return "undetermined"


def pfam_hits(pfam, prots):
    """best GA-passing Pfam hit per accession name over the protein candidates."""
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextSequence
    alpha = Alphabet.amino()
    seqs = []
    for i, p in enumerate(prots):
        if len(p) >= 20:
            seqs.append(TextSequence(name=f"p{i}".encode(), sequence=p).digitize(alpha))
    best = {}
    if not seqs:
        return best
    for hits in pyhmmer.hmmsearch(pfam, seqs, bit_cutoffs="gathering", cpus=1):
        for h in hits:
            if h.included:
                n = hits.query.name.decode() if isinstance(hits.query.name, bytes) else hits.query.name
                best[n] = max(best.get(n, 0.0), h.score)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True); ap.add_argument("--genomes", required=True)
    ap.add_argument("--lib", required=True); ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--pfam-subset", required=True); ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--genome-list", default=None, help="restrict to these ASM ids (one per line)")
    a = ap.parse_args()
    import pyhmmer
    a.out.mkdir(parents=True, exist_ok=True)
    meta = {r["ASM_ACCESSION"]: r for r in csv.DictReader(open(a.genomes), delimiter="\t")}
    opt, idio_of, role = optional_and_idiomorph(a.db)
    pfam = list(pyhmmer.plan7.HMMFile(a.pfam_subset))
    # profiles: full, genus-out for the genera that have a reference, order-out for the two orders that do
    from MATPredict.detect import classifier_build as cb  # noqa: F401
    rec_order = {p.parent.name: p.parent.parent.name for p in (a.db / "Ascomycota").glob("*/*/metadata.yaml")}
    rec_genus = {}
    for p in (a.db / "Ascomycota").glob("*/*/metadata.yaml"):
        rec_genus[p.parent.name] = yaml.safe_load(p.read_text())["organism"]["species"].split()[0]
    work = Path(tempfile.mkdtemp())
    profiles = {"full": build_profiles(a.db, set(), work)[0]}
    needed = {m.get("GENUS") for m in meta.values()}
    for g in sorted(set(rec_genus.values()) & needed):
        profiles["genus:" + g] = build_profiles(a.db, {r for r, x in rec_genus.items() if x == g}, work)[0]
    for o in ("Pleosporales", "Mycosphaerellales"):
        profiles["order:" + o] = build_profiles(a.db, {r for r, x in rec_order.items() if x == o}, work)[0]
    json.dump({k: v is not None for k, v in profiles.items()}, open(a.out / "profiles_built.json", "w"))

    rows, grows = [], []
    import time as _t
    _t0 = _t.time()
    for _i, run in enumerate(sorted(Path(a.runs).iterdir())):
        asm = run.name
        if _i % 20 == 0:
            print(f"{_i} genomes {len(rows)} loci {_t.time() - _t0:.0f}s", file=sys.stderr, flush=True)
        rep = run / "detection_report.yaml"
        m = meta.get(asm, {})
        g = dict(genome=asm, order=m.get("ORDER", ""), genus=m.get("GENUS", ""), species=m.get("SPECIES", ""),
                 taxid=m.get("NCBI_TAXID", ""), status="no_report", n_detected=0, idiomorphs="", wall=None)
        ws = run / "wall_seconds"
        if ws.exists():
            g["wall"] = ws.read_text().strip()
        if not rep.exists() or rep.stat().st_size == 0:
            err = (run / "stderr.log").read_text()[-300:].replace("\n", " ") if (run / "stderr.log").exists() else ""
            g["error"] = err
            grows.append(g); continue
        doc = yaml.load(open(rep), Loader=L)
        det = doc.get("detected") or []
        g["status"] = "called" if det else "not_detected"
        g["n_detected"] = len(det)
        g["routing"] = doc.get("routing_mode")
        g["suppressed_unpolished"] = doc.get("suppressed_unpolished")
        g["two_idiomorphs"] = json.dumps(doc.get("two_idiomorphs") or [])
        g["zygosity"] = json.dumps(doc.get("zygosity")) if doc.get("zygosity") else ""
        g["idiomorphs"] = "|".join(str(d.get("idiomorph")) for d in det)
        grows.append(g)
        if not det:
            continue
        wanted = {d["contig"] for d in det}
        for d in det:
            for ge in d.get("gene_evidence") or []:
                wanted.add(ge["contig"])
        try:
            seqs = contigs(a.lib, asm, wanted)
        except Exception as exc:  # noqa: BLE001
            seqs = {}
        for d in det:
            ev = [e for e in (d.get("gene_evidence") or [])]
            live = [e for e in ev if not e.get("lost_to")]
            by = defaultdict(list)
            for e in live:
                by[e["gene"]].append(e)
            best = {k: max(v, key=lambda e: e.get("bitscore") or 0) for k, v in by.items()}
            idio = d.get("idiomorph")
            cand = {c["idiomorph"]: c["score"] for c in d.get("idiomorph_candidates") or []}
            # winning gene: the called idiomorph's highest-bitscore live gene that is informative for it
            win, wscore = None, None
            for k, e in best.items():
                if idio_of.get(k) == idio and role.get(k) == "core_MAT":
                    if wscore is None or (e.get("bitscore") or 0) > wscore:
                        win, wscore = k, e.get("bitscore") or 0
            own = best.get(OWN_CORE.get(idio, ""), None)
            oth_idio = "MAT1-2" if idio == "MAT1-1" else "MAT1-1" if idio == "MAT1-2" else None
            oth = best.get(OWN_CORE.get(oth_idio, ""), None) if oth_idio else None

            def fx(e, pre):
                if not e:
                    return {pre + k: "" for k in ("bits", "id", "cov", "status", "edge")}
                s = seqs.get(e["contig"])
                edge = ""
                if s:
                    edge = min(int(e["start"]) - 1, len(s) - int(e["end"]))
                return {pre + "bits": e.get("bitscore"), pre + "id": e.get("identity"), pre + "cov": e.get("coverage"),
                        pre + "status": e.get("status"), pre + "edge": edge}

            row = dict(genome=asm, order=g["order"], genus=g["genus"], species=g["species"],
                       contig=d["contig"], start=d["start"], end=d["end"], family=d.get("family"),
                       idiomorph=idio, idiomorph_class=d.get("idiomorph_class"), confidence=d.get("confidence"),
                       locus_class=d.get("locus_class"), detection_pass=d.get("detection_pass"),
                       margin=d.get("idiomorph_margin"), cand11=cand.get("MAT1-1"), cand12=cand.get("MAT1-2"),
                       vote_note=json.dumps(d.get("idiomorph_vote_note")) if d.get("idiomorph_vote_note") else "",
                       resolutions=json.dumps(d.get("idiomorph_resolutions") or []),
                       ambiguous_with=json.dumps(d.get("ambiguous_with") or []),
                       genes_found="|".join(d.get("genes_found") or []), genes_missing="|".join(d.get("genes_missing") or []),
                       n_loci_genome=len(det), fragmented=d.get("fragmented"),
                       win_gene=win, win_bits=wscore, win_optional=opt.get(win) if win else None,
                       optional_found="|".join(sorted(k for k in best if opt.get(k) and role.get(k) == "core_MAT")),
                       flank_found="|".join(sorted(k for k in best if role.get(k) != "core_MAT")))
            pairs = [(r.get("winner"), r.get("loser")) for r in (d.get("idiomorph_resolutions") or [])]
            row["res_pairs"] = "|".join(f"{w}>{l}" for w, l in pairs)
            # a core gene set aside in favour of an optional gene (the Parastagonospora v0.6.0 mechanism)
            row["core_lost_to_optional"] = any(l in CORE and opt.get(w) for w, l in pairs)
            row["core_won_over_optional"] = any(w in CORE and opt.get(l) for w, l in pairs)
            lost = [e for e in ev if e.get("lost_to")]
            row["lost_genes"] = "|".join(f"{e['gene']}:{e.get('identity')}:{e.get('coverage')}>{e['lost_to']}" for e in lost)
            row.update(fx(best.get("MAT1-1-1"), "a_")); row.update(fx(best.get("MAT1-2-1"), "h_"))
            row.update(fx(own, "own_")); row.update(fx(oth, "oth_"))
            row["own_full"] = bool(own and (own.get("coverage") or 0) >= 90 and own.get("status") in MODELLED)
            row["own_present"] = own is not None
            row["both_core_in_locus"] = ("MAT1-1-1" in best and "MAT1-2-1" in best)
            # profile + pfam
            if seqs:
                prot = {}
                for k in ("MAT1-1-1", "MAT1-2-1"):
                    e = best.get(k)
                    if e and e["contig"] in seqs:
                        prot[k] = proteins(e, seqs[e["contig"]])
                if win and win not in prot and win in best and best[win]["contig"] in seqs:
                    prot[win] = proteins(best[win], seqs[best[win]["contig"]])
                genus = g["genus"]; order = g["order"]
                pkey = "genus:" + genus if ("genus:" + genus) in profiles else "full"
                for tag, key in (("prof", pkey), ("prof_order", "order:" + order if ("order:" + order) in profiles else None)):
                    hm = profiles.get(key) if key else None
                    if hm is None:
                        row[tag + "_set"] = key or ""
                        continue
                    sc = {k: prof_scores(hm, p) for k, p in prot.items()}
                    aa = sc.get("MAT1-1-1", {}).get("MAT1-1-1"); hh = sc.get("MAT1-2-1", {}).get("MAT1-2-1")
                    if aa is not None and hh is not None:
                        pc, mg = pcall(aa, hh), round(aa - hh, 1)
                    elif sc:
                        o = next(iter(sc.values())); mg = round(o["MAT1-1-1"] - o["MAT1-2-1"], 1)
                        pc = pcall(o["MAT1-1-1"], o["MAT1-2-1"])
                    else:
                        pc, mg = "no_protein", None
                    row[tag + "_set"] = key
                    row[tag + "_alpha_on_alpha"] = None if aa is None else round(aa, 1)
                    row[tag + "_hmg_on_hmg"] = None if hh is None else round(hh, 1)
                    row[tag + "_alpha_prot_on_hmg"] = None if "MAT1-1-1" not in sc else round(sc["MAT1-1-1"]["MAT1-2-1"], 1)
                    row[tag + "_hmg_prot_on_alpha"] = None if "MAT1-2-1" not in sc else round(sc["MAT1-2-1"]["MAT1-1-1"], 1)
                    row[tag + "_call"] = pc; row[tag + "_margin"] = mg
                # Pfam domains on the core proteins and on the winning optional gene
                for k, p in prot.items():
                    h = pfam_hits(pfam, p)
                    row["pfam_" + k] = ";".join(f"{n}:{round(s, 1)}" for n, s in sorted(h.items()))
            rows.append(row)
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(a.out / "loci_audit.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", restval="")
        w.writeheader(); w.writerows(rows)
    gcols = []
    for r in grows:
        for k in r:
            if k not in gcols:
                gcols.append(k)
    with open(a.out / "genome_status.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, gcols, delimiter="\t", restval="")
        w.writeheader(); w.writerows(grows)
    print(len(grows), "genomes;", len(rows), "loci")


main()
