"""Per-gene table of the loci that the partial-coverage rule touched (run on the HPCC, pixi env).

    extract_loci.py RUNS_ON_DIR JOBS_TSV NOTES_TSV PROFILES_DIR OUT_GENES.tsv OUT_LOCI.tsv GENOME_ID [GENOME_ID...]

For each genome: the locus named in vote_notes_on.tsv (contig, start); every gene_evidence entry with its model
(coverage, identity, reference, status, bitscore, lost_to), the protein translated from the exons (or the hit span
in the best of three frames), its score on the alpha-box and HMG-box profiles, internal stops, contig length, distance
of the gene to the nearest contig end, and N-gap content within 3 kb of the locus.
"""
import csv, gzip, sys
from pathlib import Path
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from score_core_profiles import load_hmms, coding_sequence, best_scores
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
runs, jobs, notes, prof, out_g, out_l, *ids = sys.argv[1:]
fa = {r[0]: r[4] for r in csv.reader(open(jobs), delimiter="\t")}
nt = {}
for r in csv.DictReader(open(notes), delimiter="\t"):
    nt[r["genome"]] = r
hmms = load_hmms(Path(prof))
gw = csv.writer(open(out_g, "w"), delimiter="\t", lineterminator="\n")
lw = csv.writer(open(out_l, "w"), delimiter="\t", lineterminator="\n")
gw.writerow("genome contig gene role start end strand identity coverage reference status method bitscore lost_to model_aa internal_stops alpha_best_of_models hmg_best_of_models dist_contig_end".split())
lw.writerow("genome contig contig_len locus_start locus_end n_loci_mat label confidence n_N_3kb genes_found".split())
for gid in ids:
    n = nt.get(gid)
    doc = yaml.load(open(f"{runs}/{gid}/detection_report.yaml"), Loader=L)
    loci = [c for c in doc["detected"] if c["family"] == "Ascomycota:MAT"]
    if n:
        loc = [c for c in loci if c["contig"] == n["contig"] and str(c["start"]) == n["start"]]
    else:
        loc = sorted(loci, key=lambda c: -len(c["genes_found"]))[:1]
    if not loc:
        print("no locus", gid); continue
    c = loc[0]
    opener = gzip.open if fa[gid].endswith(".gz") else open
    seq = None
    with opener(fa[gid], "rt") as fh:
        for r in SeqIO.parse(fh, "fasta"):
            if r.id == c["contig"]:
                seq = str(r.seq).upper(); break
    clen = len(seq)
    nN = seq[max(0, c["start"] - 3000):c["end"] + 3000].count("N")
    lw.writerow([gid, c["contig"], clen, c["start"], c["end"], len(loci), c["idiomorph"], c["confidence"], nN, ",".join(c["genes_found"])])
    for g in c["gene_evidence"]:
        prot = []
        if g.get("exons"):
            ex = tuple(sorted((int(x["start"]), int(x["end"])) for x in g["exons"]))
            dna = coding_sequence(ex, g["strand"], seq)
            prot = [str(Seq(dna[: len(dna) // 3 * 3]).translate())]
        else:
            dna = coding_sequence(((int(g["start"]), int(g["end"])),), g["strand"], seq)
            prot = [str(Seq(dna[f:][: (len(dna) - f) // 3 * 3]).translate()) for f in range(3)]
        stops = min(p.rstrip("*").count("*") for p in prot) if prot else ""
        pc = [p.replace("*", "X") for p in prot]
        sc = best_scores(hmms, pc) if g["role"] == "core_MAT" else {"MAT1-1-1": "", "MAT1-2-1": ""}
        am = g.get("alternate_model")
        if am and am.get("exons") and g["role"] == "core_MAT":
            ex2 = tuple(sorted((int(x["start"]), int(x["end"])) for x in am["exons"]))
            d2 = coding_sequence(ex2, am["strand"], seq)
            sa = best_scores(hmms, [str(Seq(d2[: len(d2) // 3 * 3]).translate()).replace("*", "X")])
            sc = {k: max(sc[k], sa[k]) for k in sc}
        gw.writerow([gid, c["contig"], g["gene"], g["role"], g["start"], g["end"], g["strand"], g.get("identity"), g.get("coverage"),
                     g.get("reference_record"), g["status"], g["method"], g.get("bitscore"), g.get("lost_to") or "",
                     len(prot[0].rstrip("*")) if g.get("exons") else "", stops,
                     round(sc["MAT1-1-1"], 1) if sc["MAT1-1-1"] != "" else "", round(sc["MAT1-2-1"], 1) if sc["MAT1-2-1"] != "" else "",
                     min(int(g["start"]) - 1, clen - int(g["end"]))])
