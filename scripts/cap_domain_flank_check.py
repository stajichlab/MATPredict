"""Domain content and flank synteny of the loci the weak-fragment cap would lower (#341).

For each locus in `review.tsv` (from `scripts/cap_annotation_review.py`), take the called idiomorph's strongest core gene (the gene the cap
ranks), translate the genome around its model in six frames (stop-free stretches of at least `--min-aa` residues; a spliced gene can be
split over ORFs) and search the ORFs with the Pfam models that name the gene's domain (`hmmsearch` gathering cutoffs, via pyhmmer):

* MAT1-1-1: PF04769 (MATalpha_HMGbox, the alpha-box);  MAT1-2-1 and MAT1-1-3: the HMG-box models PF00505, PF09011, PF14887;
* other core genes (MAT1-1-4, MAT1-1-5, MAT1-2-7 ...): no single domain is expected, so `domain` is `n/a`.

Flank synteny comes from the report's own flanking-gene hits on the same contig (SLA2, APN2, COX13): signed distance from the core gene
(negative = left of it); `flank_ok` needs at least one within `--flank-bp` (default 30,000).

A domain hit says the region encodes an HMG- or alpha-box protein, NOT that it is MAT1-2-1: an HMG-box paralog has the same domain. Flank
synteny is what ties the gene to the MAT locus. `verdict`: `domain+flank`, `domain_only`, `flank_only`, `neither`, or (for genes with
no expected domain) `flank_only` / `neither` with `domain` = `n/a`. Both are measurements of the model's surroundings, not of the call.
"""
from __future__ import annotations

import csv
import gzip
import sys
import tarfile
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from region_pfam_check import orfs  # noqa: E402

DOMAIN_ACCESSIONS = {"PF04769": "alpha-box", "PF00505": "HMG-box", "PF09011": "HMG-box", "PF14887": "HMG-box"}
EXPECTED = {"MAT1-1-1": "alpha-box", "MAT1-2-1": "HMG-box", "MAT1-1-3": "HMG-box"}
FLANK_GENES = ("SLA2", "APN2", "COX13")


def expected_domain(gene: str) -> str | None:
    return EXPECTED.get(gene)


def flank_offsets(locus: dict, gene: dict) -> list[tuple[str, int]]:
    """(flank gene, signed distance in bp) for the flanking genes on the core gene's contig; negative is to its left."""
    out = []
    for e in locus.get("gene_evidence") or []:
        if e["gene"] in FLANK_GENES and not e.get("lost_to") and e["contig"] == gene["contig"]:
            if e["end"] < gene["start"]:
                out.append((e["gene"], -(gene["start"] - e["end"])))
            elif e["start"] > gene["end"]:
                out.append((e["gene"], e["start"] - gene["end"]))
            else:
                out.append((e["gene"], 0))
    return sorted(out, key=lambda t: abs(t[1]))


def read_contig(fa: Path, contig: str) -> str:
    op = gzip.open if str(fa).endswith(".gz") else open
    buf, on = [], False
    with op(fa, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if on:
                    break
                on = line[1:].split()[0] == contig
            elif on:
                buf.append(line.strip())
    return "".join(buf)


def find_genome(genome_dir: Path, asm: str) -> Path | None:
    for p in sorted(genome_dir.glob(asm + "*")):
        n = p.name
        if n.endswith((".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna")) and ".masked" not in n:
            return p
    return None


def load_hmms(pfam: Path):
    import pyhmmer
    keep = []
    with pyhmmer.plan7.HMMFile(str(pfam)) as hf:
        for hmm in hf:
            acc = (hmm.accession or b"").decode().split(".")[0]
            if acc in DOMAIN_ACCESSIONS:
                keep.append((acc, hmm))
    return keep


def domain_hits(seq: str, offset: int, hmms, min_aa: int) -> list[dict]:
    """Pfam hits (GA cutoff) on the six-frame ORFs of `seq`; `offset` is the 1-based contig position of seq[0]."""
    import pyhmmer
    from pyhmmer.easel import Alphabet, TextSequence
    alpha = Alphabet.amino()
    orf_list = orfs(seq, min_aa)
    seqs = [TextSequence(name=f"o{i}".encode(), sequence=p).digitize(alpha) for i, (_, _, _, p) in enumerate(orf_list)]
    out = []
    for acc, hmm in hmms:
        for hits in pyhmmer.hmmsearch([hmm], seqs, bit_cutoffs="gathering"):
            for h in hits:
                if not h.included:
                    continue
                fr, o0, o1, prot = orf_list[int(h.name.decode()[1:])]
                out.append({"acc": acc, "domain": DOMAIN_ACCESSIONS[acc], "score": round(h.score, 1), "evalue": h.evalue,
                            "orf_start": offset + o0, "orf_end": offset + o1 - 1, "aa": len(prot)})
    return out


def check_locus(locus: dict, review: dict, seq_of, hmms, *, pad=600, min_aa=40, flank_bp=30000) -> dict:
    gene_name = review["strongest_own_gene"]
    gene = next(e for e in locus["gene_evidence"] if e["gene"] == gene_name and not e.get("lost_to"))
    flanks = flank_offsets(locus, gene)
    flank_ok = any(abs(d) <= flank_bp for _, d in flanks)
    want = expected_domain(gene_name)
    hits, note = [], ""
    if want:
        contig_seq = seq_of(gene["contig"])
        if contig_seq is None or len(contig_seq) < gene["end"]:
            note = "genome or contig not found"
        else:
            s = max(1, gene["start"] - pad)
            hits = domain_hits(contig_seq[s - 1:gene["end"] + pad], s, hmms, min_aa)
            hits = [h for h in hits if h["orf_end"] >= gene["start"] and h["orf_start"] <= gene["end"]]
    domain = "n/a" if not want else ("found" if any(h["domain"] == want for h in hits) else ("not_checked" if note else "absent"))
    if domain == "n/a":
        verdict = "flank_only" if flank_ok else "neither"
    elif domain == "not_checked":
        verdict = "not_checked"
    else:
        verdict = {(True, True): "domain+flank", (True, False): "domain_only", (False, True): "flank_only", (False, False): "neither"}[(domain == "found", flank_ok)]
    return {"gene": gene_name, "expected_domain": want or "", "domain": domain,
            "domain_hits": ";".join(f'{h["acc"]}:{h["score"]}:{h["orf_start"]}-{h["orf_end"]}' for h in hits),
            "flanks": ";".join(f"{g}:{d}" for g, d in flanks), "flank_ok": flank_ok, "verdict": verdict, "note": note}


def _iter_tar(stream):
    with tarfile.open(fileobj=stream, mode="r|") as tf:
        for m in tf:
            if m.isfile() and m.name.endswith("detection_report.yaml"):
                yield Path(m.name).parent.name, yaml.safe_load(tf.extractfile(m).read()) or {}


def run(docs, review_tsv: Path, genome_dir: Path, hmms, out: Path, **kw) -> dict:
    review = {(r["genome"], r["contig"], int(r["start"]), int(r["end"])): r for r in csv.DictReader(review_tsv.open(), delimiter="\t")}
    rows, counts = [], {}
    cache: dict[str, tuple[str, str] | None] = {}
    for genome, doc in docs:
        for locus in doc.get("detected") or []:
            key = (genome, locus.get("contig"), locus.get("start"), locus.get("end"))
            if key not in review or locus.get("family") != "Ascomycota:MAT":
                continue
            fa = find_genome(genome_dir, genome)
            def seq_of(contig, fa=fa):
                if fa is None:
                    return None
                if (str(fa), contig) not in cache:
                    cache[(str(fa), contig)] = read_contig(fa, contig) or None
                return cache[(str(fa), contig)]
            r = check_locus(locus, review[key], seq_of, hmms, **kw)
            r["verdict_key"] = r["verdict"]
            counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
            rows.append([genome, *key[1:], locus["idiomorph"], locus["confidence"], review[key]["review"], r["gene"], r["expected_domain"], r["domain"],
                         r["domain_hits"], r["flanks"], r["flank_ok"], r["verdict"], r["note"]])
    # the loci with no annotated gene first, then by genome
    rows.sort(key=lambda r: (0 if r[6] == "no_annotated_gene" else 1, r[0], r[1], r[2]))
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "domain_flank.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "contig", "start", "end", "idiomorph", "tier", "annotation_review", "strongest_own_gene", "expected_domain", "domain",
                    "domain_hits(acc:score:start-end)", "flanks(gene:distance_bp)", "flank_ok", "verdict", "note"])
        w.writerows(rows)
    with open(out / "summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["verdict", "loci"])
        for k, v in sorted(counts.items()):
            w.writerow([k, v])
    return counts


def main(argv: list[str]) -> int:
    src, review, genome_dir, pfam, out = argv
    docs = _iter_tar(sys.stdin.buffer if src == "-" else open(src, "rb"))
    print(run(docs, Path(review), Path(genome_dir), load_hmms(Path(pfam)), Path(out)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
