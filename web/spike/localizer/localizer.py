"""tblastn and makeblastdb stand-ins built on MMseqs2's module steps.

`mmseqs search` is a shell workflow (execvp), which cannot run in WebAssembly, so the steps it would
run are called directly, with the arguments captured from a real run. Output is the pipeline's
11-column outfmt (search.py _TBLASTN_OUTFMT). Minus-strand hits have sstart > send and sframe -1.

`run` is any callable taking an argv list and raising CalledProcessError on failure. The command
line uses subprocess; the browser bridge uses its own runner. `mmseqs` is the argv prefix for the
MMseqs2 binary (["mmseqs"] or ["node", "mmseqs-mem.js"]).
"""
import os
import shutil
import tempfile

PREFILTER_SUB = "b64:YWE6Ymxvc3VtNjIub3V0LG51Y2w6bnVjbGVvdGlkZS5vdXQ="
PREFILTER_SEED = "b64:YWE6VlRNTDgwLm91dCxudWNsOm51Y2xlb3RpZGUub3V0"
OUTFMT = ("query,target,pident,alnlen,qlen,qstart,qend,tstart,tend,evalue,bits")


def makeblastdb(argv, run=None):
    """Stand-in for makeblastdb: the pipeline only needs the genome FASTA behind the prefix."""
    src = argv[argv.index("-in") + 1]
    prefix = argv[argv.index("-out") + 1]
    shutil.copy(src, prefix + ".fa")
    return ""


def tblastn(argv, mmseqs, run, threads="1", max_seqs="5000", s="7.5", e="10"):
    """Run the MMseqs2 steps for a protein query against a translated genome; return outfmt 6 text."""
    query = argv[argv.index("-query") + 1]
    db = argv[argv.index("-db") + 1]
    with tempfile.TemporaryDirectory() as tmp:
        p = lambda name: os.path.join(tmp, name)
        q, t, res = p("q"), p("t"), p("res")
        work = os.path.join(tmp, "work")
        os.makedirs(os.path.join(work, "search"), exist_ok=True)
        orfs, pref, aln, conv = p("t_orfs_aa"), os.path.join(work, "search", "pref_0"), p("aln"), p("conv.tsv")
        m = lambda *a: run(list(mmseqs) + list(a))
        m("createdb", query, q)
        m("createdb", db + ".fa", t, "--dbtype", "2")
        m("extractorfs", t, orfs, "--min-length", "30", "--max-length", "32734", "--max-gaps", "2147483647",
          "--contig-start-mode", "2", "--contig-end-mode", "2", "--orf-start-mode", "1",
          "--forward-frames", "1,2,3", "--reverse-frames", "1,2,3", "--translation-table", "1",
          "--translate", "1", "--use-all-table-starts", "0", "--id-offset", "0", "--create-lookup", "0",
          "--threads", threads, "--compressed", "0", "-v", "3")
        m("prefilter", q, orfs, pref, "--sub-mat", PREFILTER_SUB, "--seed-sub-mat", PREFILTER_SEED,
          "-k", "0", "--target-search-mode", "0", "--k-score", "seq:2147483647,prof:2147483647",
          "--alph-size", "aa:21,nucl:5", "--max-seq-len", "65535", "--max-seqs", max_seqs,
          "--split", "0", "--split-mode", "2", "--split-memory-limit", "0", "-c", "0", "--cov-mode", "0",
          "--comp-bias-corr", "1", "--comp-bias-corr-scale", "1", "--diag-score", "1",
          "--exact-kmer-matching", "0", "--mask", "1", "--mask-prob", "0.9", "--mask-lower-case", "0",
          "--mask-n-repeat", "0", "--min-ungapped-score", "15", "--add-self-matches", "0",
          "--spaced-kmer-mode", "1", "--db-load-mode", "0", "--pca", "substitution:1.100,context:1.400",
          "--pcb", "substitution:4.100,context:5.800", "--threads", threads, "--compressed", "0", "-v", "3",
          "-s", s)
        m("align", q, orfs, pref, aln, "--sub-mat", PREFILTER_SUB, "-a", "0", "--alignment-mode", "2",
          "--alignment-output-mode", "0", "--wrapped-scoring", "0", "-e", e,
          "--min-seq-id", "0", "--min-aln-len", "0", "--seq-id-mode", "0", "--alt-ali", "0", "-c", "0",
          "--cov-mode", "0", "--max-seq-len", "65535", "--comp-bias-corr", "1", "--comp-bias-corr-scale", "1",
          "--max-rejected", "2147483647", "--max-accept", "2147483647", "--add-self-matches", "0",
          "--db-load-mode", "0", "--pca", "substitution:1.100,context:1.400", "--pcb", "substitution:4.100,context:5.800",
          "--score-bias", "0", "--realign", "0", "--realign-score-bias", "-0.2", "--realign-max-seqs", "2147483647",
          "--corr-score-weight", "0", "--gap-open", "aa:11,nucl:5", "--gap-extend", "aa:1,nucl:2", "--zdrop", "40",
          "--threads", threads, "--compressed", "0", "-v", "3")
        m("offsetalignment", q, q, t, orfs, aln, res, "--chain-alignments", "0", "--merge-query", "1",
          "--search-type", "2", "--threads", threads, "--compressed", "0", "--db-load-mode", "0", "-v", "3")
        m("convertalis", q, t, res, conv, "--format-output", OUTFMT)
        lines = []
        with open(conv) as fh:
            for line in fh:
                qn, tn, pid, alnlen, qlen, qs, qe, ts, te, ev, bits = line.rstrip("\n").split("\t")
                qlen = int(qlen)
                qcov = 100.0 * (int(qe) - int(qs) + 1) / qlen
                sframe = 1 if int(ts) < int(te) else -1
                lines.append("\t".join([qn, tn, pid, alnlen, str(qlen), f"{qcov:.3f}", ev, bits,
                                        ts, te, str(sframe)]))
        return "".join(line + "\n" for line in lines)
