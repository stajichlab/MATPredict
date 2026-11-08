#!/usr/bin/env python3
"""Study (BioProject) titles for the reads pages: one row per study accession that appears in the committed wave folders, with its title from the
PopSeq kingdom run table (study_accession, study_title). The reads site shows the title when the mouse is over a study link.

    study_titles.py PILOT_DIR POPSEQ_RUN_TABLE(.tsv/.csv, .gz) OUT.tsv

A study the run table does not have is written with an empty title (empty = unknown), so the page still links it. Standard library only.
"""
import csv, glob, gzip, io, os, re, sys

PILOT, RUNS, OUT = sys.argv[1:4]
studies = set()
for w in glob.glob(os.path.join(PILOT, "wave*")):
    f = os.path.join(w, "results_all_final.tsv")
    if re.fullmatch(r"wave\d+[a-z_]*", os.path.basename(w)) and os.path.exists(f):
        studies |= {r.get("study", "") for r in csv.DictReader(open(f, newline=""), delimiter="\t")}
studies.discard("")
fh = io.TextIOWrapper(gzip.open(RUNS, "rb"), newline="") if RUNS.endswith(".gz") else open(RUNS, newline="")
first = fh.readline(); delim = "\t" if first.count("\t") > first.count(",") else ","
cols = next(csv.reader([first], delimiter=delim))
title = {}
for r in csv.DictReader(fh, fieldnames=cols, delimiter=delim):
    s = r.get("study_accession", "")
    if s in studies and s not in title and r.get("study_title"):
        title[s] = " ".join(r["study_title"].split())
with open(OUT, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n"); w.writerow(["study", "title"])
    for s in sorted(studies): w.writerow([s, title.get(s, "")])
print(f"{len(studies)} studies; {len(title)} with a title; {len(studies) - len(title)} without")
