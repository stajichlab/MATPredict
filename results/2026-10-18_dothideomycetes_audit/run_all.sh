#!/bin/bash
# Rebuild every table from loci_audit.tsv / genome_status.tsv (+ scan/ and bioproject.tsv when present).
# PY = python with pandas and pyyaml.  Intermediates written next to the inputs: genome_audit.tsv, genome_audit_verdict.tsv
# (safe to delete; they are regenerated).  Order matters: analyze -> tables -> truth -> scan_integrate -> finalize -> the rest.
set -e
cd "$(dirname "$0")"
PY="${PY:-python3}"
"$PY" analyze.py .
"$PY" tables.py . > /dev/null
"$PY" truth.py . > /dev/null
[ -d scan ] && "$PY" scan_integrate.py . | tail -3
"$PY" finalize.py . | tail -30
"$PY" summary.py . > /dev/null
"$PY" suspects.py . > /dev/null
"$PY" paralog_check.py . > /dev/null
[ -f bioproject.tsv ] && "$PY" species_skew_by_bioproject.py . > /dev/null
"$PY" note_tables.py .
echo rebuilt
