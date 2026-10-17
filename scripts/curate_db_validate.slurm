#!/usr/bin/bash
#SBATCH -N 1 -c 1 --mem 8gb --time 00:30:00 -J curate-validate
# `matpredict curate-db validate` on a compute node with at least 8 GB (its taxonkit step has been killed on a shared
# login node). usage: sbatch scripts/curate_db_validate.slurm PHYLUM RECORD_ID
# Set REPO to the checkout whose db/candidates holds the record (default: the submit directory). The NCBI e-mail comes
# from ~/.config/matpredict/config.toml or $MATPREDICT_NCBI_EMAIL, as for any curate-db command. Set PY to the Python of an
# environment that has MATPredict's dependencies and taxonkit on its PATH (default: python).
set -euo pipefail
PHYLUM="${1:?usage: curate_db_validate.slurm PHYLUM RECORD_ID}"
RECORD="${2:?usage: curate_db_validate.slurm PHYLUM RECORD_ID}"
cd "${REPO:-$SLURM_SUBMIT_DIR}"
PYTHONPATH="$PWD/src" MATPREDICT_DB_ROOT="$PWD/db" "${PY:-python}" -m MATPredict curate-db validate --phylum "$PHYLUM" --record-id "$RECORD"
