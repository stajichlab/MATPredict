# MATPredict container image.
# Build stage: solve nothing, install the locked pixi environment (pixi.lock).
# Runtime stage: the environment, the package source and the curated database.
#
#   docker build -t matpredict .
#   docker run --rm -v "$PWD":/data \
#     matpredict detect --genome /data/genome.fna --taxid 4837 --out-dir /data/out
# Taxonomy lookups use the embedded NCBI snapshot (no network).

FROM ghcr.io/prefix-dev/pixi:0.71.3-noble AS build
WORKDIR /app
COPY pixi.toml pixi.lock pyproject.toml README.md ./
COPY src ./src
RUN pixi install --locked -e default \
 && pixi shell-hook -e default -s bash > /shell-hook.sh \
 && echo 'exec "$@"' >> /shell-hook.sh

# Offline NCBI taxonomy (curator ruling 2026-10-04): a slim all-taxa table
# (taxid, parent, rank, genetic code, scientific name; merged ids) built from a
# DATED NCBI archive, so one image = one taxonomy snapshot and runs make no
# E-utilities call. NCBI publishes no checksum for the archive; the SHA-256 is
# the one measured when the date was pinned. Change both together.
ARG TAXDUMP_DATE=2026-10-01
ARG TAXDUMP_SHA256=d744af371c0b9fc7269d80b49546b6ac4c9ddc52a2fd97bcbfd02783edb5c9eb
RUN /app/.pixi/envs/default/bin/python -c "import urllib.request,sys; urllib.request.urlretrieve(sys.argv[1], '/tmp/taxdmp.zip')" \
      "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${TAXDUMP_DATE}.zip" \
 && echo "${TAXDUMP_SHA256}  /tmp/taxdmp.zip" | sha256sum -c - \
 && PYTHONPATH=/app/src /app/.pixi/envs/default/bin/python -m MATPredict curate-db build-taxonomy \
      --taxdump /tmp/taxdmp.zip --snapshot "${TAXDUMP_DATE}" --out /app/taxonomy/ncbi_taxonomy.tsv.zst \
 && rm /tmp/taxdmp.zip

FROM ubuntu:24.04 AS runtime
ARG VERSION=dev
LABEL org.opencontainers.image.title="MATPredict" \
      org.opencontainers.image.description="Find and type fungal mating-type (MAT) loci in genome assemblies" \
      org.opencontainers.image.source="https://github.com/stajichlab/MATPredict" \
      org.opencontainers.image.version="${VERSION}"
# The environment keeps its build path (/app/.pixi/envs/default): conda
# packages and the editable MATPredict install record absolute paths.
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /shell-hook.sh /shell-hook.sh
COPY --from=build /app/src /app/src
COPY pyproject.toml README.md /app/
COPY db /app/db
COPY --from=build /app/taxonomy /app/taxonomy
# MATPREDICT_OFFLINE=1: a taxid not in the snapshot is reported (routing_error)
# instead of fetched from NCBI. Run with -e MATPREDICT_OFFLINE=0 to allow
# E-utilities for such taxids.
ENV MATPREDICT_DB_ROOT=/app/db \
    MATPREDICT_TAXONOMY=/app/taxonomy/ncbi_taxonomy.tsv.zst \
    MATPREDICT_OFFLINE=1 \
    MATPREDICT_CACHE_DIR=/data/.matpredict_cache
WORKDIR /data
ENTRYPOINT ["/bin/bash", "/shell-hook.sh", "matpredict"]
CMD ["--help"]
