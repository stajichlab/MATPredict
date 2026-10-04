# MATPredict container image.
# Build stage: solve nothing, install the locked pixi environment (pixi.lock).
# Runtime stage: the environment, the package source and the curated database.
#
#   docker build -t matpredict .
#   docker run --rm -v "$PWD":/data -e MATPREDICT_NCBI_EMAIL=you@example.org \
#     matpredict detect --genome /data/genome.fna --taxid 4837 --out-dir /data/out

FROM ghcr.io/prefix-dev/pixi:0.71.3-noble AS build
WORKDIR /app
COPY pixi.toml pixi.lock pyproject.toml README.md ./
COPY src ./src
RUN pixi install --locked -e default \
 && pixi shell-hook -e default -s bash > /shell-hook.sh \
 && echo 'exec "$@"' >> /shell-hook.sh

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
ENV MATPREDICT_DB_ROOT=/app/db \
    MATPREDICT_CACHE_DIR=/data/.matpredict_cache
WORKDIR /data
ENTRYPOINT ["/bin/bash", "/shell-hook.sh", "matpredict"]
CMD ["--help"]
