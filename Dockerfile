# syntax=docker/dockerfile:1.7

ARG BASE_IMAGE=mihem/csf_pbmc_pnp_seq@sha256:27ffba518df05e34788d58bc2abf5826b6b28e416866a406d8c9d42b388e00ef
FROM ${BASE_IMAGE}

ARG TARGETARCH

WORKDIR /csf_pbmc_pnp_seq

COPY renv.lock ./renv.lock

# Reconcile only the lockfile delta while preserving the separately installed
# Azimuth PBMC reference data package.
RUN --mount=type=cache,id=csf-pbmc-pnp-renv-r4.4-${TARGETARCH},target=/renv-cache,sharing=locked \
 R -e 'renv::restore(library = renv::paths$library(), clean = TRUE, exclude = "pbmcref.SeuratData", prompt = FALSE)'
