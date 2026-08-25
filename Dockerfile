# syntax=docker/dockerfile:1.7

FROM mihem/csf_pbmc_pnp_seq:v0.6

WORKDIR /csf_pbmc_pnp_seq

COPY renv.lock ./renv.lock

# Reconcile only the lockfile delta while preserving the separately installed
# Azimuth PBMC reference data package.
RUN R -e 'renv::restore(library = renv::paths$library(), clean = TRUE, exclude = "pbmcref.SeuratData", prompt = FALSE)'
