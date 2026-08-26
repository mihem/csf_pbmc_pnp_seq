# syntax=docker/dockerfile:1.7

# r-base 4.4.0 with ubuntu noble
FROM posit/r-base@sha256:f5b4a7e4d6dd126d3768d04b11a2ffad45df591dc94a286490c9bb2b74763d98

ARG TARGETARCH

# system libraries
RUN apt-get update --yes \
 && apt-get install --yes --no-install-recommends \
 libglpk-dev libxml2-dev git libhdf5-dev \
 libgsl-dev cmake libmagick++-dev \
 && rm -rf /var/lib/apt/lists/*

# Bootstrap with renv 1.0.8 so GitHub packages in subdirectories restore
# correctly. The project library still restores renv 1.0.7 from renv.lock.
RUN R -e 'install.packages("https://packagemanager.posit.co/cran/latest/src/contrib/Archive/renv/renv_1.0.8.tar.gz")'

RUN echo "options(Ncpus = 6)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site

WORKDIR /csf_pbmc_pnp_seq

ENV RENV_PATHS_LIBRARY=/opt/renv/library \
    RENV_PATHS_CACHE=/renv-cache \
    RENV_CONFIG_CACHE_SYMLINKS=FALSE \
    R_LIBS_USER=/opt/renv/library/linux-ubuntu-noble/R-4.4/x86_64-pc-linux-gnu

COPY renv.lock ./renv.lock

# BuildKit keeps compiled package artifacts between lockfile changes. Packages
# are copied into the image library because the cache mount is build-time only.
RUN --mount=type=cache,id=csf-pbmc-pnp-renv-r4.4-${TARGETARCH},target=/renv-cache,sharing=locked \
 R -e 'renv::restore(library = renv::paths$library(), prompt = FALSE)'

# Install the exact Azimuth PBMC reference during the image build. This avoids
# package installation and network access from the analysis pipeline itself.
RUN R -e 'SeuratData::InstallData("pbmcref", lib = renv::paths$library())'
