# r-base 4.4.0 with ubuntu noble
FROM posit/r-base@sha256:f5b4a7e4d6dd126d3768d04b11a2ffad45df591dc94a286490c9bb2b74763d98

# system libraries
# Try to only install system libraries you actually need
# Package Manager is a good resource to help discover system deps
RUN apt-get update --yes \
 && apt-get upgrade --yes \
 && apt-get install --yes \
libglpk-dev libxml2-dev git libhdf5-dev \
libgsl-dev cmake libmagick++-dev

# install R packages required 
# Change the packages list to suit your needs
RUN R -e 'install.packages("https://packagemanager.posit.co/cran/latest/src/contrib/Archive/renv/renv_1.0.8.tar.gz")'
# RUN R -e 'install.packages("https://packagemanager.posit.co/cran/latest/src/contrib/Archive/pak/pak_0.7.0.tar.gz")'

# Enable pak for renv and use 6 cores
# RUN echo "options(Ncpus = 6, renv.config.pak.enabled = TRUE)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
# pak fails to solve dependencies https://github.com/r-lib/pak/issues/724
RUN echo "options(Ncpus = 6)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site

# Copy renv files 
WORKDIR /csf_pbmc_pnp_seq
COPY renv.lock renv.lock
COPY .gitignore .gitignore

# Persist the renv library and cache directories across builds
# VOLUME /root/.local/share/renv/cache /root/.local/share/renv/library

ENV RENV_PATHS_LIBRARY=renv/library

# Restore the R environment
RUN R -e "renv::restore()"

#ARG QUARTO_VERSION="1.6.39"
#RUN curl -o quarto-linux-amd64.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb
#RUN dpkg -i quarto-linux-amd64.deb || apt-get install -f -y
#RUN rm quarto-linux-amd64.deb

