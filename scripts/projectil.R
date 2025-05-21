###############################################################################
# use ProjecTILs to project TILs from the reference map to the query
# requires running annotate.R first
###############################################################################

# load libraries ----
library(qs)
library(Seurat)
library(tidyverse)
library(ProjecTILs)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)
cd8_nk <- subset(sc_merge, subset = cluster %in% c("CD8_NK"))
cd8_nk_ciap <- subset(cd8_nk, subset = diagnosis %in% c("CIAP"))

table(cd8_nk_ciap$diagnosis)


ref <- load.reference.map()

#project reference on t cells
cd8_nk_projected <- make.projection(
    cd8_nk,
    ref = ref,
    ncores = 8,
    filter.cells = FALSE
)
#project reference on t cells
cd8_nk_ciap_projected <- make.projection(
    cd8_nk_ciap,
    ref = ref,
    ncores = 8,
    filter.cells = FALSE
)

plot.projection(
    ref,
    cd8_nk_projected,
    linesize = 0.5,
    pointsize = 0.5
)

ggsave(
    file.path("results", "projectil", "cd8_nk_projectil.pdf"),
    width = 8,
    height = 6
)

plot.projection(
    ref,
    cd8_nk_ciap_projected,
    linesize = 0.5,
    pointsize = 0.5
)

ggsave(
    file.path("results", "projectil", "cd8_nk_ciap_projectil.pdf"),
    width = 8,
    height = 6
)

