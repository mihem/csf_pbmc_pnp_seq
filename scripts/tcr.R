##################################################
# analyze T cell receptor (TCR) data
# requires running annotate.R first
##################################################

# load libraries ----
library(qs)
library(Seurat)
library(tidyverse)
library(scRepertoire)

# load filtered contig annotations ---
tcr_files <- list.files(
    file.path("raw", "tcr"),
    pattern = "filtered_contig_annotations.csv",
    full.names = TRUE,
    recursive = TRUE
)

contig_list <- lapply(tcr_files, read.csv)
names(contig_list) <- gsub(
    "raw/tcr/(.*)/outs/filtered_contig_annotations.csv",
    "\\1",
    tcr_files
)

str(contig_list, max.level = 1)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

str(sc_merge@meta.data)


