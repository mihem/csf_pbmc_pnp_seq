##################################################
# PCA plots of cluster abundance across samples
# requires running annotate.R first
##################################################

# libraries ---
library(Seurat)
library(tidyverse)
library(scMisc)
library(qs)
library(FactoMineR)
library(factoextra)

# load helper functions ---
source(file.path("scripts", "pca_helper.R")) 

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

sc_merge_csf <- subset(sc_merge, tissue %in% "CSF")
sc_merge_pbmc <- subset(sc_merge, tissue %in% "PBMC")

sc_merge_csf_main <- subset(
  sc_merge_csf,
  subset = diagnosis %in% c("CIDP", "GBS", "CIAP", "CTRL")
)

sc_merge_pbmc_main <- subset(
  sc_merge_pbmc,
  subset = diagnosis %in% c("CIDP", "GBS", "CIAP", "CTRL")
)

# PCA for CSF only
pcaSeurat1(
  object_list = sc_merge_csf_main,
  cluster = "cluster",
  sample = "sample",
  patient = "patient",
  condition = "diagnosis",
  object_name = "csf",
  color_palette = sc_merge@misc$diagnosis_col,
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)

# PCA for PBMC only
pcaSeurat1(
  object_list = sc_merge_pbmc_main,
  cluster = "cluster",
  sample = "sample",
  patient = "patient",
  condition = "diagnosis",
  object_name = "pbmc",
  color_palette = sc_merge@misc$diagnosis_col,
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)

# PCA for combined CSF and PBMC
pcaSeurat1(
  object_list = list(sc_merge_csf_main, sc_merge_pbmc_main),
  cluster = "cluster",
  sample = "sample",
  patient = "patient",
  condition = "diagnosis",
  tissue = "tissue",
  object_name = "combined_csf_pbmc",
  color_palette = sc_merge@misc$diagnosis_col,
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)

# PCA for combined CSF and PBMC
pcaSeurat1(
  object_list = list(sc_merge_csf, sc_merge_pbmc),
  cluster = "cluster",
  sample = "sample",
  patient = "patient",
  condition = "diagnosis",
  tissue = "tissue",
  object_name = "combined_csf_pbmc",
  color_palette = sc_merge@misc$diagnosis_col,
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)
