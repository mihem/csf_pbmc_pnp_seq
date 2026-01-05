##################################################
# use Azimuth to predict cell types
# requires running preprocess.R first
##################################################

# load libraries
library(tidyverse)
library(Seurat)
library(Azimuth)
library(qs)
library(scMisc)

# general settings  ----
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)

# load data ----
# read preprocessed data from preprocess.R
sc_merge <- JoinLayers(sc_merge)
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

# run Azimuth ----
sc_merge_azimuth <- RunAzimuth(
  query = sc_merge,
  reference = "pbmcref",
  assay = "RNA"
)

# prepare predictions ----
predictions_pbmcref <-
  data.frame(
    azimuth_pbmcref1 = sc_merge_azimuth$predicted.celltype.l1,
    azimuth_pbmcref1_score = sc_merge_azimuth$predicted.celltype.l1.score,
    azimuth_pbmcref2 = sc_merge_azimuth$predicted.celltype.l2,
    azimuth_pbmcref2_score = sc_merge_azimuth$predicted.celltype.l2.score,
    azimuth_pbmcref3 = sc_merge_azimuth$predicted.celltype.l3,
    azimuth_pbmcref3_score = sc_merge_azimuth$predicted.celltype.l3.score
  ) |>
  mutate(
    azimuth_pbmcref1 = ifelse(
      azimuth_pbmcref1_score < 0.4,
      "unknown",
      azimuth_pbmcref1
    )
  ) |>
  mutate(
    azimuth_pbmcref2 = ifelse(
      azimuth_pbmcref2_score < 0.4,
      "unknown",
      azimuth_pbmcref2
    )
  ) |>
  mutate(
    azimuth_pbmcref3 = ifelse(
      azimuth_pbmcref3_score < 0.4,
      "unknown",
      azimuth_pbmcref3
    )
  )

# add to seurat object ---
sc_merge <- AddMetaData(sc_merge, predictions_pbmcref)

qsave(sc_merge, file.path("objects", "sc_merge.qs"))
