##################################################
# predict cell types
# requires running annotate.R first
##################################################

# load libraries
library(tidyverse)
library(Seurat)
library(qs)
library(scMisc)

# function to map project query on ref and make predictions based on Seurat integration
mapSeurat <- function(ref, query) {
  reference_list <- list(ref = ref, query = query)
  features <- SelectIntegrationFeatures(object.list = reference_list)
  anchors <- FindTransferAnchors(
    reference = reference_list$ref,
    query = reference_list$query,
    normalization.method = "LogNormalize",
    features = features
  )
  predictions <- TransferData(
    anchorset = anchors,
    refdata = reference_list$ref$cluster
  )
  return(predictions)
}

# function to store predictions in seurat object
storePred <- function(predictions, label_col, score_col, seu_obj) {
  predictions_prep <-
    predictions |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(predicted.id, prediction.score.max, barcode) |>
    dplyr::mutate(
      predicted.id = ifelse(prediction.score.max < 0.3, "unknown", predicted.id)
    ) |>
    tibble::as_tibble() |>
    dplyr::rename(
      {{ label_col }} := predicted.id,
      {{ score_col }} := prediction.score.max
    )

  seu_obj@meta.data <-
    seu_obj@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(predictions_prep, by = "barcode") |>
    tibble::column_to_rownames(var = "barcode")

  return(seu_obj)
}

# load this dataset and subset ---
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)
cd8_nk <- subset(sc_merge, subset = cluster %in% c("CD8_NK"))

# map to tabula sapiens blood ----
cd8_terekhova <- qread(file.path(
  "objects",
  "conventional_cd8_seurat_object.qs"
))

scMisc::lss()

predictions_cd8_terehkova <- mapSeurat(ref = cd8_terekhova, query = cd8_nk)

cd8_nk <- storePred(
  predictions_cd8_terehkova,
  label_col = "cd8_terehkova_label",
  score_col = "cd8_terekhova_score",
  seu_obj = cd8_nk
)

pred_plot_cd8_nk <-
  DimPlot(
    cd8_nk,
    reduction = "umap.stacas.ss.all",
    group.by = "cd8_terehkova_label",
    raster = FALSE,
    pt.size = .5,
    alpha = .5,
    cols = colors_dutch,
    label = TRUE
  ) +
  theme_rect()

ggsave(
  plot = pred_plot_cd8_nk,
  file.path("results", "map", "map_cd8_nk_terekhova.png"),
  width = 20,
  height = 7
)
