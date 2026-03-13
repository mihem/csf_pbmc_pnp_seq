##################################################
# predict cell types
# requires running annotate.R first
##################################################

# load libraries
library(tidyverse)
library(Seurat)
library(qs)
library(scMisc)
library(writexl)

# colors from https://romanhaa.github.io/projects/scrnaseq_workflow/
colors_dutch <- c(
  '#FFC312',
  '#C4E538',
  '#12CBC4',
  '#FDA7DF',
  '#ED4C67',
  '#F79F1F',
  '#A3CB38',
  '#1289A7',
  '#D980FA',
  '#B53471',
  '#EE5A24',
  '#009432',
  '#0652DD',
  '#9980FA',
  '#833471',
  '#EA2027',
  '#006266',
  '#1B1464',
  '#5758BB',
  '#6F1E51'
)

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
# cd8_terekhova <- qread(file.path(
#   "objects",
#   "conventional_cd8_seurat_object.qs"
# ))

cd8_nk_terekhova <- qread(file.path(
  "objects",
  "conventional_cd8_nk_terekhova.qs"
))

Idents(cd8_nk_terekhova) <- cd8_nk_terekhova$cluster

cd8_nk_terekhova_small <- subset(
  cd8_nk_terekhova,
  downsample = 1000
)

predictions_cd8_nk_terehkova <- mapSeurat(
  ref = cd8_nk_terekhova_small,
  query = cd8_nk
)

cd8_nk <- storePred(
  predictions_cd8_nk_terehkova,
  label_col = "cd8_nk_terehkova_label",
  score_col = "cd8_nk_terekhova_score",
  seu_obj = cd8_nk
)

pred_plot_cd8_nk <-
  DimPlot(
    cd8_nk,
    reduction = "umap.stacas.ss.all",
    group.by = "cd8_nk_terehkova_label",
    raster = FALSE,
    pt.size = .5,
    alpha = .5,
    cols = colors_dutch,
    label = TRUE
  ) +
  theme_rect()

ggsave(
  plot = pred_plot_cd8_nk,
  file.path("results", "map", "map_cd8_nk_terekhova_cd8_nk.png"),
  width = 20,
  height = 7
)

dplyr::count(cd8_nk@meta.data, cd8_nk_terehkova_label) |>
  write_xlsx(
    file.path("results", "map", "map_cd8_nk_terekhova_cd8_nk.xlsx")
  )

qsave(
  cd8_nk,
  file.path("objects", "cd8_nk_terekhova.qs")
)
