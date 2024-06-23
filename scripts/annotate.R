##################################################
# annotate data
# requires running integrate.R first
##################################################

library(tidyverse)
library(Seurat)
library(qs)
library(scMisc)
library(presto)
library(clustree)
library(writexl)

# general settings  ----
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)

# load data ----
# read preprocessed data from preprocess.R
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

# find clusters ---
sc_merge <- FindNeighbors(
  sc_merge,
  reduction = "stacas.ss.all",
  dims = 1:30,
  assay = "RNA",
  graph.name  = c("RNA_nn", "RNA_snn")
)

for (res in seq(from = 0.2, to = 1.2, by = 0.1)) {
  sc_merge <- FindClusters(sc_merge, resolution = res)
}

qs::qsave(sc_merge, file.path("objects", "sc_merge.qs"))

# plot clustering ---
clustree <- clustree(sc_merge, prefix = "RNA_snn_res.")

ggsave(plot = clustree, file.path("results", "umap", "stacas_ss_all_clustree.png"), width = 15, height = 15)

my_cols_100 <- unname(Polychrome::createPalette(100, pals::cols25()))
resolutions <- paste0("RNA_snn_res.", seq(from = 0.2, to = 1.2, by = 0.1))

umap_list <- list()
for (res in resolutions) {
    umap_list[[res]] <-
        DimPlot(sc_merge, reduction = "umap.stacas.ss.all", pt.size = .1, raster = FALSE, alpha = 0.1, group.by = res, cols = my_cols_100, label = TRUE) +
        theme_rect() +
        NoLegend()
}

umap_list <- patchwork::wrap_plots(umap_list, ncol = 2)
ggsave(plot = umap_list, file.path("results", "umap", "stacas_ss_all_resolutions.png"), width = 16, height = 40)

#  use resolution 1

# find markers helper function
findMarkers <- function(ident1, ident2 = NULL, object, only_pos, min_pct, logfc_threshold, assay = assay) {
  result <- Seurat::FindMarkers(object, ident.1 = ident1, ident.2 = ident2, min.pct = min_pct, logfc.threshold = logfc_threshold, only.pos = only_pos, assay = assay) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    dplyr::arrange(desc(avg_log2FC))
  return(result)
}


Idents(sc_merge) <- sc_merge$RNA_snn_res.1
topmarkers <-
  lapply(
    # unique(sc_merge@misc$cluster_order),
    unique(sc_merge$RNA_snn_res.1),
    function(x) {
      message("Processing cluster ", x)
      try(findMarkers(ident1 = x, object = sc_merge, only_pos = TRUE, min_pct = 0.1, logfc_threshold = 0.25, assay = "RNA"))
    }
  )

names(topmarkers) <- unique(sc_merge$RNA_snn_res.1)
write_xlsx(topmarkers, file.path("results", "de" , "topmarkers.xlsx"))
# names(topmarkers) <- sc_merge@misc$cluster_order

# dot plots ---
DefaultAssay(sc_merge) <- "RNA"

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "cellmarkers_seed",
  dot_min = 0.01,
  height = 8,
  width = 15
)

# load Seurat

dplyr::count(sc_merge@meta.data, RNA_snn_res.1) |>
  dplyr::arrange(desc(n))

table(sc_merge$RNA_snn_res.1, sc_merge$diagnosis)
