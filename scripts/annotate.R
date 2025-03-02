##################################################
# annotate data
# requires running integrate.R first
##################################################

# load libraries
library(tidyverse)
library(Seurat)
library(qs)
library(scMisc)
library(presto)
library(clustree)
library(readxl)
library(writexl)
library(cerebroAppLite)
library(Polychrome)

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
  graph.name = c("RNA_nn", "RNA_snn")
)

for (res in seq(from = 0.2, to = 1.2, by = 0.1)) {
  sc_merge <- FindClusters(sc_merge, resolution = res)
}

sc_merge <- FindClusters(sc_merge, resolution = 1)

qs::qsave(sc_merge, file.path("objects", "sc_merge.qs"))

# plot clustering ---
clustree <- clustree(sc_merge, prefix = "RNA_snn_res.")

ggsave(
  plot = clustree,
  file.path("results", "umap", "stacas_ss_all_clustree.png"),
  width = 15,
  height = 15
)

my_cols_100 <- unname(Polychrome::createPalette(100, pals::cols25()))
resolutions <- paste0("RNA_snn_res.", seq(from = 0.2, to = 1.2, by = 0.1))

umap_list <- list()
for (res in resolutions) {
  umap_list[[res]] <-
    DimPlot(
      sc_merge,
      reduction = "umap.stacas.ss.all",
      pt.size = .1,
      raster = FALSE,
      alpha = 0.1,
      group.by = res,
      cols = my_cols_100,
      label = TRUE
    ) +
    theme_rect() +
    NoLegend()
}

umap_list <- patchwork::wrap_plots(umap_list, ncol = 2)
ggsave(
  plot = umap_list,
  file.path("results", "umap", "stacas_ss_all_resolutions.png"),
  width = 16,
  height = 40
)

#  use resolution 1

# umap ---
umap_cluster <-
  DimPlot(
    sc_merge,
    reduction = "umap.stacas.ss.all",
    pt.size = .1,
    raster = FALSE,
    alpha = 0.1,
    group.by = "RNA_snn_res.1",
    cols = my_cols_100,
    label = TRUE
  ) +
  theme_rect() +
  NoLegend()

ggsave(
  plot = umap_cluster,
  file.path("results", "umap", "stacas_ss_all_res1_cluster.png"),
  width = 10,
  height = 10
)

# split by cluster
umap_split_cluster <-
  DimPlot(
    sc_merge,
    reduction = "umap.stacas.ss.all",
    pt.size = .1,
    raster = FALSE,
    alpha = 0.1,
    split.by = "RNA_snn_res.1",
    cols = my_cols_100,
    label = TRUE
  ) +
  theme_rect() +
  NoLegend()

ggsave(
  plot = umap_split_cluster,
  file.path("results", "umap", "stacas_ss_all_split_res1.png"),
  width = 50,
  height = 5,
  limitsize = FALSE
)

# manually select misclustered cell with cerebro ----
sc_cerebro <-
  DietSeurat(
    sc_merge,
    dimreducs = c("umap.stacas.ss.all")
  )

# somehow the scale.data is still present, only data required for gene expression
sc_cerebro[["RNA"]]$scale.data <- NULL
sc_cerebro[["RNA"]]$counts <- NULL

sc_cerebro@meta.data <-
  sc_cerebro@meta.data |>
  tibble::rownames_to_column("barcode") |>
  select(
    barcode,
    sample,
    tissue,
    diagnosis,
    nCount_RNA,
    nFeature_RNA,
    percent_mt,
    RNA_snn_res.1
  ) |>
  tibble::column_to_rownames(var = "barcode")

cerebroAppLite::exportFromSeurat(
  object = sc_cerebro,
  file = file.path("objects", "sc_cerebro.crb"),
  experiment_name = "PNP",
  groups = "RNA_snn_res.1",
  organism = "hg",
  nUMI = "nCount_RNA",
  nGene = "nFeature_RNA",
  use_delayed_array = FALSE
)

cerebroAppLite::launchCerebro(
  crb_file_to_load = file.path("objects", "sc_cerebro.crb"),
)

# load manually selected cells ---
outliers_clusters <- c("cl9", "cl12", "cl23", "cl24")
outliers <- map_dfr(
  outliers_clusters,
  function(x) {
    read_csv(file.path("lookup", paste0(x, "_outliers.csv"))) |>
      mutate(cluster = paste0(RNA_snn_res.1, "_outliers")) |>
      select(cell_barcode, cluster)
  }
)

sc_merge@meta.data <-
  sc_merge@meta.data |>
  tibble::rownames_to_column("cell_barcode") |>
  dplyr::left_join(outliers) |>
  dplyr::mutate(cluster = dplyr::coalesce(cluster, RNA_snn_res.1)) |>
  tibble::column_to_rownames(var = "cell_barcode")

# annotate clusters ----
Idents(sc_merge) <- factor(
  sc_merge$cluster,
  levels = sort(unique(sc_merge$cluster))
)

annotations <- readxl::read_xlsx(file.path("lookup", "annotations.xlsx")) |>
  mutate(cluster = as.character(number)) |>
  arrange(cluster)

names(annotations$final) <- levels(sc_merge)
sc_merge <- RenameIdents(sc_merge, annotations$final)

# save custom cluster and condition order and colors in seurat object
cluster_order <-
  read_xlsx(file.path("lookup", "annotations.xlsx")) |>
  pull(cluster_order)

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


sc_merge@misc$cluster_order <- cluster_order[!is.na(cluster_order)]
sc_merge@misc$cluster_col <- setNames(
  unname(Polychrome::createPalette(
    length(sc_merge@misc$cluster_order),
    colors_dutch
  )),
  sc_merge@misc$cluster_order
)

sc_merge@misc$diagnosis_order <- c(
  "CTRL",
  "CIAP",
  "GBS",
  "CIDP",
  "MAG",
  "MFS",
  "PNC",
  "CAN",
  "PPN"
)

sc_merge@misc$diagnosis_col <- setNames(
  pals::cols25(length(sc_merge@misc$diagnosis_order)),
  sc_merge@misc$diagnosis_order
)

sc_merge@misc$tissue_diagnosis_order <- c(
  paste0("CSF_", sc_merge@misc$diagnosis_order),
  paste0("PBMC_", sc_merge@misc$diagnosis_order)
)

sc_merge@misc$tissue_diagnosis_col <- setNames(
  pals::cols25(length(sc_merge@misc$tissue_diagnosis_order)),
  sc_merge@misc$tissue_diagnosis_order
)

sc_merge$tissue_diagnosis <- factor(paste0(sc_merge$tissue, "_", sc_merge$diagnosis), levels = sc_merge@misc$tissue_diagnosis_order)
sc_merge$tissue_group <- paste0(sc_merge$tissue, "_", sc_merge$group)

#sanity check
setdiff(sc_merge@misc$cluster_order, levels(sc_merge))
setdiff(levels(sc_merge), sc_merge@misc$cluster_order)
any(duplicated(sc_merge@misc$cluster_order))

# save ordered annotations in meta.data
sc_merge$cluster <- factor(
  Idents(sc_merge),
  levels = sc_merge@misc$cluster_order
)
Idents(sc_merge) <- sc_merge$cluster
DefaultAssay(sc_merge) <- "RNA"

# final annotated UMAP plot by cluster
umap <-
  DimPlot(
    sc_merge,
    reduction = "umap.stacas.ss.all",
    pt.size = .1,
    raster = FALSE,
    alpha = 0.1,
    group.by = "cluster",
    cols = sc_merge@misc$cluster_col,
    label = TRUE
  ) +
  theme_rect() +
  NoLegend() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("")

ggsave(
  plot = umap,
  file.path("results", "umap", "umap_annotated.png"),
  width = 10,
  height = 10
)

qs::qsave(sc_merge, file.path("objects", "sc_merge.qs"))


# find markers helper function
findMarkers <- function(
  ident1,
  ident2 = NULL,
  object,
  only_pos,
  min_pct,
  logfc_threshold,
  assay = assay
) {
  result <- Seurat::FindMarkers(
    object,
    ident.1 = ident1,
    ident.2 = ident2,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    only.pos = only_pos,
    assay = assay
  ) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::relocate(gene, avg_log2FC, p_val, p_val_adj) |>
    dplyr::arrange(desc(avg_log2FC))
  return(result)
}

topmarkers <-
  lapply(
    sc_merge@misc$cluster_order,
    function(x) {
      message("Processing cluster ", x)
      try(findMarkers(
        ident1 = x,
        object = sc_merge,
        only_pos = TRUE,
        min_pct = 0.1,
        logfc_threshold = 0.25,
        assay = "RNA"
      ))
    }
  )

names(topmarkers) <- sc_merge@misc$cluster_order
write_xlsx(topmarkers, file.path("results", "de", "topmarkers.xlsx"))

# dot plots ---
DefaultAssay(sc_merge) <- "RNA"

Idents(sc_merge) <- sc_merge$cluster

dotPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "cellmarkers_seed",
  dot_min = 0.01,
  height = 8,
  width = 12
)

# feature plots ---
scMisc::fPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "cellmarkers_seed",
  reduction = "umap.stacas.ss.all",
  order = FALSE,
  width = 24
)

# feature plots ---
scMisc::fPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "bc_gerd",
  reduction = "umap.stacas.ss.all",
  order = FALSE,
  width = 32,
  height = 32
)

# feature plots ---
scMisc::fPlot(
  path = file.path("lookup", "markers.csv"),
  object = sc_merge,
  par = "dying_cells",
  reduction = "umap.stacas.ss.all",
  order = FALSE,
  width = 40,
  height = 20
)
