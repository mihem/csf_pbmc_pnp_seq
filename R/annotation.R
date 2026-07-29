cluster_integrated_object <- function(object, config) {
  configure_runtime(config)
  object <- Seurat::FindNeighbors(
    object,
    reduction = "stacas.ss.all",
    dims = seq_len(as.integer(config$annotation$dimensions)),
    assay = "RNA",
    graph.name = c("RNA_nn", "RNA_snn")
  )
  for (resolution in unlist(config$annotation$resolutions)) {
    object <- Seurat::FindClusters(
      object,
      resolution = resolution,
      random.seed = as.integer(config$annotation$cluster_seed)
    )
  }
  object
}

read_manual_overrides <- function(config, manual_files) {
  stopifnot(length(manual_files) == 7L)
  outlier_files <- file.path(
    "lookup",
    paste0(c("cl9", "cl12", "cl23", "cl24"), "_outliers.csv")
  )
  outliers <- purrr::map_dfr(outlier_files, function(path) {
    readr::read_csv(path, show_col_types = FALSE) |>
      dplyr::transmute(
        cell_barcode,
        cluster = paste0(RNA_snn_res.1, "_outliers")
      )
  })
  nk <- readr::read_csv(config$paths$nk_overrides, show_col_types = FALSE) |>
    dplyr::transmute(cell_barcode, cluster = "23_nk_cells")
  stopifnot(!anyDuplicated(outliers$cell_barcode), !anyDuplicated(nk$cell_barcode))
  list(outliers = outliers, nk = nk)
}

apply_annotations <- function(object, overrides, config) {
  metadata <- object@meta.data |>
    tibble::rownames_to_column("cell_barcode") |>
    dplyr::left_join(
      dplyr::rename(overrides$outliers, cluster_outlier = cluster),
      by = "cell_barcode"
    ) |>
    dplyr::mutate(
      cluster = dplyr::coalesce(cluster_outlier, as.character(RNA_snn_res.1))
    ) |>
    dplyr::select(-cluster_outlier) |>
    dplyr::left_join(
      dplyr::rename(overrides$nk, cluster_nk = cluster),
      by = "cell_barcode"
    ) |>
    dplyr::mutate(cluster = dplyr::coalesce(cluster_nk, cluster)) |>
    dplyr::select(-cluster_nk)

  annotations <- readxl::read_xlsx(config$paths$annotations) |>
    dplyr::transmute(cluster = as.character(number), final, cluster_order)
  stopifnot(!anyDuplicated(annotations$cluster))
  metadata <- metadata |>
    dplyr::left_join(
      dplyr::select(annotations, cluster, final),
      by = "cluster"
    )
  stopifnot(!anyNA(metadata$final))
  metadata_out <- as.data.frame(
    dplyr::select(metadata, -cluster, -cell_barcode)
  )
  rownames(metadata_out) <- metadata$cell_barcode
  object@meta.data <- metadata_out

  cluster_order <- annotations$cluster_order[!is.na(annotations$cluster_order)]
  object$cluster <- factor(metadata$final, levels = cluster_order)
  Seurat::Idents(object) <- object$cluster

  cluster_colors <- c(
    "#FFC01C", "#C6E635", "#16CDC5", "#FBA7DD", "#ED4F69",
    "#F8A10D", "#A3CB32", "#1686A3", "#D881F9", "#B02A6C",
    "#E45300", "#22A145", "#1655DD", "#9B81FB", "#813270",
    "#F62E38", "#16696D", "#4F3B7F", "#5858BC", "#773260",
    "#E2DDC2", "#22FE2A", "#F800D1", "#86600D", "#F29488",
    "#D8D6FF", "#FB009F"
  )
  stopifnot(length(cluster_order) == length(cluster_colors))
  object@misc$cluster_order <- unname(cluster_order)
  object@misc$cluster_col <- stats::setNames(
    cluster_colors,
    cluster_order
  )
  object@misc$diagnosis_order <- c(
    "CTRL", "CIAP", "GBS", "CIDP", "MAG", "MFS", "PNC", "CAN", "PPN"
  )
  object@misc$diagnosis_col <- stats::setNames(
    pals::cols25(length(object@misc$diagnosis_order)),
    object@misc$diagnosis_order
  )
  object@misc$tissue_diagnosis_order <- c(
    paste0("CSF_", object@misc$diagnosis_order),
    paste0("PBMC_", object@misc$diagnosis_order)
  )
  object@misc$tissue_diagnosis_col <- stats::setNames(
    pals::cols25(length(object@misc$tissue_diagnosis_order)),
    object@misc$tissue_diagnosis_order
  )
  object$tissue_diagnosis <- factor(
    paste0(object$tissue, "_", object$diagnosis),
    levels = object@misc$tissue_diagnosis_order
  )
  object$tissue_group <- paste0(object$tissue, "_", object$group)
  Seurat::DefaultAssay(object) <- "RNA"
  object
}

export_cerebro_checkpoint <- function(object, path) {
  ensure_parent_dir(path)
  cerebro <- Seurat::DietSeurat(object, dimreducs = "umap.stacas.ss.all")
  cerebro[["RNA"]]$scale.data <- NULL
  cerebro[["RNA"]]$counts <- NULL
  cerebro@meta.data <- cerebro@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(
      barcode, sample, tissue, diagnosis, nCount_RNA,
      nFeature_RNA, percent_mt, RNA_snn_res.1
    ) |>
    tibble::column_to_rownames("barcode")
  cerebroAppLite::exportFromSeurat(
    object = cerebro,
    file = path,
    experiment_name = "PNP",
    groups = "RNA_snn_res.1",
    organism = "hg",
    nUMI = "nCount_RNA",
    nGene = "nFeature_RNA",
    use_delayed_array = FALSE
  )
  path
}

write_annotation_umap <- function(object, path) {
  ensure_parent_dir(path)
  plot <- Seurat::DimPlot(
    object,
    reduction = "umap.stacas.ss.all",
    pt.size = 0.1,
    raster = FALSE,
    alpha = 0.1,
    group.by = "cluster",
    cols = object@misc$cluster_col,
    label = TRUE
  ) +
    scMisc::theme_rect() +
    Seurat::NoLegend() +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    ggplot2::ggtitle("")
  ggplot2::ggsave(path, plot = plot, width = 6, height = 6)
  path
}

write_annotation_dotplot <- function(
  object, marker_file, marker_column, path, width, height
) {
  markers <- utils::read.csv(marker_file, check.names = FALSE)[[marker_column]]
  markers <- unique(markers[!is.na(markers) & nzchar(markers)])
  stopifnot(length(markers) > 0L, all(markers %in% rownames(object)))
  Seurat::DefaultAssay(object) <- "RNA"
  Seurat::Idents(object) <- object$cluster
  plot <- Seurat::DotPlot(
    object,
    features = markers,
    dot.scale = 10,
    scale.by = "size",
    dot.min = 0.01,
    scale = TRUE
  ) +
    viridis::scale_color_viridis(option = "viridis") +
    ggplot2::scale_size(range = c(0, 10)) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1, face = "italic"
      )
    ) +
    ggplot2::labs(x = NULL, y = NULL)
  ensure_parent_dir(path)
  ggplot2::ggsave(path, plot, width = width, height = height)
  path
}
