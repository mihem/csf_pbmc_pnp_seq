cerebro_result_dir <- function() {
  file.path("results", "targets", "cerebro")
}

prepare_cerebro_export <- function(object) {
  metadata <- c(
    "cluster", "diagnosis", "tissue", "nCount_RNA", "nFeature_RNA",
    "percent_mt", "sample", "patient"
  )
  stopifnot(
    inherits(object, "Seurat"),
    "RNA" %in% names(object@assays),
    "umap.stacas.ss.all" %in% names(object@reductions),
    all(metadata %in% colnames(object[[]]))
  )

  export <- Seurat::DietSeurat(
    object,
    layers = "data",
    assays = "RNA",
    dimreducs = "umap.stacas.ss.all"
  )
  export[["RNA"]]$scale.data <- NULL
  export[["RNA"]]$counts <- NULL
  export@meta.data <- export@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(barcode, dplyr::all_of(metadata)) |>
    tibble::column_to_rownames("barcode")
  stopifnot(
    identical(colnames(export), rownames(export[[]])),
    identical(SeuratObject::Layers(export[["RNA"]]), "data")
  )
  export
}

add_cerebro_marker_genes <- function(object) {
  object <- cerebroAppLite::getMarkerGenes(
    object,
    # Cell-surface annotation uses a live Ensembl query and is not reproducible.
    organism = "offline",
    groups = "cluster",
    only_pos = TRUE,
    assay = "RNA",
    min_pct = 0.1,
    thresh_logFC = 0.25,
    thresh_p_val = 0.05,
    verbose = FALSE
  )
  markers <- object@misc$marker_genes$cerebro_seurat$cluster
  stopifnot(is.data.frame(markers), nrow(markers) > 0L, "avg_log2FC" %in% names(markers))
  object@misc$marker_genes$cerebro_seurat$cluster <- dplyr::arrange(
    markers, dplyr::desc(avg_log2FC)
  )
  object
}

write_cerebro_object <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  qs::qsave(object, path)
  stopifnot(file.exists(path), file.info(path)$size > 0L)
  path
}

export_cerebro_artifact <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  cerebroAppLite::exportFromSeurat(
    object = object,
    file = path,
    experiment_name = "IN-seq",
    groups = c("cluster", "sample", "tissue", "diagnosis"),
    organism = "hg",
    nUMI = "nCount_RNA",
    nGene = "nFeature_RNA",
    use_delayed_array = FALSE
  )
  stopifnot(file.exists(path), file.info(path)$size > 0L)
  path
}

write_cerebro_export_status <- function(object_file, cerebro_file, path) {
  stopifnot(all(file.exists(c(object_file, cerebro_file))))
  status <- tibble::tibble(
    status = "exported",
    artifact = c("seurat_checkpoint", "cerebro_export"),
    file = c(object_file, cerebro_file),
    bytes = as.numeric(file.info(c(object_file, cerebro_file))$size),
    next_step = c(
      "Retained reproducible checkpoint; no manual action required",
      "Open manually with cerebroAppLite::launchCerebro(crb_file_to_load = file)"
    )
  )
  stopifnot(all(status$bytes > 0L))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(status, path)
  path
}
