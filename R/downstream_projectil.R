projectil_result_dir <- function() {
  file.path("results", "targets", "projectil")
}

assess_projectil_readiness <- function(object, reference_path) {
  package_available <- requireNamespace("ProjecTILs", quietly = TRUE)
  reference_exists <- file.exists(reference_path)
  query_valid <- inherits(object, "Seurat") &&
    "cluster" %in% colnames(object[[]]) &&
    "diagnosis" %in% colnames(object[[]])
  query_cells <- if (query_valid) {
    sum(as.character(object$cluster) == "CD8TEM_3")
  } else {
    0L
  }
  reference_valid <- FALSE
  shared_genes <- 0L
  if (reference_exists) {
    reference <- qs::qread(reference_path)
    reference_valid <- inherits(reference, "Seurat") &&
      all(c("RNA", "integrated") %in% names(reference@assays)) &&
      all(c("pca", "umap") %in% names(reference@reductions))
    if (query_valid && reference_valid) {
      shared_genes <- length(intersect(rownames(object), rownames(reference)))
    }
  }
  ready <- package_available && reference_valid && query_cells > 0L && shared_genes > 0L

  tibble::tibble(
    status = if (ready) "ready_for_manual_specification" else "blocked",
    package_available = package_available,
    reference_exists = reference_exists,
    reference_valid = reference_valid,
    query_cd8tem_3_cells = query_cells,
    shared_genes = shared_genes,
    reason = if (ready) {
      paste(
        "Inputs are runnable, but no projection is emitted until diagnosis",
        "contrasts and plotting semantics replace the broken exploratory script."
      )
    } else {
      "One or more package, reference, query-cell, or feature-overlap requirements failed."
    }
  )
}

write_projectil_readiness <- function(status, path) {
  stopifnot(
    nrow(status) == 1L,
    status$status %in% c("ready_for_manual_specification", "blocked")
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(status, path)
  path
}

projectil_terekhova_colors <- function() {
  c(
    "Naive" = "#854c90", "Naive-IFN" = "#a3a5aa", "Trm" = "#bfa0d1",
    "Tmem KLRC2+" = "#4fa6a2", "Tem GZMK+" = "#406f9d",
    "Tem GZMB+" = "#e89e34", "NKT-like" = "#99ad6e",
    "Temra" = "#dd587d", "Proliferative" = "#fc8f7e",
    "HLA-DR+" = "#d22b90", "Tcm CCR4+" = "#eb8123",
    "Tcm CCR4-" = "#99ad6e"
  )
}

run_projectil_cd8tem3 <- function(object, reference_path, ncores = 2L) {
  suppressPackageStartupMessages(library("ProjecTILs", character.only = TRUE))
  reference <- qs::qread(reference_path)
  cells <- colnames(object)[as.character(object$cluster) == "CD8TEM_3"]
  stopifnot(length(cells) > 0L)
  query <- object[, cells]
  projected <- ProjecTILs::make.projection(
    query,
    ref = reference,
    ncores = as.integer(ncores),
    filter.cells = FALSE
  )
  list(reference = reference, projected = projected)
}

write_projectil_plots <- function(projection) {
  root <- projectil_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  diagnoses <- c(all = NA_character_, CIAP = "CIAP", CIDP = "CIDP", GBS = "GBS")
  colors <- projectil_terekhova_colors()
  paths <- vapply(names(diagnoses), function(name) {
    projected <- projection$projected
    diagnosis <- diagnoses[[name]]
    if (!is.na(diagnosis)) {
      cells <- colnames(projected)[as.character(projected$diagnosis) == diagnosis]
      stopifnot(length(cells) > 0L)
      projected <- projected[, cells]
    }
    plot <- ProjecTILs::plot.projection(
      projection$reference,
      projected,
      linesize = 0.5,
      pointsize = 1,
      ref.size = 0.1,
      cols = colors
    ) +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
          color = "black", linewidth = 1, fill = NA
        ),
        aspect.ratio = 0.7
      ) +
      ggplot2::labs(
        x = "UMAP1", y = "UMAP2",
        title = if (name == "all") "CD8TEM_3" else paste("CD8TEM_3", diagnosis)
      )
    path <- file.path(
      root,
      paste0("cd8tem_3_", tolower(name), "_terekhova_projectil.pdf")
    )
    ggplot2::ggsave(path, plot, width = 7, height = 4)
    path
  }, character(1))
  unname(paths)
}

sural_projectil_result_dir <- function() {
  file.path(projectil_result_dir(), "sural")
}

run_sural_projectil_cd8tem3 <- function(object, reference_path, ncores = 2L) {
  suppressPackageStartupMessages(library("ProjecTILs", character.only = TRUE))
  reference_input <- qs::qread(reference_path)
  variable_features <- SeuratObject::VariableFeatures(
    reference_input, assay = "RNA"
  )
  cells <- colnames(object)[as.character(object$cluster) == "CD8TEM_3"]
  stopifnot(
    inherits(reference_input, "Seurat"),
    identical(names(reference_input@assays), "RNA"),
    identical(SeuratObject::Layers(reference_input[["RNA"]]), "data"),
    length(variable_features) > 0L,
    all(c("sample", "ic_cluster") %in% colnames(reference_input[[]])),
    "umap.rpca" %in% names(reference_input@reductions),
    ncol(SeuratObject::Embeddings(reference_input, "umap.rpca")) == 2L,
    !anyNA(reference_input$ic_cluster),
    setequal(
      names(reference_input@misc$ic_cluster_col),
      levels(factor(reference_input$ic_cluster))
    ),
    length(intersect(variable_features, rownames(object))) /
      length(variable_features) >= 0.8,
    length(cells) > 0L
  )
  reference <- ProjecTILs::make.reference(
    reference_input,
    assay = "RNA",
    atlas.name = "sural_immune",
    annotation.column = "ic_cluster",
    recalculate.umap = FALSE,
    dimred = "umap.rpca",
    ndim = 30L,
    color.palette = reference_input@misc$ic_cluster_col,
    store.markers = FALSE,
    seed = 42
  )
  reference[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = SeuratObject::Embeddings(reference, "umap.rpca"),
    key = "UMAP_",
    assay = "integrated"
  )
  reference[["umap.rpca"]] <- NULL
  projected <- ProjecTILs::make.projection(
    object[, cells],
    ref = reference,
    ncores = as.integer(ncores),
    filter.cells = FALSE
  )
  predicted <- ProjecTILs::cellstate.predict(
    reference,
    query = projected,
    reduction = "pca",
    ndim = 30L,
    k = 5L,
    min.confidence = 0.2
  )
  list(reference = reference, projected = predicted)
}

write_sural_projectil_outputs <- function(projection) {
  root <- sural_projectil_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  reference <- projection$reference
  projected <- projection$projected
  stopifnot(
    all(c("functional.cluster", "functional.cluster.conf") %in%
      colnames(projected[[]])),
    !anyNA(projected$functional.cluster)
  )

  prediction_path <- file.path(root, "cd8tem_3_cell_predictions.csv")
  summary_path <- file.path(root, "cd8tem_3_prediction_summary.csv")
  plot_path <- file.path(root, "cd8tem_3_all_sural_projectil.pdf")
  predictions <- projected[[]] |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::select(
      "cell_id", "tissue", "diagnosis", "sample", "patient",
      "functional.cluster", "functional.cluster.conf"
    ) |>
    dplyr::rename(
      predicted_sural_cluster = "functional.cluster",
      prediction_confidence = "functional.cluster.conf"
    )
  readr::write_csv(predictions, prediction_path)
  summary <- predictions |>
    dplyr::count(.data$predicted_sural_cluster, name = "cells") |>
    dplyr::mutate(fraction = .data$cells / sum(.data$cells)) |>
    dplyr::arrange(dplyr::desc(.data$cells))
  readr::write_csv(summary, summary_path)

  plot <- ProjecTILs::plot.projection(
    reference,
    projected,
    linesize = 0.15,
    pointsize = 0.35,
    ref.size = 0.1,
    cols = reference@misc$atlas.palette
  ) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "black", linewidth = 1, fill = NA
      )
    ) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2", title = "CD8TEM_3")
  ggplot2::ggsave(plot_path, plot, width = 8, height = 6)

  c(plot_path, prediction_path, summary_path)
}
