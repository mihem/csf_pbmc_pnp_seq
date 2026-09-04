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
  assignment_path <- file.path(root, "cd8tem_3_sural_cluster_assignments.pdf")
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
  tnk_clusters <- c("CD4", "Treg", "MAIT", "CD4_CD8", "CD8", "NK_CD8", "NK")
  reference_frequencies <- reference[[]] |>
    dplyr::filter(.data$functional.cluster %in% tnk_clusters) |>
    dplyr::count(.data$functional.cluster, name = "reference_cells") |>
    dplyr::mutate(
      reference_fraction = .data$reference_cells / sum(.data$reference_cells)
    )
  summary <- predictions |>
    dplyr::count(.data$predicted_sural_cluster, name = "cells") |>
    dplyr::mutate(fraction = .data$cells / sum(.data$cells)) |>
    dplyr::left_join(
      reference_frequencies,
      by = c("predicted_sural_cluster" = "functional.cluster")
    ) |>
    dplyr::mutate(
      enrichment = .data$fraction / .data$reference_fraction
    ) |>
    dplyr::arrange(dplyr::desc(.data$enrichment))
  stopifnot(!anyNA(summary$enrichment))
  readr::write_csv(summary, summary_path)

  assignment_plot <- summary |>
    dplyr::mutate(
      predicted_sural_cluster = factor(
        .data$predicted_sural_cluster,
        levels = rev(.data$predicted_sural_cluster)
      ),
      label = sprintf("%.2fx (%.1f%%)", .data$enrichment, 100 * .data$fraction)
    ) |>
    ggplot2::ggplot(ggplot2::aes(
      x = .data$enrichment,
      y = .data$predicted_sural_cluster,
      fill = .data$predicted_sural_cluster
    )) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label), hjust = -0.1, size = 3.5
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.25))
    ) +
    ggplot2::scale_fill_manual(
      values = reference@misc$atlas.palette, guide = "none"
    ) +
    ggplot2::labs(
      x = "Observed / expected assignment", y = NULL,
      title = "Normalized sural cluster assignment",
      subtitle = "Expected from T/NK reference cluster abundance"
    ) +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(assignment_path, assignment_plot, width = 6, height = 4)

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

  c(plot_path, assignment_path, prediction_path, summary_path)
}
