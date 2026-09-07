sural_study_transfer_result_dir <- function() {
  file.path("results", "targets", "sural_study_transfer")
}

read_sural_study_transfer_query <- function(path) {
  object <- qs::qread(path)
  required_metadata <- c(
    "sample", "ic_cluster", "sex", "age", "center",
    "level0", "level1", "level2", "diagnosis"
  )
  stopifnot(
    inherits(object, "Seurat"),
    identical(names(object@assays), "RNA"),
    identical(SeuratObject::Layers(object[["RNA"]]), "counts"),
    "umap.rpca" %in% names(object@reductions),
    ncol(SeuratObject::Embeddings(object, "umap.rpca")) == 2L,
    all(required_metadata %in% colnames(object[[]])),
    !anyNA(object$sample),
    !anyNA(object$ic_cluster),
    !anyNA(object$diagnosis),
    identical(as.character(object$diagnosis), as.character(object$level2))
  )
  object
}

select_balanced_study_reference_cells <- function(
  object, max_cells_per_patient_group, max_cells_per_cluster_tissue, seed
) {
  stopifnot(
    all(c("patient", "tissue", "cluster") %in% colnames(object[[]])),
    max_cells_per_patient_group >= 1L,
    max_cells_per_cluster_tissue >= 1L
  )
  metadata <- object[[]] |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::transmute(
      cell_id = .data$cell_id,
      patient = as.character(.data$patient),
      tissue = as.character(.data$tissue),
      cluster = as.character(.data$cluster)
    )
  stopifnot(
    !anyNA(metadata$patient),
    !anyNA(metadata$tissue),
    !anyNA(metadata$cluster)
  )

  withr::with_seed(seed, metadata$random_order <- stats::runif(nrow(metadata)))
  selected <- metadata |>
    dplyr::arrange(.data$random_order) |>
    dplyr::group_by(.data$patient, .data$tissue, .data$cluster) |>
    dplyr::slice_head(n = max_cells_per_patient_group) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$random_order) |>
    dplyr::group_by(.data$tissue, .data$cluster) |>
    dplyr::slice_head(n = max_cells_per_cluster_tissue) |>
    dplyr::ungroup()
  stopifnot(
    nrow(selected) > 0L,
    !anyDuplicated(selected$cell_id),
    setequal(unique(metadata$cluster), unique(selected$cluster))
  )
  selected$cell_id
}

prepare_sural_study_transfer_features <- function(
  reference, query, nfeatures
) {
  reference <- Seurat::FindVariableFeatures(
    reference,
    selection.method = "vst",
    nfeatures = max(as.integer(nfeatures * 1.5), nfeatures),
    verbose = FALSE
  )
  features <- SeuratObject::VariableFeatures(reference)
  technical <- grepl(
    "^(MT-|RPL|RPS|HBA[0-9]|HBB|IG[HKL]|TR[ABDG][VJCD])",
    features
  )
  features <- features[!technical & features %in% rownames(query)]
  features <- utils::head(features, nfeatures)
  stopifnot(length(features) >= min(2000L, nfeatures))
  features
}

summarize_sural_study_transfer <- function(predictions) {
  label_summary <- predictions |>
    dplyr::count(
      .data$predicted_cluster, name = "cells", sort = TRUE
    ) |>
    dplyr::mutate(fraction = .data$cells / sum(.data$cells))
  cluster_summary <- predictions |>
    dplyr::group_by(.data$ic_cluster, .data$predicted_cluster) |>
    dplyr::summarise(
      cells = dplyr::n(),
      median_prediction_score = stats::median(.data$prediction_score_max),
      median_prediction_margin = stats::median(.data$prediction_score_margin),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$ic_cluster) |>
    dplyr::mutate(fraction_within_ic_cluster = .data$cells / sum(.data$cells)) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$ic_cluster, dplyr::desc(.data$cells))
  diagnosis_summary <- predictions |>
    dplyr::count(
      .data$diagnosis, .data$predicted_cluster, name = "cells"
    ) |>
    dplyr::group_by(.data$diagnosis) |>
    dplyr::mutate(fraction_within_diagnosis = .data$cells / sum(.data$cells)) |>
    dplyr::ungroup()
  sample_summary <- predictions |>
    dplyr::count(
      .data$sample, .data$diagnosis, .data$predicted_cluster,
      name = "cells"
    ) |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(fraction_within_sample = .data$cells / sum(.data$cells)) |>
    dplyr::ungroup()
  cd8tem3_summary <- predictions |>
    dplyr::group_by(.data$ic_cluster) |>
    dplyr::summarise(
      total_cells = dplyr::n(),
      predicted_cd8tem_3 = sum(.data$predicted_cluster == "CD8TEM_3"),
      raw_cd8tem_3_winner = sum(.data$predicted_cluster_raw == "CD8TEM_3"),
      score_at_least_0_3 = sum(.data$prediction_score_cd8tem_3 >= 0.3),
      score_at_least_0_5 = sum(.data$prediction_score_cd8tem_3 >= 0.5),
      median_cd8tem_3_score = stats::median(.data$prediction_score_cd8tem_3),
      maximum_cd8tem_3_score = max(.data$prediction_score_cd8tem_3),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      predicted_fraction = .data$predicted_cd8tem_3 / .data$total_cells
    ) |>
    dplyr::arrange(dplyr::desc(.data$predicted_cd8tem_3))

  list(
    label_summary = label_summary,
    cluster_summary = cluster_summary,
    diagnosis_summary = diagnosis_summary,
    sample_summary = sample_summary,
    cd8tem3_summary = cd8tem3_summary
  )
}

run_sural_study_label_transfer <- function(
  study, query, max_cells_per_patient_group = 50L,
  max_cells_per_cluster_tissue = 500L, nfeatures = 3000L,
  dimensions = 30L, confidence_threshold = 0.5,
  margin_threshold = 0.1, seed = 20260906L
) {
  stopifnot(
    inherits(study, "Seurat"),
    inherits(query, "Seurat"),
    all(c("patient", "tissue", "cluster") %in% colnames(study[[]])),
    all(c("sample", "ic_cluster", "diagnosis") %in% colnames(query[[]])),
    confidence_threshold > 0,
    confidence_threshold < 1,
    margin_threshold >= 0,
    margin_threshold < 1,
    dimensions >= 2L
  )
  reference_cells <- select_balanced_study_reference_cells(
    study,
    max_cells_per_patient_group,
    max_cells_per_cluster_tissue,
    seed
  )
  reference_metadata <- study[[]][reference_cells, , drop = FALSE]
  reference_counts <- SeuratObject::LayerData(
    study[["RNA"]], layer = "counts"
  )[, reference_cells, drop = FALSE]
  reference <- Seurat::CreateSeuratObject(
    counts = reference_counts,
    assay = "RNA",
    meta.data = reference_metadata
  )
  SeuratObject::DefaultAssay(reference) <- "RNA"
  SeuratObject::DefaultAssay(query) <- "RNA"
  reference <- Seurat::NormalizeData(reference, verbose = FALSE)
  query <- Seurat::NormalizeData(query, verbose = FALSE)
  features <- prepare_sural_study_transfer_features(
    reference, query, as.integer(nfeatures)
  )
  reference <- Seurat::ScaleData(
    reference, features = features, verbose = FALSE
  )
  reference <- Seurat::RunPCA(
    reference,
    features = features,
    npcs = as.integer(dimensions),
    verbose = FALSE,
    seed.use = seed
  )
  anchors <- Seurat::FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = "LogNormalize",
    reduction = "pcaproject",
    reference.reduction = "pca",
    features = features,
    dims = seq_len(dimensions),
    verbose = FALSE
  )
  transferred <- Seurat::TransferData(
    anchorset = anchors,
    refdata = as.character(reference$cluster),
    dims = seq_len(dimensions),
    verbose = FALSE
  )
  score_columns <- setdiff(
    grep("^prediction[.]score[.]", names(transferred), value = TRUE),
    "prediction.score.max"
  )
  stopifnot(
    length(score_columns) == length(unique(reference$cluster)),
    "prediction.score.CD8TEM_3" %in% score_columns,
    identical(rownames(transferred), colnames(query))
  )
  score_matrix <- as.matrix(transferred[, score_columns, drop = FALSE])
  score_order <- t(apply(score_matrix, 1L, order, decreasing = TRUE))
  second_index <- score_order[, 2L]
  second_score <- score_matrix[
    cbind(seq_len(nrow(score_matrix)), second_index)
  ]
  second_label <- sub(
    "^prediction[.]score[.]", "", score_columns[second_index]
  )
  cd8tem3_column <- match("prediction.score.CD8TEM_3", score_columns)
  competitor_score <- apply(
    score_matrix[, -cd8tem3_column, drop = FALSE], 1L, max
  )
  metadata <- query[[]] |>
    tibble::rownames_to_column("cell_id")
  umap <- as.data.frame(SeuratObject::Embeddings(query, "umap.rpca")) |>
    tibble::rownames_to_column("cell_id")
  names(umap)[2:3] <- c("UMAP_1", "UMAP_2")
  predictions <- transferred |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::rename(
      predicted_cluster_raw = "predicted.id",
      prediction_score_max = "prediction.score.max",
      prediction_score_cd8tem_3 = "prediction.score.CD8TEM_3"
    ) |>
    dplyr::mutate(
      predicted_cluster_raw = as.character(.data$predicted_cluster_raw),
      second_predicted_cluster = second_label,
      second_prediction_score = second_score,
      prediction_score_margin =
        .data$prediction_score_max - .data$second_prediction_score,
      cd8tem3_score_margin =
        .data$prediction_score_cd8tem_3 - competitor_score,
      high_confidence =
        .data$prediction_score_max >= confidence_threshold &
        .data$prediction_score_margin >= margin_threshold,
      predicted_cluster = dplyr::if_else(
        .data$high_confidence, .data$predicted_cluster_raw, "unknown"
      ),
      .after = "predicted_cluster_raw"
    ) |>
    dplyr::left_join(metadata, by = "cell_id", relationship = "one-to-one") |>
    dplyr::left_join(umap, by = "cell_id", relationship = "one-to-one")
  stopifnot(
    nrow(predictions) == ncol(query),
    !anyNA(predictions$ic_cluster),
    !anyNA(predictions$UMAP_1),
    !anyNA(predictions$prediction_score_cd8tem_3)
  )
  summaries <- summarize_sural_study_transfer(predictions)
  reference_summary <- reference[[]] |>
    dplyr::count(.data$tissue, .data$cluster, name = "cells") |>
    dplyr::arrange(.data$tissue, .data$cluster)

  list(
    predictions = predictions,
    label_summary = summaries$label_summary,
    cluster_summary = summaries$cluster_summary,
    diagnosis_summary = summaries$diagnosis_summary,
    sample_summary = summaries$sample_summary,
    cd8tem3_summary = summaries$cd8tem3_summary,
    reference_summary = reference_summary,
    features = tibble::tibble(feature = features),
    parameters = tibble::tibble(
      parameter = c(
        "max_cells_per_patient_group", "max_cells_per_cluster_tissue",
        "nfeatures", "dimensions", "confidence_threshold",
        "margin_threshold", "seed"
      ),
      value = as.character(c(
        max_cells_per_patient_group, max_cells_per_cluster_tissue,
        nfeatures, dimensions, confidence_threshold,
        margin_threshold, seed
      ))
    ),
    study_colors = study@misc$cluster_col,
    sural_colors = query@misc$ic_cluster_col
  )
}

write_sural_study_transfer_workbook <- function(tnk_transfer, cd8_transfer) {
  path <- file.path(
    sural_study_transfer_result_dir(), "sural_study_label_transfer.xlsx"
  )
  ensure_parent_dir(path)
  writexl::write_xlsx(
    list(
      tnk_predictions = tnk_transfer$predictions,
      tnk_label_summary = tnk_transfer$label_summary,
      tnk_cluster_summary = tnk_transfer$cluster_summary,
      tnk_diagnosis_summary = tnk_transfer$diagnosis_summary,
      tnk_sample_summary = tnk_transfer$sample_summary,
      tnk_cd8tem3_summary = tnk_transfer$cd8tem3_summary,
      tnk_reference = tnk_transfer$reference_summary,
      tnk_features = tnk_transfer$features,
      cd8_predictions = cd8_transfer$predictions,
      cd8_label_summary = cd8_transfer$label_summary,
      cd8_cluster_summary = cd8_transfer$cluster_summary,
      cd8_diagnosis_summary = cd8_transfer$diagnosis_summary,
      cd8_sample_summary = cd8_transfer$sample_summary,
      cd8_cd8tem3_summary = cd8_transfer$cd8tem3_summary,
      cd8_reference = cd8_transfer$reference_summary,
      cd8_features = cd8_transfer$features,
      parameters = tnk_transfer$parameters
    ),
    path
  )
  path
}

make_sural_cd8tem3_evidence_plot <- function(
  transfer, background, scope_title
) {
  predictions <- transfer$predictions |>
    dplyr::mutate(
      call = dplyr::case_when(
        .data$predicted_cluster == "CD8TEM_3" ~ "High confidence",
        .data$predicted_cluster_raw == "CD8TEM_3" ~
          "Lower confidence winner",
        TRUE ~ NA_character_
      ),
      call = factor(
        .data$call,
        levels = c("Lower confidence winner", "High confidence")
      )
    )
  winners <- dplyr::filter(predictions, !is.na(.data$call))
  base <- ggplot2::ggplot(
    background, ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)
  ) +
    ggplot2::geom_point(color = "grey90", size = 0.1, alpha = 0.3) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_classic()
  score_plot <- base +
    ggplot2::geom_point(
      data = predictions,
      ggplot2::aes(color = .data$prediction_score_cd8tem_3),
      size = 0.1, alpha = 0.5
    ) +
    viridis::scale_color_viridis(option = "magma", direction = -1) +
    ggplot2::labs(title = "Absolute score", color = "CD8TEM_3\nscore")
  margin_limit <- max(abs(predictions$cd8tem3_score_margin))
  margin_plot <- base +
    ggplot2::geom_point(
      data = predictions,
      ggplot2::aes(color = .data$cd8tem3_score_margin),
      size = 0.1, alpha = 0.5
    ) +
    ggplot2::scale_color_gradient2(
      low = "#2166AC", mid = "grey92", high = "#B2182B",
      midpoint = 0, limits = c(-margin_limit, margin_limit)
    ) +
    ggplot2::labs(
      title = "Score minus best competitor",
      color = "CD8TEM_3\nmargin"
    )
  call_plot <- base +
    ggplot2::geom_point(
      data = winners,
      ggplot2::aes(color = .data$call),
      size = 0.8, alpha = 0.95
    ) +
    ggplot2::scale_color_manual(values = c(
      "Lower confidence winner" = "#FDAE61",
      "High confidence" = "#D7191C"
    )) +
    ggplot2::labs(title = "Final call", color = NULL)
  high <- sum(predictions$predicted_cluster == "CD8TEM_3")
  winner_count <- sum(predictions$predicted_cluster_raw == "CD8TEM_3")
  patchwork::wrap_plots(score_plot, margin_plot, call_plot, nrow = 1L) +
    patchwork::plot_annotation(
      title = paste("CD8TEM_3 evidence:", scope_title),
      subtitle = paste(
        high, "high-confidence calls among", winner_count,
        "raw winners; positive margins indicate CD8TEM_3 is the top label"
      )
    )
}

make_sural_transfer_heatmap <- function(transfer, title) {
  data <- transfer$cluster_summary |>
    dplyr::mutate(
      label = dplyr::if_else(
        .data$fraction_within_ic_cluster >= 0.05,
        scales::percent(.data$fraction_within_ic_cluster, accuracy = 1),
        ""
      )
    )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$predicted_cluster,
      y = .data$ic_cluster,
      fill = .data$fraction_within_ic_cluster
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 3) +
    viridis::scale_fill_viridis(
      option = "magma", labels = scales::label_percent()
    ) +
    ggplot2::labs(
      title = title,
      x = "Predicted study T/NK cluster",
      y = "Existing sural cluster",
      fill = "Within-cluster\nfraction"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

make_sural_transfer_label_umap <- function(
  transfer, background, title
) {
  label_colors <- c(transfer$study_colors, unknown = "#BDBDBD")
  ggplot2::ggplot(
    background, ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)
  ) +
    ggplot2::geom_point(color = "grey92", size = 0.1, alpha = 0.25) +
    ggplot2::geom_point(
      data = transfer$predictions,
      ggplot2::aes(color = .data$predicted_cluster),
      size = 0.2, alpha = 0.75
    ) +
    ggplot2::scale_color_manual(values = label_colors, drop = TRUE) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title, x = "UMAP 1", y = "UMAP 2",
      color = "Predicted cluster"
    ) +
    ggplot2::theme_classic()
}

write_sural_study_transfer_plots <- function(
  tnk_transfer, cd8_transfer, query
) {
  root <- sural_study_transfer_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  background <- as.data.frame(
    SeuratObject::Embeddings(query, "umap.rpca")
  )
  names(background) <- c("UMAP_1", "UMAP_2")
  paths <- file.path(root, c(
    "all_tnk_cd8tem3_evidence_umap.pdf",
    "cd8_nkcd8_cd8tem3_evidence_umap.pdf",
    "all_tnk_transfer_heatmap.pdf",
    "cd8_nkcd8_transfer_heatmap.pdf",
    "all_tnk_predicted_labels_umap.pdf",
    "cd8_nkcd8_predicted_labels_umap.pdf"
  ))
  plots <- list(
    make_sural_cd8tem3_evidence_plot(
      tnk_transfer, background, "all T/NK cells"
    ),
    make_sural_cd8tem3_evidence_plot(
      cd8_transfer, background, "CD8 and NK_CD8 cells"
    ),
    make_sural_transfer_heatmap(
      tnk_transfer, "All sural T/NK cells mapped to study T/NK clusters"
    ),
    make_sural_transfer_heatmap(
      cd8_transfer, "Sural CD8/NK_CD8 cells mapped to study T/NK clusters"
    ),
    make_sural_transfer_label_umap(
      tnk_transfer, background, "Predicted labels: all sural T/NK cells"
    ),
    make_sural_transfer_label_umap(
      cd8_transfer, background, "Predicted labels: sural CD8/NK_CD8 cells"
    )
  )
  widths <- c(15, 15, 11, 9, 10, 10)
  heights <- c(5.5, 5.5, 6.5, 4.5, 7, 7)
  for (index in seq_along(paths)) {
    ggplot2::ggsave(
      paths[[index]], plots[[index]],
      width = widths[[index]], height = heights[[index]]
    )
  }
  paths
}
