deg_comparisons <- function() {
  tibble::tribble(
    ~comparison, ~condition_1, ~condition_2, ~tissue,
    "cidp_ctrl_csf", "CIDP", "CTRL", "CSF",
    "gbs_ctrl_csf", "GBS", "CTRL", "CSF",
    "cidp_ctrl_pbmc", "CIDP", "CTRL", "PBMC",
    "gbs_ctrl_pbmc", "GBS", "CTRL", "PBMC"
  )
}

subset_deg_object <- function(object, condition_1, condition_2, tissue) {
  required <- c("cluster", "diagnosis", "patient", "tissue")
  stopifnot(
    inherits(object, "Seurat"),
    all(required %in% colnames(object[[]])),
    condition_1 != condition_2,
    length(tissue) == 1L
  )

  keep <- object$diagnosis %in% c(condition_1, condition_2) &
    object$tissue == tissue
  stopifnot(any(keep), !anyNA(keep))
  result <- subset(object, cells = colnames(object)[keep])
  result$diagnosis <- droplevels(result$diagnosis)
  Seurat::Idents(result) <- result$cluster

  stopifnot(
    identical(unique(as.character(result$tissue)), tissue),
    setequal(unique(as.character(result$diagnosis)), c(condition_1, condition_2)),
    !anyNA(result$patient),
    !anyNA(result$cluster)
  )
  result
}

make_deg_cluster_pseudobulk <- function(object) {
  result <- Libra::to_pseudobulk(
    object,
    cell_type_col = "cluster",
    label_col = "diagnosis",
    replicate_col = "patient"
  )
  stopifnot(
    is.list(result),
    length(result) > 0L,
    !is.null(names(result)),
    !anyDuplicated(names(result)),
    all(vapply(result, ncol, integer(1)) > 0L)
  )
  result
}

deg_cluster_names <- function(pseudobulk, cluster_order) {
  cluster_order[cluster_order %in% names(pseudobulk)]
}

deg_cluster_eligibility <- function(object, pseudobulk) {
  clusters <- levels(object$cluster)
  if (is.null(clusters)) {
    clusters <- sort(unique(as.character(object$cluster)))
  }
  tibble::tibble(
    cluster = clusters,
    modeled = clusters %in% names(pseudobulk),
    status = ifelse(
      modeled,
      "modeled",
      "excluded_by_Libra_pseudobulk_eligibility"
    )
  )
}

deg_model_metadata <- function(counts, lookup) {
  stopifnot(
    is.matrix(counts) || inherits(counts, "Matrix"),
    !is.null(colnames(counts)),
    all(c("patient", "diagnosis", "sex", "age") %in% colnames(lookup)),
    !anyDuplicated(lookup$patient)
  )
  metadata <- tibble::tibble(
    patient = sub(":.+", "", colnames(counts))
  ) |>
    dplyr::left_join(lookup, by = "patient")
  stopifnot(
    nrow(metadata) == ncol(counts),
    !anyNA(metadata$patient),
    !anyNA(metadata$diagnosis),
    !anyNA(metadata$sex),
    !anyNA(metadata$age)
  )
  metadata
}

run_deg_limma <- function(counts, lookup, condition_1, condition_2) {
  dge <- edgeR::DGEList(counts = counts, group = colnames(counts))
  keep <- rowSums(edgeR::cpm(dge) > 1) > 2
  stopifnot(any(keep))
  dge <- edgeR::calcNormFactors(dge[keep, ], method = "TMM")

  metadata <- deg_model_metadata(dge$counts, lookup)
  design <- stats::model.matrix(
    ~ 0 + diagnosis + sex + age,
    data = metadata
  )
  contrast_name <- paste0(
    "diagnosis", condition_1, "-diagnosis", condition_2
  )
  contrast <- limma::makeContrasts(contrasts = contrast_name, levels = design)

  fit <- limma::voomWithQualityWeights(dge, design, plot = FALSE) |>
    limma::lmFit(design = design, block = NULL) |>
    limma::contrasts.fit(contrast) |>
    limma::eBayes(robust = TRUE)

  limma::topTable(fit, n = Inf, adjust.method = "BH") |>
    tibble::rownames_to_column("gene") |>
    tibble::as_tibble() |>
    dplyr::arrange(dplyr::desc(logFC)) |>
    dplyr::rename(avg_log2FC = logFC, p_val_adj = adj.P.Val)
}

run_deg_combined <- function(object, lookup, condition_1, condition_2) {
  pseudobulk <- Libra::to_pseudobulk(
    object,
    cell_type_col = "orig.ident",
    label_col = "diagnosis",
    replicate_col = "patient"
  )
  stopifnot(is.list(pseudobulk), "PNP" %in% names(pseudobulk))
  run_deg_limma(
    pseudobulk[["PNP"]],
    lookup,
    condition_1,
    condition_2
  )
}

collect_deg_cluster_results <- function(clusters, results) {
  stopifnot(
    length(clusters) == length(results),
    length(clusters) > 0L,
    !anyDuplicated(clusters),
    all(vapply(results, is.data.frame, logical(1)))
  )
  stats::setNames(results, clusters)
}

deg_output_paths <- function(comparison, root = "results/targets/de") {
  prefix <- file.path(root, paste0("de_", comparison))
  c(
    cluster_workbook = paste0(prefix, "_cluster.xlsx"),
    combined_workbook = paste0(prefix, "_combined.xlsx"),
    plot = paste0(prefix, ".pdf"),
    eligibility = paste0(prefix, "_cluster_eligibility.csv")
  )
}

write_deg_cluster_workbook <- function(results, path) {
  ensure_parent_dir(path)
  stopifnot(length(results) > 0L, !anyDuplicated(names(results)))
  writexl::write_xlsx(results, path)
  path
}

write_deg_combined_workbook <- function(result, path) {
  ensure_parent_dir(path)
  stopifnot(is.data.frame(result), nrow(result) > 0L)
  writexl::write_xlsx(list(combined = result), path)
  path
}

write_deg_eligibility <- function(eligibility, path) {
  ensure_parent_dir(path)
  readr::write_csv(eligibility, path)
  path
}

write_deg_count_plot <- function(results, colors, title, path) {
  ensure_parent_dir(path)
  counts <- tibble::tibble(
    cluster = names(results),
    n = vapply(
      results,
      function(result) sum(result$p_val_adj < 0.05, na.rm = TRUE),
      integer(1)
    )
  ) |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::slice_head(n = 10L)
  counts$cluster <- factor(
    counts$cluster,
    levels = rev(counts$cluster)
  )
  stopifnot(all(as.character(counts$cluster) %in% names(colors)))

  plot <- ggplot2::ggplot(
    counts,
    ggplot2::aes(x = cluster, y = n, fill = cluster)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "", y = "", title = title)
  ggplot2::ggsave(path, plot = plot, width = 3, height = 3)
  path
}
