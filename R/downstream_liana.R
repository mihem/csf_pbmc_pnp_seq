liana_result_dir <- function() {
  file.path("results", "targets", "liana")
}

liana_olink_markers <- function() {
  c(
    "CCL2", "CCL3", "CXCL8", "GZMA", "KLRD1", "SIRT2", "TNFSF14",
    "TREM2", "VEGFA", "IFNG", "CXCL12", "CXCL10", "TNFRSF9", "CCL7",
    "TNFSF12", "IL1B", "IL18"
  )
}

liana_selected_ligands <- function() {
  c("CCL7", "CXCL8", "CCL2", "CCL3", "CXCL10", "IFNG")
}

make_liana_resource <- function(markers = liana_olink_markers()) {
  resource <- liana::select_resource("Consensus")$Consensus |>
    dplyr::filter(
      .data$source_genesymbol %in% .env$markers |
        .data$target_genesymbol %in% .env$markers
    ) |>
    dplyr::arrange(
      .data$source_genesymbol,
      .data$target_genesymbol,
      .data$source,
      .data$target
    )
  stopifnot(
    nrow(resource) > 0L,
    all(c("source_genesymbol", "target_genesymbol") %in% names(resource)),
    any(resource$source_genesymbol %in% markers),
    any(resource$target_genesymbol %in% markers)
  )
  resource
}

prepare_liana_object <- function(object, seed = 123L, max_cells = 1000L) {
  stopifnot(
    inherits(object, "Seurat"),
    "cluster" %in% colnames(object[[]]),
    "RNA" %in% names(object@assays),
    !anyNA(object$cluster),
    max_cells > 0L
  )
  SeuratObject::DefaultAssay(object) <- "RNA"
  SeuratObject::Idents(object) <- object$cluster
  result <- withr::with_seed(seed, subset(object, downsample = max_cells))
  cluster_counts <- table(SeuratObject::Idents(result))
  stopifnot(
    ncol(result) > 0L,
    all(cluster_counts <= max_cells),
    setequal(names(cluster_counts), unique(as.character(object$cluster)))
  )
  result
}

run_liana_consensus <- function(object, resource, seed = 123L) {
  stopifnot(
    inherits(object, "Seurat"),
    is.data.frame(resource),
    nrow(resource) > 0L
  )
  withr::with_seed(seed, liana::liana_wrap(
    object,
    method = c("natmi", "connectome", "logfc", "sca"),
    resource = "custom",
    external_resource = resource,
    expression_pct = 0.1,
    parallel = FALSE,
    interactions = resource,
    verbose = FALSE
  ))
}

aggregate_liana_results <- function(results) {
  result <- liana::liana_aggregate(results, verbose = FALSE)
  stopifnot(
    is.data.frame(result),
    nrow(result) > 0L,
    all(c(
      "source", "target", "ligand.complex", "receptor.complex",
      "aggregate_rank"
    ) %in% names(result))
  )
  dplyr::arrange(
    result,
    .data$aggregate_rank,
    .data$source,
    .data$target,
    .data$ligand.complex,
    .data$receptor.complex
  )
}

write_liana_qs <- function(value, path) {
  ensure_parent_dir(path)
  qs::qsave(value, path, preset = "high")
  path
}

write_liana_aggregate <- function(results, path) {
  ensure_parent_dir(path)
  readr::write_csv(results, path)
  path
}

write_liana_all_dotplot <- function(results, path) {
  plot <- liana::liana_dotplot(results) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 14, face = "plain"),
      strip.text = ggplot2::element_text(size = 14)
    )
  ensure_parent_dir(path)
  ggplot2::ggsave(path, plot, width = 100, height = 10, limitsize = FALSE)
  path
}

selected_liana_results <- function(
  results,
  ligands = liana_selected_ligands(),
  target = "CD8TEM_3"
) {
  selected <- dplyr::filter(
    results,
    .data$ligand.complex %in% .env$ligands,
    .data$target == .env$target
  )
  stopifnot(nrow(selected) > 0L)
  selected
}

write_liana_selected_dotplot <- function(results, path) {
  plot <- liana::liana_dotplot(results) +
    ggplot2::coord_flip() +
    ggplot2::scale_x_discrete(limits = rev) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 10, angle = 90, hjust = 1, vjust = 0.5,
        colour = "black", face = "plain"
      ),
      axis.text.y = ggplot2::element_text(colour = "black", size = 14),
      strip.text = ggplot2::element_text(size = 14, angle = 90)
    )
  ensure_parent_dir(path)
  ggplot2::ggsave(path, plot, width = 6, height = 4)
  path
}

write_liana_selected_chord <- function(results, path) {
  selected <- dplyr::filter(
    results,
    .data$ligand.complex %in% liana_selected_ligands()
  )
  stopifnot(nrow(selected) > 0L)
  ensure_parent_dir(path)
  grDevices::pdf(path, width = 20, height = 20)
  on.exit(grDevices::dev.off(), add = TRUE)
  liana::chord_freq(selected)
  path
}

liana_interaction_keys <- function(results) {
  paste(
    results$source,
    results$target,
    results$ligand.complex,
    results$receptor.complex,
    sep = "\r"
  )
}

validate_liana_results <- function(results, legacy_path) {
  stopifnot(file.exists(legacy_path), is.list(results))
  legacy <- qs::qread(legacy_path)
  methods <- c("natmi", "connectome", "logfc", "sca")
  stopifnot(
    identical(names(results), methods),
    identical(names(legacy), methods)
  )

  dplyr::bind_rows(lapply(methods, function(method) {
    current <- results[[method]]
    baseline <- legacy[[method]]
    key_columns <- c("source", "target", "ligand.complex", "receptor.complex")
    stopifnot(
      all(key_columns %in% names(current)),
      all(key_columns %in% names(baseline)),
      nrow(current) > 0L,
      nrow(baseline) > 0L
    )
    current_keys <- unique(liana_interaction_keys(current))
    baseline_keys <- unique(liana_interaction_keys(baseline))
    overlap <- length(intersect(current_keys, baseline_keys))
    union_count <- length(union(current_keys, baseline_keys))
    tibble::tibble(
      method = method,
      current_interactions = length(current_keys),
      legacy_interactions = length(baseline_keys),
      shared_interactions = overlap,
      interaction_jaccard = overlap / union_count,
      schema_compatible = all(names(baseline) %in% names(current))
    )
  }))
}

write_liana_validation <- function(validation, path) {
  ensure_parent_dir(path)
  readr::write_csv(validation, path)
  path
}
