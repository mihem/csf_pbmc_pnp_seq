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
