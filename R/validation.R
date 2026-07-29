read_legacy_object <- function(path) {
  qs::qread(path, nthreads = 6)
}

assert_same_cells_features <- function(object, baseline, stage) {
  stopifnot(
    identical(dim(object), dim(baseline)),
    identical(rownames(object), rownames(baseline)),
    identical(colnames(object), colnames(baseline))
  )
  tibble::tibble(stage, check = "cells_and_features", passed = TRUE)
}

validate_sc_list <- function(object, path) {
  baseline <- read_legacy_object(path)
  stopifnot(
    identical(names(object), names(baseline)),
    identical(vapply(object, ncol, numeric(1)), vapply(baseline, ncol, numeric(1))),
    identical(lapply(object, colnames), lapply(baseline, colnames))
  )
  tibble::tibble(stage = "sc_list", check = "sample_cells", passed = TRUE)
}

validate_preprocessed_merge <- function(object, path) {
  baseline <- read_legacy_object(path)
  result <- assert_same_cells_features(object, baseline, "sc_merge_pre")
  columns <- c(
    "orig.ident", "nCount_RNA", "nFeature_RNA", "pseudonym", "patient",
    "percent_mt", "tissue", "sample", "sex", "diagnosis",
    "incat_at_lumbar_puncture", "incat_follow_up",
    "onls_at_lumbar_puncture", "onls_follow_up", "icu"
  )
  stopifnot(
    all(columns %in% colnames(object[[]])),
    all(vapply(columns, function(column) {
      isTRUE(all.equal(
        object[[column, drop = TRUE]],
        baseline[[column, drop = TRUE]],
        check.attributes = FALSE,
        tolerance = 1e-12
      ))
    }, logical(1))),
    identical(SeuratObject::Layers(object), SeuratObject::Layers(baseline))
  )
  for (layer in SeuratObject::Layers(baseline)) {
    stopifnot(identical(
      SeuratObject::LayerData(object, assay = "RNA", layer = layer),
      SeuratObject::LayerData(baseline, assay = "RNA", layer = layer)
    ))
  }
  result
}

validate_preprocessed_reductions <- function(object, path) {
  baseline <- read_legacy_object(path)
  result <- assert_same_cells_features(object, baseline, "sc_preprocessed")
  pca <- Seurat::Embeddings(object, "pca")
  pca_baseline <- Seurat::Embeddings(baseline, "pca")
  component_correlations <- diag(stats::cor(pca, pca_baseline))
  stopifnot(
    min(abs(component_correlations)) > 0.99,
    length(intersect(
      Seurat::VariableFeatures(object),
      Seurat::VariableFeatures(baseline)
    )) >= 1998L,
    isTRUE(all.equal(
      Seurat::Embeddings(object, "nmf"),
      Seurat::Embeddings(baseline, "nmf"),
      tolerance = 1e-12,
      check.attributes = FALSE
    ))
  )
  result
}

validate_azimuth <- function(object, path) {
  baseline <- read_legacy_object(path)
  result <- assert_same_cells_features(object, baseline, "sc_azimuth")
  for (level in 1:3) {
    source_label <- paste0("predicted.celltype.l", level)
    source_score <- paste0(source_label, ".score")
    target_label <- paste0("azimuth_pbmcref", level)
    target_score <- paste0(target_label, "_score")
    expected <- as.character(baseline[[source_label, drop = TRUE]])
    score <- baseline[[source_score, drop = TRUE]]
    expected[score < 0.4] <- "unknown"
    stopifnot(
      identical(as.character(object[[target_label, drop = TRUE]]), expected),
      isTRUE(all.equal(
        object[[target_score, drop = TRUE]], score,
        check.attributes = FALSE, tolerance = 1e-12
      ))
    )
  }
  result
}

validate_integration <- function(object, path) {
  baseline <- read_legacy_object(path)
  result <- assert_same_cells_features(object, baseline, "sc_integrated")
  for (reduction in c("stacas.ss.all", "umap.stacas.ss.all")) {
    stopifnot(
      reduction %in% Seurat::Reductions(object),
      identical(
        dim(Seurat::Embeddings(object, reduction)),
        dim(Seurat::Embeddings(baseline, reduction))
      ),
      all(is.finite(Seurat::Embeddings(object, reduction)))
    )
  }
  set.seed(42L)
  cells <- sample(colnames(object), 5000L)
  compare_reduction <- function(reduction) {
    x <- scale(Seurat::Embeddings(object, reduction)[cells, , drop = FALSE])
    y <- scale(Seurat::Embeddings(baseline, reduction)[cells, , drop = FALSE])
    rotation <- svd(crossprod(x, y))
    aligned <- x %*% rotation$u %*% t(rotation$v)
    neighbors_x <- RANN::nn2(x, k = 31L)$nn.idx[, -1, drop = FALSE]
    neighbors_y <- RANN::nn2(y, k = 31L)$nn.idx[, -1, drop = FALSE]
    jaccard <- vapply(seq_len(nrow(neighbors_x)), function(index) {
      length(intersect(neighbors_x[index, ], neighbors_y[index, ])) /
        length(union(neighbors_x[index, ], neighbors_y[index, ]))
    }, numeric(1))
    c(correlation = stats::cor(c(aligned), c(y)), jaccard = stats::median(jaccard))
  }
  stacas <- compare_reduction("stacas.ss.all")
  umap <- compare_reduction("umap.stacas.ss.all")
  stopifnot(
    stacas[["correlation"]] > 0.99,
    stacas[["jaccard"]] > 0.9,
    umap[["correlation"]] > 0.85,
    umap[["jaccard"]] > 0.5
  )
  result
}

validate_annotation <- function(object, path) {
  baseline <- read_legacy_object(path)
  result <- assert_same_cells_features(object, baseline, "sc_annotated")
  stopifnot(
    identical(levels(object$cluster), levels(baseline$cluster)),
    identical(table(object$cluster), table(baseline$cluster)),
    identical(as.character(object$cluster), as.character(baseline$cluster)),
    identical(as.character(Seurat::Idents(object)), as.character(object$cluster)),
    identical(object@misc, baseline@misc)
  )
  result
}
