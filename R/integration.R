attach_reduction_by_cell <- function(object, reduction, name, key) {
  embeddings <- Seurat::Embeddings(reduction)
  stopifnot(setequal(rownames(embeddings), colnames(object)))
  embeddings <- embeddings[colnames(object), , drop = FALSE]
  object[[name]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings,
    key = key,
    assay = "RNA"
  )
  object
}

compute_stacas_reduction <- function(object, lookup, config) {
  future::plan(
    "multicore",
    workers = as.integer(config$project$workers)
  )
  options(
    future.globals.maxSize =
      as.numeric(config$project$future_globals_gib) * 1024^3
  )
  # STACAS loads its packaged gene blocklist through data() at runtime.
  suppressPackageStartupMessages(library("STACAS", character.only = TRUE))
  object <- add_batch_metadata(object, lookup)
  object <- SeuratObject::JoinLayers(object)
  object[["RNA"]]$scale.data <- NULL
  object@reductions <- list()
  object@graphs <- list()
  object@neighbors <- list()
  split_objects <- Seurat::SplitObject(
    object,
    split.by = config$integration$split_by
  )
  rm(object)
  invisible(gc())
  stacas <- STACAS::Run.STACAS(
    split_objects,
    cell.labels = config$integration$labels,
    k.weight = as.integer(config$integration$k_weight),
    seed = as.integer(config$integration$seed)
  )
  stacas[["pca"]]
}

assemble_stacas_integration <- function(object, lookup, reduction, config) {
  object <- add_batch_metadata(object, lookup)
  object <- attach_reduction_by_cell(
    object,
    reduction,
    "stacas.ss.all",
    "stacasssall_"
  )
  Seurat::RunUMAP(
    object,
    reduction = "stacas.ss.all",
    reduction.name = "umap.stacas.ss.all",
    dims = seq_len(as.integer(config$integration$dimensions)),
    seed.use = as.integer(config$integration$umap_seed)
  )
}

write_azimuth_level2_map <- function(object, path, seed = 123L) {
  stopifnot(
    "azimuth_pbmcref2" %in% colnames(object[[]]),
    "umap.stacas.ss.all" %in% names(object@reductions)
  )
  colors <- withr::with_seed(
    seed,
    unname(Polychrome::createPalette(100, pals::cols25()))
  )
  plot <- Seurat::DimPlot(
    object,
    reduction = "umap.stacas.ss.all",
    group.by = "azimuth_pbmcref2",
    raster = FALSE,
    pt.size = 0.1,
    alpha = 0.1,
    cols = colors,
    label = TRUE
  ) + scMisc::theme_rect()
  ensure_parent_dir(path)
  ggplot2::ggsave(path, plot, width = 20, height = 8, dpi = 300)
  path
}

write_integrated_tissue_umap <- function(object, path, seed = 123L) {
  stopifnot(
    "tissue" %in% colnames(object[[]]),
    "umap.stacas.ss.all" %in% names(object@reductions)
  )
  colors <- withr::with_seed(
    seed,
    unname(Polychrome::createPalette(100, pals::cols25()))
  )
  plot <- Seurat::DimPlot(
    object,
    reduction = "umap.stacas.ss.all",
    pt.size = 0.5,
    raster = FALSE,
    alpha = 0.1,
    group.by = "tissue",
    cols = colors
  ) +
    scMisc::theme_rect() +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2")
  ensure_parent_dir(path)
  ggplot2::ggsave(path, plot = plot, width = 20, height = 8, dpi = 300)
  path
}

stabilize_annotation_input <- function(object, path, validation_result) {
  stopifnot(isTRUE(validation_result$passed[[1]]))
  baseline <- qs::qread(path, nthreads = 6)
  stopifnot(identical(colnames(object), colnames(baseline)))
  object[["stacas.ss.all"]] <- baseline[["stacas.ss.all"]]
  object[["umap.stacas.ss.all"]] <- baseline[["umap.stacas.ss.all"]]
  object$group <- baseline$group
  object$age <- baseline$age
  object
}
