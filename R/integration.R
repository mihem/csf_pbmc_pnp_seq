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

apply_annotation_checkpoint <- function(object, path) {
  checkpoint <- qs::qread(path, nthreads = 6)
  stopifnot(
    identical(names(checkpoint), c("stacas.ss.all", "umap.stacas.ss.all")),
    inherits(checkpoint[["stacas.ss.all"]], "DimReduc"),
    inherits(checkpoint[["umap.stacas.ss.all"]], "DimReduc"),
    identical(
      rownames(Seurat::Embeddings(checkpoint[["stacas.ss.all"]])),
      colnames(object)
    ),
    identical(
      rownames(Seurat::Embeddings(checkpoint[["umap.stacas.ss.all"]])),
      colnames(object)
    )
  )
  object[["stacas.ss.all"]] <- checkpoint[["stacas.ss.all"]]
  object[["umap.stacas.ss.all"]] <- checkpoint[["umap.stacas.ss.all"]]
  object
}

prepare_annotation_input <- function(object, lookup, checkpoint_path) {
  object |>
    add_batch_metadata(lookup) |>
    apply_annotation_checkpoint(checkpoint_path)
}
