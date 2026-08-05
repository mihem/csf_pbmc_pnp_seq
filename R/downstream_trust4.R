sural_trust4_result_dir <- function() {
  file.path("results", "targets", "trust4")
}

discover_sural_trust4_inputs <- function(raw_dir, sample_map) {
  files <- list.files(
    raw_dir,
    pattern = "_TRUST4_output_barcode_report\\.tsv$",
    full.names = TRUE
  )
  files <- sort(files)
  library_id <- sub(
    "_TRUST4_output_barcode_report\\.tsv$", "", basename(files)
  )
  sample_map <- unlist(sample_map, use.names = TRUE)

  stopifnot(
    dir.exists(raw_dir),
    length(files) > 0L,
    !anyDuplicated(library_id),
    !is.null(names(sample_map)),
    setequal(library_id, names(sample_map))
  )

  tibble::tibble(
    library_id = library_id,
    sample = unname(sample_map[library_id]),
    file = files
  )
}

sural_trust4_input_files <- function(manifest) {
  stopifnot(
    identical(colnames(manifest), c("library_id", "sample", "file")),
    all(file.exists(manifest$file))
  )
  unname(manifest$file)
}

read_sural_metadata_umap <- function(path, reduction) {
  object <- qs::qread(path)
  stopifnot(
    inherits(object, "Seurat"),
    reduction %in% names(object@reductions)
  )

  metadata <- object@meta.data |>
    tibble::rownames_to_column("cell_id")
  if (!"sample" %in% names(metadata)) {
    metadata$sample <- sub("_.*$", "", metadata$cell_id)
  }
  umap <- as.data.frame(object@reductions[[reduction]]@cell.embeddings) |>
    tibble::rownames_to_column("cell_id")
  stopifnot(
    !anyDuplicated(metadata$cell_id),
    !anyDuplicated(umap$cell_id),
    setequal(metadata$cell_id, umap$cell_id),
    ncol(umap) == 3L
  )
  names(umap)[2:3] <- c("UMAP_1", "UMAP_2")

  dplyr::left_join(metadata, umap, by = "cell_id")
}

classify_trust4_receptor <- function(cell_type) {
  dplyr::case_when(
    cell_type == "B" ~ "BCR",
    grepl("T$", cell_type) ~ "TCR",
    TRUE ~ "Other"
  )
}

map_sural_trust4_barcodes <- function(
  manifest, input_files, cells, cluster_column
) {
  stopifnot(
    length(input_files) == nrow(manifest),
    identical(unname(manifest$file), input_files),
    all(c("cell_id", "sample", cluster_column) %in% names(cells))
  )

  reports <- purrr::map2_dfr(
    input_files,
    seq_along(input_files),
    function(path, index) {
      report <- readr::read_tsv(
        path,
        show_col_types = FALSE,
        name_repair = "minimal"
      )
      names(report)[1L] <- "barcode"
      report |>
        dplyr::mutate(
          library_id = manifest$library_id[[index]],
          sample = manifest$sample[[index]],
          .before = 1L
        )
    }
  ) |>
    dplyr::rename(trust4_cell_type = "cell_type") |>
    dplyr::mutate(
      barcode = sub("-[0-9]+$", "", .data$barcode),
      receptor = classify_trust4_receptor(.data$trust4_cell_type)
    )

  cell_index <- cells |>
    dplyr::transmute(
      cell_id,
      sample = as.character(.data$sample),
      barcode = sub(
        "-[0-9]+$", "", sub("^[^_]+_", "", .data$cell_id)
      ),
      seurat_cluster = as.character(.data[[cluster_column]]),
      UMAP_1,
      UMAP_2
    )

  stopifnot(
    !anyDuplicated(reports[c("library_id", "barcode")]),
    !anyDuplicated(cell_index[c("sample", "barcode")])
  )

  mapped <- reports |>
    dplyr::left_join(
      cell_index,
      by = c("sample", "barcode"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(
      match_status = ifelse(is.na(.data$cell_id), "unmatched", "matched"),
      .after = "barcode"
    )

  matched_cells <- mapped |>
    dplyr::filter(.data$match_status == "matched") |>
    dplyr::distinct(
      .data$cell_id, .data$receptor, .data$seurat_cluster,
      .data$UMAP_1, .data$UMAP_2
    )
  summary <- matched_cells |>
    dplyr::count(
      receptor, seurat_cluster, name = "cell_count", .drop = FALSE
    ) |>
    dplyr::arrange(.data$receptor, dplyr::desc(.data$cell_count))
  match_summary <- mapped |>
    dplyr::count(library_id, sample, receptor, match_status) |>
    tidyr::pivot_wider(
      names_from = match_status,
      values_from = n,
      values_fill = 0
    ) |>
    dplyr::arrange(.data$library_id, .data$receptor)

  list(
    records = mapped,
    matched_cells = matched_cells,
    cluster_summary = summary,
    match_summary = match_summary
  )
}

write_sural_trust4_table <- function(
  mapped, file_name = "trust4_cells_with_seurat_clusters.xlsx"
) {
  path <- file.path(sural_trust4_result_dir(), file_name)
  ensure_parent_dir(path)
  writexl::write_xlsx(
    list(
      matched_cells = dplyr::filter(
        mapped$records, .data$match_status == "matched"
      ),
      all_trust4_records = mapped$records,
      cluster_summary = mapped$cluster_summary,
      match_summary = mapped$match_summary
    ),
    path
  )
  path
}

write_sural_trust4_umap <- function(
  cells, mapped, file_name = "trust4_cells_umap.png",
  title = "Cells with TRUST4 receptor calls"
) {
  path <- file.path(sural_trust4_result_dir(), file_name)
  ensure_parent_dir(path)
  highlighted <- mapped$matched_cells |>
    dplyr::mutate(
      receptor = factor(.data$receptor, levels = c("TCR", "BCR", "Other"))
    )
  colors <- c(TCR = "#D1495B", BCR = "#00798C", Other = "#F2A541")
  plot <- ggplot2::ggplot(
    cells,
    ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2)
  ) +
    ggplot2::geom_point(color = "grey75", size = 0.08, alpha = 0.2) +
    ggplot2::geom_point(
      data = highlighted,
      ggplot2::aes(color = .data$receptor),
      size = 1.1,
      alpha = 0.9
    ) +
    ggplot2::scale_color_manual(values = colors, drop = TRUE) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      x = "UMAP 1",
      y = "UMAP 2",
      color = NULL
    ) +
    ggplot2::theme_classic()
  ggplot2::ggsave(path, plot, width = 7, height = 6, dpi = 300)
  path
}

write_sural_trust4_cluster_summary <- function(
  mapped, file_name = "trust4_receptor_cells_by_cluster.pdf",
  title = "TRUST4 receptor-positive cells by Seurat cluster"
) {
  path <- file.path(sural_trust4_result_dir(), file_name)
  ensure_parent_dir(path)
  data <- mapped$cluster_summary |>
    dplyr::filter(.data$receptor %in% c("TCR", "BCR")) |>
    dplyr::mutate(
      seurat_cluster = stats::reorder(
        .data$seurat_cluster, .data$cell_count, FUN = sum
      )
    )
  colors <- c(TCR = "#D1495B", BCR = "#00798C")
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$seurat_cluster,
      y = .data$cell_count,
      fill = .data$receptor
    )
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$cell_count),
      vjust = -0.3,
      size = 3
    ) +
    ggplot2::facet_wrap(~receptor, ncol = 1L, scales = "free") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_pretty(),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "Number of cells"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1, vjust = 1
      )
    )
  ggplot2::ggsave(path, plot, width = 9, height = 6)
  path
}
