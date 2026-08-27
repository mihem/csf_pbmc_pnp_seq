tcr_cluster_names <- function() {
  c(
    "CD4naive_1", "CD4TCM_1", "CD4TCM_2", "CD4TEM", "Treg", "MAIT",
    "gdT", "CD4CTL", "CD8naive", "CD8TCM", "CD8TEM_1", "CD8TEM_2",
    "CD8TEM_3", "NKCD56bright_1", "NKCD56bright_2", "NKCD56dim"
  )
}

discover_tcr_inputs <- function(raw_tcr_dir = file.path("raw", "tcr")) {
  files <- list.files(
    raw_tcr_dir,
    pattern = "filtered_contig_annotations\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  files <- sort(files)
  library_ids <- basename(dirname(dirname(files)))
  expected_ids <- c(
    paste0("CSF_pool_", seq_len(9L)),
    paste0("PBMC_pool_", seq_len(9L)),
    "PNP38"
  )

  stopifnot(
    dir.exists(raw_tcr_dir),
    length(files) == 19L,
    !anyDuplicated(library_ids),
    setequal(library_ids, expected_ids),
    all(basename(dirname(files)) == "outs")
  )

  tibble::tibble(library_id = library_ids, file = files)
}

tcr_input_files <- function(manifest) {
  stopifnot(
    identical(colnames(manifest), c("library_id", "file")),
    nrow(manifest) == 19L,
    all(file.exists(manifest$file))
  )
  unname(manifest$file)
}

import_tcr_contigs <- function(manifest, input_files, donor_assignments, lookup) {
  stopifnot(
    length(input_files) == nrow(manifest),
    identical(unname(manifest$file), input_files),
    all(setdiff(manifest$library_id, "PNP38") %in% names(donor_assignments)),
    !anyDuplicated(lookup$pseudonym)
  )

  contigs <- lapply(input_files, utils::read.csv)
  names(contigs) <- manifest$library_id
  patient_lookup <- lookup[, c("pseudonym", "patient")] |>
    dplyr::add_row(pseudonym = "unassigned", patient = "unassigned") |>
    dplyr::add_row(pseudonym = "doublet", patient = "doublet")

  multiplexed_ids <- setdiff(manifest$library_id, "PNP38")
  for (library_id in multiplexed_ids) {
    assignments <- donor_assignments[[library_id]]
    stopifnot(!anyDuplicated(assignments$cell))
    assignment_index <- match(contigs[[library_id]]$barcode, assignments$cell)
    pseudonym <- assignments$donor_id[assignment_index]
    patient_index <- match(pseudonym, patient_lookup$pseudonym)
    contigs[[library_id]]$patient <- patient_lookup$patient[patient_index]
  }

  split_contigs <- lapply(multiplexed_ids, function(library_id) {
    pieces <- split(contigs[[library_id]], contigs[[library_id]]$patient)
    pieces[c("doublet", "unassigned")] <- NULL
    pieces[!is.na(names(pieces))]
  })
  names(split_contigs) <- multiplexed_ids
  split_contigs <- purrr::list_flatten(split_contigs)
  split_contigs[["PNP38"]] <- contigs[["PNP38"]]

  names(split_contigs) <- sub(
    pattern = "([^_]+).*_([^_]+)$",
    replacement = "\\1_\\2",
    names(split_contigs)
  )
  names(split_contigs)[names(split_contigs) == "PNP38"] <- "CSF_P14"
  split_contigs <- split_contigs[order(names(split_contigs))]

  stopifnot(
    length(split_contigs) == 61L,
    !anyDuplicated(names(split_contigs)),
    identical(names(split_contigs), sort(names(split_contigs)))
  )
  split_contigs
}

combine_tcr_contigs <- function(contigs) {
  sample_ids <- sub(".*_([^_]+)$", "\\1", names(contigs))
  tissues <- sub("([^_]+)_.*", "\\1", names(contigs))
  scRepertoire::combineTCR(
    contigs,
    samples = tissues,
    ID = sample_ids,
    removeNA = TRUE,
    removeMulti = TRUE
  )
}

annotate_tcr_cells <- function(sc_annotated, combined_tcr) {
  cells <- colnames(sc_annotated)[
    as.character(sc_annotated$cluster) %in% tcr_cluster_names()
  ]
  object <- subset(sc_annotated, cells = cells)
  object <- scRepertoire::combineExpression(
    combined_tcr,
    object,
    cloneCall = "aa",
    group.by = "sample",
    proportion = FALSE,
    cloneSize = c(
      Single = 1,
      Small = 5,
      Medium = 20,
      Large = 100,
      Hyperexpanded = 500
    )
  )
  object
}

tcr_result_dir <- function() {
  file.path("results", "targets", "tcr")
}

make_tcr_report_tables <- function(sc_tcr) {
  metadata <- sc_tcr[[]]
  ctaa_sample <- dplyr::count(metadata, .data$CTaa, .data$sample) |>
    tidyr::drop_na("CTaa") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "n", values_fill = 0) |>
    dplyr::mutate(total = rowSums(dplyr::pick(where(is.numeric)))) |>
    dplyr::arrange(dplyr::desc(.data$total), .data$CTaa) |>
    dplyr::select(-"total")
  ctaa_sample <- ctaa_sample[, c("CTaa", sort(setdiff(names(ctaa_sample), "CTaa")))]

  clone_size_count <- metadata |>
    dplyr::filter(!is.na(.data$CTaa)) |>
    dplyr::group_by(.data$cloneSize) |>
    dplyr::summarise(
      n_clones = dplyr::n_distinct(.data$CTaa), n_cells = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$cloneSize)

  shared <- metadata |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::select(
      "sample", "patient", "CTaa", "tissue_diagnosis", "cloneSize",
      "clonalFrequency"
    ) |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::group_by(.data$CTaa) |>
    dplyr::filter(dplyr::n() > 1L) |>
    dplyr::arrange(.data$CTaa, .data$sample) |>
    dplyr::ungroup()
  shared_summary <- shared |>
    dplyr::group_by(.data$CTaa) |>
    dplyr::mutate(
      sample_count = dplyr::n(),
      patients_count = dplyr::n_distinct(.data$patient),
      sample_num = paste0("sample_", dplyr::row_number())
    ) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(
      id_cols = c("CTaa", "sample_count", "patients_count"),
      names_from = "sample_num",
      values_from = c(
        "sample", "patient", "tissue_diagnosis", "cloneSize",
        "clonalFrequency"
      ),
      names_sep = "_"
    ) |>
    dplyr::arrange(
      dplyr::desc(.data$patients_count), dplyr::desc(.data$sample_count),
      .data$CTaa
    )

  shared_csf_pbmc <- metadata |>
    dplyr::select("CTaa", "tissue", "diagnosis", "clonalFrequency", "patient") |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::group_by(.data$CTaa) |>
    dplyr::filter(dplyr::n() > 1L) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(
      names_from = "tissue", values_from = "clonalFrequency",
      names_prefix = "freq_"
    ) |>
    tidyr::drop_na("freq_CSF", "freq_PBMC") |>
    dplyr::filter(.data$freq_CSF > 1 | .data$freq_PBMC > 1) |>
    dplyr::mutate(
      log2_ratio = log2(.data$freq_CSF / .data$freq_PBMC),
      enriched_in = dplyr::case_when(
        .data$log2_ratio > 1 ~ "CSF",
        .data$log2_ratio < -1 ~ "PBMC",
        TRUE ~ "Similar"
      )
    )

  shared_between_patients <- metadata |>
    dplyr::select(
      "patient", "CTaa", "tissue", "diagnosis", "cloneSize",
      "clonalFrequency"
    ) |>
    tidyr::drop_na() |>
    dplyr::group_by(.data$CTaa) |>
    dplyr::reframe(
      patients = paste(sort(unique(.data$patient)), collapse = ","),
      n_patients = dplyr::n_distinct(.data$patient),
      diagnosis = paste(sort(unique(.data$diagnosis)), collapse = ","),
      tissue = paste(sort(unique(.data$tissue)), collapse = ","),
      clone_size = unique(.data$cloneSize),
      clonal_frequency = unique(.data$clonalFrequency)
    ) |>
    dplyr::mutate(shared = ifelse(.data$n_patients > 1L, "shared", "unique")) |>
    dplyr::filter(.data$n_patients > 1L) |>
    dplyr::arrange(dplyr::desc(.data$clonal_frequency), .data$CTaa)

  list(
    CTaa_sample = ctaa_sample,
    cloneSize_count = clone_size_count,
    shared_clones_summary = shared_summary,
    shared_clones_csf_pbmc_abundance = shared_csf_pbmc,
    shared_clones_between_patients = shared_between_patients
  )
}

write_tcr_report_tables <- function(tables) {
  root <- tcr_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(root, paste0(names(tables), ".xlsx"))
  purrr::walk2(tables, paths, writexl::write_xlsx)
  unname(paths)
}

tcr_clone_size_colors <- function(object) {
  levels <- levels(object$cloneSize)
  stats::setNames(rev(viridis::turbo(length(levels))), levels)
}

make_tcr_abundance_table <- function(metadata, x_axis) {
  absolute <- table(metadata$cloneSize, metadata[[x_axis]]) |>
    as.data.frame.matrix() |>
    tibble::rownames_to_column("cell")
  percentage <- absolute |>
    dplyr::mutate(dplyr::across(
      tidyselect::where(is.numeric),
      ~ round(.x / sum(.x) * 100, 2)
    ))
  list(absolute = absolute, percentage = percentage)
}

write_tcr_abundance_tables <- function(sc_tcr) {
  root <- abundance_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  metadata <- sc_tcr[[]]
  tables <- list(
    cluster = make_tcr_abundance_table(metadata, "cluster"),
    sample = make_tcr_abundance_table(metadata, "sample"),
    tissue_group = make_tcr_abundance_table(metadata, "tissue_group")
  )
  paths <- file.path(
    root,
    paste0("abundance_tbl_sc_tcr_", names(tables), ".xlsx")
  )
  for (index in seq_along(tables)) {
    writexl::write_xlsx(tables[[index]], paths[[index]])
  }
  paths
}

make_tcr_stacked_plot <- function(metadata, x_axis, x_order, colors) {
  data <- table(metadata$cloneSize, metadata[[x_axis]]) |>
    as.data.frame.matrix() |>
    tibble::rownames_to_column("cloneSize") |>
    tidyr::pivot_longer(-cloneSize, names_to = "group", values_to = "count")
  data$cloneSize <- factor(data$cloneSize, levels = names(colors))
  data$group <- factor(data$group, levels = x_order)
  data <- dplyr::filter(data, .data$count != 0)
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$group, y = .data$count, fill = .data$cloneSize)
  ) +
    ggplot2::geom_col(color = "black", linewidth = 0.1, position = "fill") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL)) +
    ggplot2::labs(x = NULL, y = "Proportion of cells") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.3)
    )
}

write_tcr_abundance_plots <- function(sc_tcr) {
  root <- abundance_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  metadata <- sc_tcr[[]]
  diagnosis_order <- sc_tcr@misc$diagnosis_order
  main <- metadata$diagnosis %in% c("CTRL", "CIAP", "CIDP", "GBS")
  configs <- list(
    list("sc_tcr_cluster", rep(TRUE, nrow(metadata)), "cluster",
         sc_tcr@misc$cluster_order, 5, 3),
    list("sc_tcr_sample", rep(TRUE, nrow(metadata)), "sample",
         unique(metadata$sample), 10, 3),
    list("sc_tcr_tissue_group", rep(TRUE, nrow(metadata)), "tissue_group",
         unique(metadata$tissue_group), 5, 5),
    list("sc_tcr_csf_diagnosis", metadata$tissue == "CSF", "diagnosis",
         diagnosis_order, 5, 3),
    list("sc_tcr_pbmc_diagnosis", metadata$tissue == "PBMC", "diagnosis",
         diagnosis_order, 5, 3),
    list(
      "sc_tcr_csf_cd8tem_3_groups_diagnosis",
      metadata$tissue == "CSF" & metadata$cluster == "CD8TEM_3" &
        metadata$diagnosis %in% c("CIAP", "CIDP", "GBS"),
      "diagnosis", diagnosis_order, 3.2, 3
    ),
    list(
      "sc_tcr_pbmc_cd8tem_3_groups_diagnosis",
      metadata$tissue == "PBMC" & metadata$cluster == "CD8TEM_3" &
        metadata$diagnosis %in% c("CIAP", "CIDP", "GBS"),
      "diagnosis", diagnosis_order, 3.2, 3
    ),
    list("sc_tcr_main_groups_csf_diagnosis", main & metadata$tissue == "CSF",
         "diagnosis", diagnosis_order, 4, 4),
    list("sc_tcr_main_groups_pbmc_diagnosis", main & metadata$tissue == "PBMC",
         "diagnosis", diagnosis_order, 4, 4),
    list("sc_tcr_main_groups_tissue_diagnosis", main, "tissue_diagnosis",
         sc_tcr@misc$tissue_diagnosis_order, 5, 5)
  )
  colors <- tcr_clone_size_colors(sc_tcr)
  paths <- vapply(configs, function(config) {
    path <- file.path(root, paste0("stacked_barplot_", config[[1L]], ".pdf"))
    plot <- make_tcr_stacked_plot(
      metadata[config[[2L]], , drop = FALSE],
      config[[3L]],
      config[[4L]],
      colors
    )
    ggplot2::ggsave(path, plot, width = config[[5L]], height = config[[6L]])
    path
  }, character(1))
  unname(paths)
}

write_tcr_basic_plots <- function(combined_tcr, sc_tcr, seed = 42L) {
  root <- tcr_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  save_plot <- function(name, plot, width, height) {
    path <- file.path(root, name)
    ggplot2::ggsave(path, plot = plot, width = width, height = height)
    path
  }
  sample_colors <- Polychrome::createPalette(
    length(combined_tcr), Polychrome::palette36.colors(36L)
  )
  names(sample_colors) <- names(combined_tcr)
  quant <- scRepertoire::clonalQuant(combined_tcr, cloneCall = "aa", scale = FALSE) +
    ggplot2::scale_fill_manual(values = sample_colors) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  length_plot <- scRepertoire::clonalLength(combined_tcr, cloneCall = "aa")

  clone_labels <- levels(sc_tcr$cloneSize)
  clone_cols <- stats::setNames(rev(viridis::turbo(length(clone_labels))), clone_labels)
  umap <- Seurat::DimPlot(
    sc_tcr, group.by = "cloneSize", pt.size = 0.1,
    reduction = "umap.stacas.ss.all", raster = FALSE, cols = clone_cols
  ) + scMisc::theme_rect() + ggplot2::labs(x = "UMAP1", y = "UMAP2", title = NULL)
  umap[[1]]$layers[[1]]$aes_params$alpha <- dplyr::case_when(
    is.na(sc_tcr$cloneSize) ~ 0.01,
    sc_tcr$cloneSize == "Single (0 < X <= 1)" ~ 0.02,
    sc_tcr$cloneSize == "Small (1 < X <= 5)" ~ 0.2,
    sc_tcr$cloneSize == "Medium (5 < X <= 20)" ~ 0.3,
    sc_tcr$cloneSize == "Large (20 < X <= 100)" ~ 0.5,
    sc_tcr$cloneSize == "Hyperexpanded (X > 100)" ~ 1,
    TRUE ~ 0.01
  )

  c(
    save_plot("tcr_quant_abs.pdf", quant, 7, 7),
    save_plot("tcr_length_aa.pdf", length_plot, 7, 7),
    save_plot("umap_cloneSize.pdf", umap, 10, 7)
  )
}

make_tcr_shared_clone_correlation_plot <- function(
  shared_clones,
  lookup,
  variable,
  x_label,
  colors,
  jitter_width,
  log_x = FALSE,
  seed = 42L
) {
  stopifnot(
    variable %in% names(lookup),
    all(c("patient", "diagnosis", "log2_ratio") %in% names(shared_clones))
  )
  clone_data <- shared_clones |>
    dplyr::left_join(
      dplyr::select(lookup, "patient", value = tidyselect::all_of(variable)),
      by = "patient"
    ) |>
    dplyr::filter(!is.na(.data$value))
  if (log_x) clone_data <- dplyr::filter(clone_data, .data$value > 0)

  patient_data <- clone_data |>
    dplyr::group_by(.data$patient, .data$diagnosis, .data$value) |>
    dplyr::summarise(log2_ratio = mean(.data$log2_ratio), .groups = "drop")
  test <- stats::cor.test(
    patient_data$value,
    patient_data$log2_ratio,
    method = "spearman",
    exact = FALSE
  )
  subtitle <- sprintf(
    "Patient-level Spearman rho = %.2f, p = %.2g, n = %d",
    unname(test$estimate), test$p.value, nrow(patient_data)
  )

  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = clone_data,
      ggplot2::aes(
        x = .data$value,
        y = .data$log2_ratio,
        color = .data$diagnosis
      ),
      size = 1,
      alpha = 0.5,
      position = ggplot2::position_jitter(
        width = jitter_width, height = 0, seed = seed
      )
    ) +
    ggplot2::geom_point(
      data = patient_data,
      ggplot2::aes(
        x = .data$value,
        y = .data$log2_ratio,
        color = .data$diagnosis
      ),
      size = 2.5
    ) +
    ggplot2::geom_smooth(
      data = patient_data,
      ggplot2::aes(x = .data$value, y = .data$log2_ratio, group = 1),
      method = "lm", se = TRUE, linewidth = 0.5, color = "black"
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = x_label,
      y = "log2(CSF / blood frequency)",
      subtitle = subtitle,
      caption = "Small points: clonotypes; large points and fit: patient means"
    ) +
    ggplot2::theme_classic()
  if (log_x) plot <- plot + ggplot2::scale_x_log10()
  plot
}

write_tcr_shared_clone_plots <- function(sc_tcr, lookup, tables) {
  root <- tcr_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  summary <- tables$shared_clones_summary
  total_unique <- dplyr::n_distinct(stats::na.omit(sc_tcr$CTaa))
  shared_plot <- tibble::tibble(
    category = c("unique", "shared_between_patients", "shared_within_patient"),
    count = c(
      total_unique - nrow(summary),
      sum(summary$patients_count > 1),
      sum(summary$patients_count == 1)
    )
  ) |>
    ggplot2::ggplot(ggplot2::aes(x = "", y = .data$count, fill = .data$category)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$count),
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::coord_polar("y") +
    ggplot2::theme_void()
  enrichment_plot <- tables$shared_clones_csf_pbmc_abundance |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$diagnosis,
        y = .data$log2_ratio,
        fill = .data$diagnosis
      )
    ) +
    ggplot2::geom_violin(alpha = 0.5, width = 1.5) +
    ggplot2::geom_boxplot(alpha = 0.5, width = 0.2) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey30"
    ) +
    ggplot2::scale_fill_manual(values = sc_tcr@misc$diagnosis_col) +
    ggplot2::scale_y_continuous(breaks = seq(-10, 10)) +
    ggplot2::labs(
      x = "",
      y = "log2(CSF / blood frequency)",
      title = "Shared Clone Enrichment",
      subtitle = "Positive = enriched in CSF, Negative = enriched in blood"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")

  duration_plot <- make_tcr_shared_clone_correlation_plot(
    tables$shared_clones_csf_pbmc_abundance,
    lookup,
    "disease_duration_in_months",
    "Disease duration (months, log scale)",
    sc_tcr@misc$diagnosis_col,
    jitter_width = 0.06,
    log_x = TRUE
  )
  protein_plot <- make_tcr_shared_clone_correlation_plot(
    tables$shared_clones_csf_pbmc_abundance,
    lookup,
    "csf_protein",
    "CSF protein (mg/L)",
    sc_tcr@misc$diagnosis_col,
    jitter_width = 10
  )
  paths <- file.path(
    root,
    c(
      "tcr_shared_clones_summary.pdf",
      "shared_clones_enrichment_ratio_expanded.pdf",
      "shared_clones_disease_duration_correlation.pdf",
      "shared_clones_csf_protein_correlation.pdf"
    )
  )
  ggplot2::ggsave(paths[[1L]], shared_plot, width = 5, height = 5)
  ggplot2::ggsave(paths[[2L]], enrichment_plot, width = 5, height = 3)
  ggplot2::ggsave(paths[[3L]], duration_plot, width = 5, height = 3)
  ggplot2::ggsave(paths[[4L]], protein_plot, width = 5, height = 3)
  paths
}

write_tcr_albumin_quotient_plot <- function(sc_tcr, lookup, tables) {
  plot <- make_tcr_shared_clone_correlation_plot(
    tables$shared_clones_csf_pbmc_abundance,
    lookup,
    "albumin_quotient",
    "Albumin quotient (QAlb)",
    sc_tcr@misc$diagnosis_col,
    jitter_width = 0.2
  )
  path <- file.path(
    tcr_result_dir(), "shared_clones_albumin_quotient_correlation.pdf"
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = 5, height = 3)
  path
}

write_tcr_csf_cell_count_plot <- function(sc_tcr, lookup, tables) {
  plot <- make_tcr_shared_clone_correlation_plot(
    tables$shared_clones_csf_pbmc_abundance,
    lookup,
    "cell_count",
    "CSF cell count (cells/µL)",
    sc_tcr@misc$diagnosis_col,
    jitter_width = 0.1
  )
  path <- file.path(
    tcr_result_dir(), "shared_clones_csf_cell_count_correlation.pdf"
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = 5, height = 3)
  path
}

write_tcr_analysis_plots <- function(sc_tcr, lookup, tables, seed = 42L) {
  root <- tcr_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  save_plot <- function(name, plot, width, height) {
    path <- file.path(root, name)
    withr::with_seed(seed, ggplot2::ggsave(path, plot = plot, width = width, height = height))
    path
  }
  main <- subset(sc_tcr, subset = diagnosis %in% c("CTRL", "CIAP", "CIDP", "GBS"))
  main$tissue_diagnosis <- droplevels(main$tissue_diagnosis)
  overlay <- function(object, facet = NULL) {
    plot <- if (is.null(facet)) {
      scRepertoire::clonalOverlay(
        object, reduction = "umap.stacas.ss.all", cutpoint = 20,
        bins = 25, pt.size = 0.1, pt.alpha = 0.1
      )
    } else {
      scRepertoire::clonalOverlay(
        object, reduction = "umap.stacas.ss.all", cutpoint = 20,
        bins = 25, pt.size = 0.1, pt.alpha = 0.1, facet.by = facet
      )
    }
    plot +
      ggplot2::scale_color_manual(values = object@misc$cluster_col) +
      scMisc::theme_rect() + Seurat::NoLegend() +
      ggplot2::labs(x = "UMAP1", y = "UMAP2")
  }
  overlay_all <- overlay(sc_tcr)
  overlay_group <- overlay(sc_tcr, "tissue_group")
  overlay_dx <- overlay(main) + ggplot2::facet_wrap(~tissue_diagnosis, nrow = 2)

  overlap_dx <- scRepertoire::clonalOverlap(
    main, cloneCall = "aa", method = "overlap", group.by = "tissue_diagnosis"
  )
  overlap_sample <- scRepertoire::clonalOverlap(
    main, cloneCall = "aa", method = "overlap", group.by = "sample"
  ) + ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

  total_unique <- dplyr::n_distinct(stats::na.omit(sc_tcr$CTaa))
  summary <- tables$shared_clones_summary
  shared_plot <- tibble::tibble(
    category = c("unique", "shared_between_patients", "shared_within_patient"),
    count = c(
      total_unique - nrow(summary), sum(summary$patients_count > 1),
      sum(summary$patients_count == 1)
    )
  ) |>
    ggplot2::ggplot(ggplot2::aes(x = "", y = .data$count, fill = .data$category)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$count),
      position = ggplot2::position_stack(vjust = 0.5)
    ) +
    ggplot2::coord_polar("y") + ggplot2::theme_void()
  enrichment <- tables$shared_clones_csf_pbmc_abundance |>
    ggplot2::ggplot(
      ggplot2::aes(x = .data$diagnosis, y = .data$log2_ratio, fill = .data$diagnosis)
    ) +
    ggplot2::geom_violin(alpha = 0.5, width = 1.5) +
    ggplot2::geom_boxplot(alpha = 0.5, width = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::scale_fill_manual(values = sc_tcr@misc$diagnosis_col) +
    ggplot2::scale_y_continuous(breaks = seq(-10, 10)) +
    ggplot2::labs(x = NULL, y = "log2(CSF / blood frequency)") +
    ggplot2::theme_classic() + ggplot2::theme(legend.position = "none")

  distribution <- scRepertoire::clonalSizeDistribution(
    sc_tcr, cloneCall = "aa", method = "ward.D2", group.by = "sample"
  )
  pbmc <- subset(sc_tcr, subset = tissue == "PBMC")
  diversity_data <- scRepertoire::clonalDiversity(
    pbmc, cloneCall = "aa", group.by = "patient", x.axis = "diagnosis",
    exportTable = TRUE, metric = "ace"
  ) |>
    dplyr::rename(value = "diagnosis") |>
    dplyr::mutate(metric = "ACE") |>
    dplyr::left_join(dplyr::select(lookup, "patient", "diagnosis"), by = "patient")
  diversity <- ggplot2::ggplot(
    diversity_data,
    ggplot2::aes(x = .data$diagnosis, y = .data$value, fill = .data$diagnosis)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    ggplot2::geom_jitter(shape = 21, width = 0.2) +
    ggplot2::scale_fill_manual(values = pbmc@misc$diagnosis_col) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Diversity Score", fill = "Diagnosis") +
    ggplot2::theme_bw() + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      legend.position = "none"
    )
  Seurat::Idents(main) <- main$cluster
  startrac <- scRepertoire::StartracDiversity(
    main, cloneCall = "strict", chain = "both", type = "tissue",
    group.by = "patient"
  ) + ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

  mean_data <- sc_tcr[[]] |>
    dplyr::group_by(.data$patient, .data$tissue, .data$diagnosis) |>
    dplyr::summarise(
      clonalFrequency = mean(.data$clonalFrequency, na.rm = TRUE), .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = "tissue", values_from = "clonalFrequency",
      names_prefix = "clonalFreq_"
    ) |>
    tidyr::drop_na("clonalFreq_CSF", "clonalFreq_PBMC")
  mean_correlation <- ggplot2::ggplot(
    mean_data,
    ggplot2::aes(
      x = .data$clonalFreq_CSF, y = .data$clonalFreq_PBMC,
      color = .data$diagnosis
    )
  ) + ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = sc_tcr@misc$diagnosis_col) +
    ggplot2::labs(
      x = "Mean Clonal Frequency (CSF)", y = "Mean Clonal Frequency (PBMC)",
      title = "Mean Clonal Frequency per Patient"
    ) + ggplot2::theme_classic()
  correlation_input <- sc_tcr[[]] |>
    dplyr::filter(!is.na(.data$clonalFrequency), !is.na(.data$CTaa)) |>
    dplyr::select("patient", "tissue", "diagnosis", "clonalFrequency", "CTaa") |>
    dplyr::distinct()
  matched <- dplyr::inner_join(
    dplyr::filter(correlation_input, .data$tissue == "CSF"),
    dplyr::filter(correlation_input, .data$tissue == "PBMC"),
    by = c("patient", "diagnosis", "CTaa"), suffix = c("_CSF", "_PBMC")
  ) |>
    dplyr::distinct()
  all_correlation <- ggplot2::ggplot(
    matched,
    ggplot2::aes(
      x = .data$clonalFrequency_CSF, y = .data$clonalFrequency_PBMC,
      color = .data$diagnosis
    )
  ) + ggplot2::geom_jitter(size = 1, alpha = 0.6) +
    ggplot2::scale_color_manual(values = sc_tcr@misc$diagnosis_col) +
    ggplot2::labs(
      x = "Clonal Frequency (CSF)", y = "Clonal Frequency (PBMC)",
      title = "Matched Clones Across Tissues"
    ) + ggplot2::theme_classic()

  protein_shared <- tables$shared_clones_csf_pbmc_abundance |>
    dplyr::left_join(dplyr::select(lookup, "patient", "csf_protein"), by = "patient")
  protein_shared_plot <- ggplot2::ggplot(
    protein_shared,
    ggplot2::aes(x = .data$csf_protein, y = .data$log2_ratio, color = .data$diagnosis)
  ) + ggplot2::geom_point(size = 1) +
    ggplot2::geom_smooth(
      ggplot2::aes(group = 1), method = "lm", se = TRUE, linewidth = 0.5,
      color = "black"
    ) + ggplot2::scale_color_manual(values = sc_tcr@misc$diagnosis_col) +
    ggplot2::labs(x = "CSF Protein (mg/L)", y = "log2(CSF / blood frequency)") +
    ggplot2::theme_classic()
  protein_data <- sc_tcr[[]] |>
    dplyr::filter(.data$tissue == "CSF") |>
    dplyr::left_join(dplyr::select(lookup, "patient", "csf_protein"), by = "patient") |>
    dplyr::group_by(.data$sample) |>
    dplyr::summarise(
      clonalFrequency = mean(.data$clonalFrequency, na.rm = TRUE),
      csf_protein = mean(.data$csf_protein, na.rm = TRUE),
      diagnosis = dplyr::first(.data$diagnosis), .groups = "drop"
    )
  test <- stats::cor.test(protein_data$csf_protein, protein_data$clonalFrequency)
  protein_plot <- ggplot2::ggplot(
    protein_data,
    ggplot2::aes(
      x = .data$csf_protein, y = .data$clonalFrequency, color = .data$diagnosis
    )
  ) + ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = sc_tcr@misc$diagnosis_col) +
    ggplot2::geom_smooth(method = "lm", color = "blue") +
    ggplot2::labs(
      x = "CSF Protein (mg/L)", y = "Mean Clonal Frequency",
      subtitle = sprintf("r = %.3g, p = %.3g", test$estimate, test$p.value)
    ) + ggplot2::theme_classic()

  c(
    save_plot("tcr_clonal_overlay.pdf", overlay_all, 10, 10),
    save_plot("tcr_clonal_overlay_group.pdf", overlay_group, 5, 7),
    save_plot("tcr_clonal_overlay_diagnosis.pdf", overlay_dx, 10, 7),
    save_plot("tcr_alluvial_sample_tissue_diagnosis.pdf", make_tcr_alluvial_plot(sc_tcr), 15, 12),
    save_plot("tcr_alluvial_sample_tissue_diagnosis_main_groups.pdf", make_tcr_alluvial_plot(main), 15, 12),
    save_plot("tcr_clonal_overlap_tissue_diagnosis.pdf", overlap_dx, 15, 10),
    save_plot("tcr_clonal_overlap_tissue_sample.pdf", overlap_sample, 20, 17),
    save_plot("tcr_shared_clones_summary.pdf", shared_plot, 5, 5),
    save_plot("shared_clones_enrichment_ratio_expanded.pdf", enrichment, 5, 3),
    save_plot("tcr_clonal_distribution.pdf", distribution, 10, 7),
    save_plot("tcr_clonal_diversity.pdf", diversity, 10, 7),
    save_plot("tcr_startrac_diversity.pdf", startrac, 4, 5),
    save_plot("tcr_csf_pbmc_correlation_mean.pdf", mean_correlation, 8, 6),
    save_plot("tcr_csf_pbmc_correlation_all.pdf", all_correlation, 8, 6),
    save_plot("shared_clones_csf_protein_correlation.pdf", protein_shared_plot, 5, 3),
    save_plot("tcr_clonality_csf_protein_correlation.pdf", protein_plot, 8, 6)
  )
}

make_tcr_alluvial_plot <- function(object) {
  top_clones <- object[[]] |>
    dplyr::count(.data$CTaa, .data$tissue, .data$diagnosis) |>
    tidyr::drop_na() |>
    dplyr::arrange(dplyr::desc(.data$n)) |>
    dplyr::slice_max(.data$n, n = 10L, with_ties = FALSE)
  object$CTaa_top <- ifelse(
    object$CTaa %in% top_clones$CTaa,
    object$CTaa,
    NA
  )
  scRepertoire::alluvialClones(
    object,
    cloneCall = "aa",
    y.axes = c("tissue", "patient", "diagnosis", "cluster"),
    color = "CTaa_top"
  ) +
    ggplot2::scale_fill_manual(values = scales::hue_pal()(10L))
}

write_tcr_alluvial_plots <- function(sc_tcr, seed = 42L) {
  root <- tcr_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  main <- subset(sc_tcr, subset = diagnosis %in% c("CTRL", "CIAP", "CIDP", "GBS"))
  paths <- file.path(
    root,
    c(
      "tcr_alluvial_sample_tissue_diagnosis.pdf",
      "tcr_alluvial_sample_tissue_diagnosis_main_groups.pdf"
    )
  )
  plots <- list(make_tcr_alluvial_plot(sc_tcr), make_tcr_alluvial_plot(main))
  for (index in seq_along(paths)) {
    withr::with_seed(
      seed + index - 1L,
      ggplot2::ggsave(paths[[index]], plots[[index]], width = 15, height = 12)
    )
  }
  paths
}

write_tcr_invariant_plots <- function(sc_tcr, seed = 42L) {
  root <- tcr_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  object <- scRepertoire::annotateInvariant(sc_tcr, type = "MAIT", species = "human")
  object <- scRepertoire::annotateInvariant(object, type = "iNKT", species = "human")
  make_plot <- function(type) {
    score <- paste0(type, ".score")
    plot <- Seurat::DimPlot(
      object, reduction = "umap.stacas.ss.all", group.by = score,
      pt.size = 0.1, cols = c("0" = "grey", "1" = "red"), raster = FALSE
    ) + Seurat::NoLegend() + scMisc::theme_rect() +
      ggplot2::labs(x = "UMAP1", y = "UMAP2", title = paste(type, "Cells Highlighted"))
    plot[[1]]$layers[[1]]$aes_params$alpha <- ifelse(
      object[[score, drop = TRUE]] == 1, 1, 0.01
    )
    plot
  }
  paths <- file.path(root, c("highlighted_MAIT_cells.pdf", "highlighted_iNKT_cells.pdf"))
  withr::with_seed(seed, {
    ggplot2::ggsave(paths[[1L]], make_plot("MAIT"), width = 5, height = 5)
    ggplot2::ggsave(paths[[2L]], make_plot("iNKT"), width = 5, height = 5)
  })
  paths
}
