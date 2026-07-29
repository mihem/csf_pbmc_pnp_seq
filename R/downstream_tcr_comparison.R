tcr_comparison_input_manifest <- function(
    root = file.path("raw", "sukenikova"),
    reactive = file.path("lookup", "sukenikova_nature_gbs_supp_table_2.csv")
) {
  files <- sort(list.files(root, pattern = "\\.tsv$", full.names = TRUE))
  stopifnot(dir.exists(root), length(files) == 22L, file.exists(reactive))
  c(files, reactive)
}

read_sukenikova_bulk <- function(files) {
  purrr::map_dfr(files, function(path) {
    filename <- basename(path)
    parts <- stringr::str_match(
      filename, "^(PT\\d+)_(AC|REC)_(blood|CSF|biopsy)_(.+)\\.tsv$"
    )
    data <- janitor::clean_names(readr::read_tsv(path, show_col_types = FALSE))
    if ("frame_type" %in% names(data)) {
      data <- dplyr::filter(data, .data$frame_type == "In")
    } else {
      data <- dplyr::filter(
        data, !is.na(.data$productive_frequency), .data$productive_frequency > 0
      )
    }
    count_column <- if ("templates" %in% names(data)) "templates" else "seq_reads"
    data |>
      dplyr::mutate(
        v_raw = dplyr::coalesce(.data$v_gene, .data$v_gene_ties),
        j_raw = dplyr::coalesce(.data$j_gene, .data$j_gene_ties),
        v = stringr::str_extract(
          stringr::str_replace_all(.data$v_raw, "TCRBV", "TRBV"), "TRBV[0-9\\-]+"
        ),
        j = stringr::str_extract(
          stringr::str_replace_all(.data$j_raw, "TCRBJ", "TRBJ"), "TRBJ[0-9\\-]+"
        )
      ) |>
      dplyr::filter(
        !is.na(.data$amino_acid), .data$amino_acid != "",
        dplyr::between(nchar(.data$amino_acid), 5L, 25L),
        !stringr::str_detect(.data$amino_acid, "\\*|_|X")
      ) |>
      dplyr::group_by(.data$amino_acid, .data$v) |>
      dplyr::summarise(
        count = sum(.data[[count_column]], na.rm = TRUE),
        j = dplyr::first(stats::na.omit(.data$j), default = NA_character_),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        frequency = .data$count / sum(.data$count), file_name = filename,
        suk_patient = parts[1L, 2L], suk_timepoint = parts[1L, 3L],
        suk_tissue = parts[1L, 4L], suk_celltype = parts[1L, 5L]
      )
  })
}

prepare_tcr_comparison <- function(sc_tcr, input_files, seed = 42L) {
  stopifnot(length(input_files) == 23L, all(file.exists(input_files)))
  bulk <- read_sukenikova_bulk(input_files[grepl("\\.tsv$", input_files)])
  reactive <- readr::read_csv(
    input_files[grepl("\\.csv$", input_files)], show_col_types = FALSE
  )
  metadata <- sc_tcr[[]] |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::filter(!is.na(.data$CTaa)) |>
    dplyr::mutate(
      TRB_CDR3 = sub(".*_", "", .data$CTaa),
      TRBV = sub("\\..*", "", sub(".*_", "", .data$CTgene))
    )
  suk <- bulk |>
    dplyr::group_by(.data$suk_patient, .data$amino_acid, .data$v) |>
    dplyr::summarise(
      total_count = sum(.data$count),
      tissues = paste(sort(unique(stats::na.omit(.data$suk_tissue))), collapse = ","),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$suk_patient) |>
    dplyr::slice_max(.data$total_count, n = 1000L, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(cell_id = paste0("suk_", .data$suk_patient, "_", dplyr::row_number())) |>
    dplyr::transmute(
      cell_id = .data$cell_id, CDR3b = .data$amino_acid,
      sample = paste0("Suk_", .data$suk_patient), patient = .data$suk_patient,
      tissue = .data$tissues, diagnosis = "GBS_Sukenikova", source = "Sukenikova"
    )
  heming <- metadata |>
    dplyr::transmute(
      cell_id = .data$cell_id, CDR3b = .data$TRB_CDR3, sample = .data$sample,
      patient = .data$patient, tissue = .data$tissue,
      diagnosis = as.character(.data$diagnosis), source = "Heming"
    ) |>
    dplyr::filter(!is.na(.data$CDR3b), .data$CDR3b != "")
  combined <- dplyr::bind_rows(heming, suk)
  gliph <- withr::with_seed(seed, turboGliph::gliph_combined(
    cdr3_sequences = data.frame(CDR3b = combined$CDR3b),
    local_similarities = TRUE, local_method = "fisher",
    global_similarities = TRUE, global_method = "fisher", motif_length = 3,
    gccutoff = 1, vgene_match = "none", n_cores = 1, result_folder = "",
    lcminove = 10
  ))
  clusters <- gliph$cluster_properties
  membership <- purrr::map_dfr(seq_len(nrow(clusters)), function(index) {
    members <- strsplit(as.character(clusters$members[[index]]), "\\s+")[[1L]]
    tibble::tibble(
      cluster = index, CDR3b = members[members != ""],
      cluster_size = clusters$cluster_size[[index]]
    )
  })
  cluster_meta <- membership |>
    dplyr::left_join(dplyr::distinct(combined), by = "CDR3b") |>
    dplyr::mutate(tissue_standard = dplyr::case_when(
      .data$tissue == "CSF" ~ "CSF",
      .data$tissue == "PBMC" | stringr::str_detect(.data$tissue, "blood") ~ "PBMC",
      .data$tissue == "biopsy" ~ "Biopsy",
      TRUE ~ .data$tissue
    ))
  totals <- dplyr::count(
    dplyr::filter(cluster_meta, !is.na(.data$diagnosis)), .data$diagnosis,
    name = "total"
  )
  diagnoses <- setdiff(unique(stats::na.omit(cluster_meta$diagnosis)), "CTRL")
  diagnosis_enrichment <- purrr::map_dfr(unique(cluster_meta$cluster), function(cluster) {
    cluster_data <- dplyr::filter(cluster_meta, .data$cluster == .env$cluster)
    if (nrow(cluster_data) < 3L) return(tibble::tibble())
    purrr::map_dfr(diagnoses, function(diagnosis) {
      in_dx <- sum(cluster_data$diagnosis == diagnosis, na.rm = TRUE)
      in_ctrl <- sum(cluster_data$diagnosis == "CTRL", na.rm = TRUE)
      total_dx <- totals$total[match(diagnosis, totals$diagnosis)]
      total_ctrl <- totals$total[match("CTRL", totals$diagnosis)]
      test <- stats::fisher.test(matrix(
        c(in_dx, in_ctrl, total_dx - in_dx, total_ctrl - in_ctrl), nrow = 2L
      ))
      tibble::tibble(
        cluster = cluster, diagnosis = diagnosis, cluster_size = nrow(cluster_data),
        odds_ratio = unname(test$estimate), p_value = test$p.value
      )
    })
  }) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p_value, "BH"))
  tissue_enrichment <- purrr::map_dfr(unique(cluster_meta$cluster), function(cluster) {
    data <- dplyr::filter(cluster_meta, .data$cluster == .env$cluster)
    n_csf <- sum(data$tissue_standard == "CSF", na.rm = TRUE)
    n_pbmc <- sum(data$tissue_standard == "PBMC", na.rm = TRUE)
    if (nrow(data) < 3L || n_csf + n_pbmc < 2L) return(tibble::tibble())
    total_csf <- sum(cluster_meta$tissue_standard == "CSF", na.rm = TRUE)
    total_pbmc <- sum(cluster_meta$tissue_standard == "PBMC", na.rm = TRUE)
    test <- stats::fisher.test(matrix(
      c(n_csf, total_csf - n_csf, n_pbmc, total_pbmc - n_pbmc), nrow = 2L
    ))
    tibble::tibble(
      cluster = cluster, n_csf = n_csf, n_pbmc = n_pbmc,
      tissue_or = unname(test$estimate), tissue_p = test$p.value
    )
  }) |>
    dplyr::mutate(
      tissue_padj = stats::p.adjust(.data$tissue_p, "BH"),
      tissue_bias = dplyr::case_when(
        .data$tissue_padj < 0.1 & .data$tissue_or > 1 ~ "CSF-enriched",
        .data$tissue_padj < 0.1 & .data$tissue_or < 1 ~ "PBMC-enriched",
        TRUE ~ "No bias"
      )
    )
  myelin <- dplyr::transmute(
    reactive, CDR3b = .data$cdr3b_aa, specificity = .data$specificity,
    pns_myelin_antigen = .data$pns_myelin_antigen,
    gliph_cluster_id = .data$gliph_cluster_id,
    public_clonotype = .data$public_clonotype, pt = .data$pt
  )
  exact <- metadata |>
    dplyr::select("cell_id", "TRB_CDR3", "sample", "patient", "tissue", "diagnosis") |>
    dplyr::distinct() |>
    dplyr::inner_join(myelin, by = c("TRB_CDR3" = "CDR3b"))
  myelin_clusters <- cluster_meta$cluster[cluster_meta$CDR3b %in% myelin$CDR3b]
  cluster_matches <- cluster_meta |>
    dplyr::filter(
      .data$cluster %in% .env$myelin_clusters, .data$source == "Heming",
      !.data$CDR3b %in% myelin$CDR3b
    ) |>
    dplyr::distinct(.data$cluster, .data$CDR3b, .data$patient, .data$tissue, .data$diagnosis)
  neuropathy <- c("CIAP", "CIDP", "GBS", "MAG", "MFS", "PNC", "CAN", "PPN")
  myelin_enrichment <- purrr::map_dfr(neuropathy, function(diagnosis) {
    dx_total <- sum(metadata$diagnosis == diagnosis)
    ctrl_total <- sum(metadata$diagnosis == "CTRL")
    dx_match <- sum(exact$diagnosis == diagnosis) + sum(cluster_matches$diagnosis == diagnosis)
    ctrl_match <- sum(exact$diagnosis == "CTRL") + sum(cluster_matches$diagnosis == "CTRL")
    test <- stats::fisher.test(matrix(
      c(dx_match, ctrl_match, dx_total - dx_match, ctrl_total - ctrl_match), nrow = 2L
    ))
    tibble::tibble(
      diagnosis = diagnosis, n_myelin_reactive = dx_match, n_total = dx_total,
      fraction = dx_match / dx_total, odds_ratio = unname(test$estimate),
      p_value = test$p.value
    )
  }) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p_value, "BH"))
  list(
    diagnosis_enrichment = diagnosis_enrichment,
    tissue_enrichment = tissue_enrichment,
    myelin_enrichment = myelin_enrichment,
    exact_matches = exact,
    diagnosis_colors = c(sc_tcr@misc$diagnosis_col, GBS_Sukenikova = "#FF7F00")
  )
}

tcr_comparison_theme <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) + ggplot2::theme(
    plot.title = ggplot2::element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    strip.text = ggplot2::element_text(size = base_size, face = "bold"),
    strip.background = ggplot2::element_rect(fill = "grey95")
  )
}

write_tcr_comparison_plots <- function(result, seed = 42L) {
  root <- file.path(tcr_result_dir(), "comparison")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  save <- function(name, plot, width, height) {
    path <- file.path(root, name)
    withr::with_seed(seed, ggplot2::ggsave(path, plot, width = width, height = height))
    path
  }
  overview <- result$diagnosis_enrichment |>
    dplyr::filter(.data$p_adj < 0.1, .data$odds_ratio > 1) |>
    dplyr::left_join(
      dplyr::select(result$tissue_enrichment, "cluster", "tissue_bias"), by = "cluster"
    ) |>
    dplyr::mutate(
      tissue_bias = tidyr::replace_na(.data$tissue_bias, "No bias"),
      neg_log10_fdr = -log10(pmax(.data$p_adj, 1e-300))
    )
  labels <- overview |>
    dplyr::group_by(.data$diagnosis) |>
    dplyr::slice_max(.data$cluster_size, n = 2L, with_ties = FALSE) |>
    dplyr::ungroup()
  overview_plot <- ggplot2::ggplot(
    overview, ggplot2::aes(x = .data$diagnosis, y = .data$cluster_size)
  ) +
    ggplot2::geom_jitter(
      ggplot2::aes(
        fill = .data$diagnosis, size = .data$neg_log10_fdr,
        shape = .data$tissue_bias
      ),
      width = 0.2, alpha = 0.75, stroke = 0.4
    ) +
    ggplot2::scale_fill_manual(values = result$diagnosis_colors, guide = "none") +
    ggplot2::scale_shape_manual(values = c(
      "CSF-enriched" = 24, "No bias" = 21, "PBMC-enriched" = 25
    )) +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "-log10(FDR)") +
    ggrepel::geom_text_repel(
      data = labels, ggplot2::aes(label = paste0("C", .data$cluster)),
      size = 2.5, max.overlaps = 15, seed = seed
    ) + tcr_comparison_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "Disease-enriched GLIPH2 specificity groups overview", x = NULL,
      y = "Cluster size (# clonotypes)"
    )
  tissue <- overview |>
    dplyr::count(.data$diagnosis, .data$tissue_bias, name = "n_clusters") |>
    dplyr::mutate(tissue_bias = factor(
      .data$tissue_bias, c("CSF-enriched", "No bias", "PBMC-enriched")
    ))
  tissue_plot <- ggplot2::ggplot(
    tissue,
    ggplot2::aes(x = .data$diagnosis, y = .data$n_clusters, fill = .data$tissue_bias)
  ) + ggplot2::geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c(
      "CSF-enriched" = "#E41A1C", "No bias" = "grey70", "PBMC-enriched" = "#377EB8"
    )) + tcr_comparison_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "Tissue compartmentalization of disease-enriched clusters",
      x = "Diagnosis", y = "Number of enriched clusters"
    )
  myelin <- dplyr::mutate(
    result$myelin_enrichment,
    label = ifelse(.data$p_adj < 0.05, "*", ""),
    diagnosis = factor(.data$diagnosis, levels = unique(.data$diagnosis))
  )
  myelin_plot <- ggplot2::ggplot(
    myelin,
    ggplot2::aes(x = .data$diagnosis, y = .data$fraction * 100, fill = .data$diagnosis)
  ) + ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(.data$label, "\n", sprintf("OR=%.1f", .data$odds_ratio))),
      vjust = -0.3, size = 2.8
    ) + ggplot2::scale_fill_manual(values = result$diagnosis_colors) +
    tcr_comparison_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "Myelin-reactive TCR fraction by diagnosis",
      subtitle = "Exact + cluster-level matches to Sukenikova myelin-reactive clones",
      x = NULL, y = "Myelin-reactive TCRs (%)"
    ) + ggplot2::expand_limits(y = max(myelin$fraction * 100) * 1.4)
  table_data <- result$exact_matches |>
    dplyr::arrange(.data$tissue, .data$diagnosis, .data$TRB_CDR3) |>
    dplyr::select("TRB_CDR3", "patient", "tissue", "diagnosis", "pt", "specificity")
  table_path <- file.path(root, "sukenikova_reactive_table.pdf")
  grDevices::pdf(table_path, width = 12, height = 8)
  grid::grid.draw(gridExtra::tableGrob(table_data, rows = NULL))
  grDevices::dev.off()
  c(
    save("fig_gliph_enriched_clusters_overview_bubble.pdf", overview_plot, 9, 6.5),
    save("fig_gliph_tissue_bias_fraction_by_diagnosis.pdf", tissue_plot, 8, 5.5),
    save("fig_gliph_myelin_reactive_fraction_by_diagnosis.pdf", myelin_plot, 7, 5.5),
    table_path
  )
}
