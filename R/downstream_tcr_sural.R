sural_tcr_comparison_result_dir <- function() {
  file.path(sural_trust4_result_dir(), "cross_tissue_tcr")
}

valid_trb_cdr3 <- function(sequence) {
  !is.na(sequence) &
    dplyr::between(nchar(sequence), 5L, 30L) &
    grepl("^C[A-Z]+[FW]$", sequence) &
    !grepl("[X*_]", sequence)
}

extract_csf_pbmc_trb <- function(tcr_contigs, patients) {
  selected <- tcr_contigs[
    sub(".*_", "", names(tcr_contigs)) %in% patients
  ]
  stopifnot(
    length(selected) == length(patients) * 2L,
    all(grepl("^(CSF|PBMC)_P[0-9]+$", names(selected)))
  )

  purrr::imap_dfr(selected, function(contigs, sample_id) {
    required <- c(
      "barcode", "is_cell", "high_confidence", "chain", "v_gene",
      "d_gene", "j_gene", "c_gene", "productive", "cdr3", "umis",
      "reads"
    )
    stopifnot(all(required %in% names(contigs)))
    is_true <- function(value) {
      as.character(value) %in% c("TRUE", "True", "true")
    }

    contigs |>
      dplyr::filter(
        .data$chain == "TRB",
        is_true(.data$is_cell),
        is_true(.data$high_confidence),
        is_true(.data$productive)
      ) |>
      dplyr::transmute(
        patient = sub(".*_", "", sample_id),
        tissue = sub("_.*", "", sample_id),
        sample = sample_id,
        cell_id = paste(sample_id, .data$barcode, sep = "_"),
        TRB_CDR3aa = toupper(trimws(.data$cdr3)),
        V = .data$v_gene,
        D = .data$d_gene,
        J = .data$j_gene,
        C = .data$c_gene,
        umi_count = as.numeric(.data$umis),
        read_count = as.numeric(.data$reads),
        source = "10x VDJ",
        immune_subcluster = NA_character_
      ) |>
      dplyr::filter(valid_trb_cdr3(.data$TRB_CDR3aa)) |>
      dplyr::distinct(
        .data$patient, .data$tissue, .data$cell_id, .data$TRB_CDR3aa,
        .keep_all = TRUE
      )
  })
}

extract_sural_trb <- function(
  trust4_records, patient_map, primary_only = FALSE
) {
  patient_map <- unlist(patient_map, use.names = TRUE)
  chain_columns <- c(
    "chain1", "chain2", "secondary_chain1", "secondary_chain2"
  )
  stopifnot(
    all(c("library_id", "sample", "barcode", chain_columns) %in%
      names(trust4_records)),
    setequal(unique(trust4_records$library_id), names(patient_map))
  )

  trust4_records |>
    dplyr::filter(.data$receptor == "TCR") |>
    tidyr::pivot_longer(
      dplyr::all_of(chain_columns),
      names_to = "chain_slot",
      values_to = "chain"
    ) |>
    dplyr::filter(
      !primary_only | .data$chain_slot %in% c("chain1", "chain2")
    ) |>
    tidyr::separate_longer_delim("chain", delim = ";") |>
    dplyr::filter(!is.na(.data$chain), .data$chain != "*") |>
    tidyr::separate_wider_delim(
      "chain",
      delim = ",",
      names = c(
        "V", "D", "J", "C", "cdr3_nt", "TRB_CDR3aa",
        "read_count", "consensus_id", "similarity", "complete"
      ),
      too_few = "align_start",
      too_many = "merge"
    ) |>
    dplyr::mutate(
      patient = unname(patient_map[.data$library_id]),
      tissue = "Sural",
      sample = .data$library_id,
      cell_id = paste(.data$library_id, .data$barcode, sep = "_"),
      TRB_CDR3aa = toupper(trimws(.data$TRB_CDR3aa)),
      read_count = suppressWarnings(as.numeric(.data$read_count)),
      umi_count = NA_real_,
      source = "TRUST4",
      immune_subcluster = .data$seurat_cluster
    ) |>
    dplyr::filter(
      grepl("^TRBV", .data$V),
      valid_trb_cdr3(.data$TRB_CDR3aa)
    ) |>
    dplyr::select(
      "patient", "tissue", "sample", "cell_id", "TRB_CDR3aa",
      "V", "D", "J", "C", "umi_count", "read_count", "source",
      "library_id", "barcode", "match_status", "immune_subcluster"
    ) |>
    dplyr::distinct(
      .data$patient, .data$tissue, .data$cell_id, .data$TRB_CDR3aa,
      .keep_all = TRUE
    )
}

prepare_sural_tcr_comparison <- function(
  tcr_contigs, trust4_records, patient_map
) {
  patient_map <- unlist(patient_map, use.names = TRUE)
  patients <- unname(patient_map)
  stopifnot(!anyDuplicated(patients), length(patients) == 3L)

  csf_pbmc <- extract_csf_pbmc_trb(tcr_contigs, patients)
  sural <- extract_sural_trb(trust4_records, patient_map)
  cell_clones <- dplyr::bind_rows(csf_pbmc, sural) |>
    dplyr::mutate(
      patient = factor(.data$patient, levels = patients),
      tissue = factor(.data$tissue, levels = c("CSF", "PBMC", "Sural"))
    ) |>
    dplyr::arrange(.data$patient, .data$tissue, .data$TRB_CDR3aa)

  clone_counts <- cell_clones |>
    dplyr::group_by(.data$patient, .data$tissue, .data$TRB_CDR3aa) |>
    dplyr::summarise(
      cell_count = dplyr::n_distinct(.data$cell_id),
      read_count = sum(.data$read_count, na.rm = TRUE),
      V = paste(sort(unique(stats::na.omit(.data$V))), collapse = ";"),
      J = paste(sort(unique(stats::na.omit(.data$J))), collapse = ";"),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$patient, .data$tissue) |>
    dplyr::mutate(frequency = .data$cell_count / sum(.data$cell_count)) |>
    dplyr::ungroup()

  tracking <- clone_counts |>
    dplyr::select(
      "patient", "TRB_CDR3aa", "tissue", "cell_count", "frequency"
    ) |>
    tidyr::pivot_wider(
      names_from = "tissue",
      values_from = c("cell_count", "frequency"),
      names_glue = "{.value}_{tissue}",
      values_fill = 0
    ) |>
    dplyr::mutate(
      n_tissues =
        (.data$cell_count_CSF > 0) +
        (.data$cell_count_PBMC > 0) +
        (.data$cell_count_Sural > 0),
      presence = dplyr::case_when(
        .data$cell_count_CSF > 0 & .data$cell_count_PBMC > 0 &
          .data$cell_count_Sural > 0 ~ "CSF + blood + sural",
        .data$cell_count_CSF > 0 & .data$cell_count_PBMC > 0 ~
          "CSF + blood",
        .data$cell_count_CSF > 0 & .data$cell_count_Sural > 0 ~
          "CSF + sural",
        .data$cell_count_PBMC > 0 & .data$cell_count_Sural > 0 ~
          "blood + sural",
        .data$cell_count_CSF > 0 ~ "CSF only",
        .data$cell_count_PBMC > 0 ~ "blood only",
        TRUE ~ "sural only"
      )
    ) |>
    dplyr::arrange(
      .data$patient, dplyr::desc(.data$n_tissues), .data$TRB_CDR3aa
    )

  presence_levels <- c(
    "CSF only", "blood only", "sural only", "CSF + blood",
    "CSF + sural", "blood + sural", "CSF + blood + sural"
  )
  intersection_summary <- tracking |>
    dplyr::count(.data$patient, .data$presence, name = "clonotype_count") |>
    tidyr::complete(
      patient = factor(patients, levels = patients),
      presence = presence_levels,
      fill = list(clonotype_count = 0L)
    ) |>
    dplyr::mutate(
      presence = factor(.data$presence, levels = presence_levels)
    ) |>
    dplyr::arrange(.data$patient, .data$presence)

  sural_shared <- tracking |>
    dplyr::filter(
      .data$cell_count_Sural > 0,
      .data$cell_count_CSF > 0 | .data$cell_count_PBMC > 0
    )
  shared_cell_details <- cell_clones |>
    dplyr::semi_join(
      sural_shared,
      by = c("patient", "TRB_CDR3aa")
    )

  list(
    patient_mapping = tibble::tibble(
      sural_library = names(patient_map), patient = patients
    ),
    cell_clones = cell_clones,
    clone_counts = clone_counts,
    tracking = tracking,
    sural_shared = sural_shared,
    shared_cell_details = shared_cell_details,
    intersection_summary = intersection_summary
  )
}

write_sural_tcr_comparison_workbook <- function(comparison) {
  path <- file.path(
    sural_tcr_comparison_result_dir(),
    "tcr_tracking_csf_blood_sural.xlsx"
  )
  ensure_parent_dir(path)
  writexl::write_xlsx(
    list(
      patient_mapping = comparison$patient_mapping,
      clone_counts = comparison$clone_counts,
      clone_tracking = comparison$tracking,
      sural_shared_clones = comparison$sural_shared,
      shared_cell_details = comparison$shared_cell_details,
      intersection_summary = comparison$intersection_summary
    ),
    path
  )
  path
}

write_sural_tcr_intersection_plot <- function(comparison) {
  path <- file.path(
    sural_tcr_comparison_result_dir(),
    "trb_shared_tissue_intersections.pdf"
  )
  ensure_parent_dir(path)
  data <- comparison$intersection_summary |>
    dplyr::filter(grepl("\\+", .data$presence))
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$presence,
      y = .data$clonotype_count,
      fill = .data$presence
    )
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$clonotype_count),
      vjust = -0.3,
      size = 3.2
    ) +
    ggplot2::facet_wrap(~patient, nrow = 1L) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_pretty(),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "CSF + blood" = "#6C5CE7",
        "CSF + sural" = "#D1495B",
        "blood + sural" = "#00798C",
        "CSF + blood + sural" = "#E09F3E"
      ),
      guide = "none"
    ) +
    ggplot2::labs(
      title = "Exact TRB clonotypes shared across tissues",
      subtitle = "Exact CDR3 amino-acid matches within each patient",
      x = NULL,
      y = "Number of clonotypes"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(path, plot, width = 10, height = 4.5)
  path
}

write_sural_tcr_tracking_plot <- function(comparison) {
  path <- file.path(
    sural_tcr_comparison_result_dir(),
    "sural_shared_trb_clonotype_tracking.pdf"
  )
  ensure_parent_dir(path)
  data <- comparison$sural_shared |>
    dplyr::select(
      "patient", "TRB_CDR3aa",
      dplyr::starts_with("cell_count_")
    ) |>
    tidyr::pivot_longer(
      dplyr::starts_with("cell_count_"),
      names_to = "tissue",
      values_to = "cell_count",
      names_prefix = "cell_count_"
    ) |>
    dplyr::mutate(
      tissue = factor(.data$tissue, levels = c("CSF", "PBMC", "Sural")),
      patient = droplevels(.data$patient)
    )
  stopifnot(nrow(data) > 0L)
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$tissue, y = .data$TRB_CDR3aa)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$cell_count,
        fill = .data$tissue,
        alpha = .data$cell_count > 0
      ),
      shape = 21,
      color = "black",
      stroke = 0.25
    ) +
    ggplot2::facet_wrap(~patient, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = c(CSF = "#D1495B", PBMC = "#3264A8", Sural = "#E09F3E"),
      guide = "none"
    ) +
    ggplot2::scale_alpha_manual(
      values = c(`TRUE` = 0.9, `FALSE` = 0.08),
      guide = "none"
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 7)) +
    ggplot2::labs(
      title = "Sural TRB clonotypes tracked across tissues",
      subtitle = "P28 has no exact sural-shared TRB clonotypes",
      x = NULL,
      y = "TRB CDR3 amino-acid sequence",
      size = "Cells"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  ggplot2::ggsave(path, plot, width = 8, height = 5.5)
  path
}

extract_csf_pbmc_tra <- function(tcr_contigs, patients) {
  selected <- tcr_contigs[
    sub(".*_", "", names(tcr_contigs)) %in% patients
  ]
  stopifnot(length(selected) == length(patients) * 2L)

  purrr::imap_dfr(selected, function(contigs, sample_id) {
    is_true <- function(value) {
      as.character(value) %in% c("TRUE", "True", "true")
    }
    contigs |>
      dplyr::filter(
        .data$chain == "TRA",
        is_true(.data$is_cell),
        is_true(.data$high_confidence),
        is_true(.data$productive)
      ) |>
      dplyr::transmute(
        patient = sub(".*_", "", sample_id),
        tissue = sub("_.*", "", sample_id),
        sample = sample_id,
        cell_id = paste(sample_id, .data$barcode, sep = "_"),
        clonotype = toupper(trimws(.data$cdr3)),
        source = "10x VDJ"
      ) |>
      dplyr::filter(valid_trb_cdr3(.data$clonotype)) |>
      dplyr::distinct(
        .data$patient, .data$tissue, .data$cell_id, .data$clonotype,
        .keep_all = TRUE
      )
  })
}

extract_sural_tra <- function(
  trust4_records, patient_map, primary_only = FALSE
) {
  patient_map <- unlist(patient_map, use.names = TRUE)
  chain_columns <- c(
    "chain1", "chain2", "secondary_chain1", "secondary_chain2"
  )
  trust4_records |>
    dplyr::filter(.data$receptor == "TCR") |>
    tidyr::pivot_longer(
      dplyr::all_of(chain_columns),
      names_to = "chain_slot",
      values_to = "chain"
    ) |>
    dplyr::filter(
      !primary_only | .data$chain_slot %in% c("chain1", "chain2")
    ) |>
    tidyr::separate_longer_delim("chain", delim = ";") |>
    dplyr::filter(!is.na(.data$chain), .data$chain != "*") |>
    tidyr::separate_wider_delim(
      "chain",
      delim = ",",
      names = c(
        "V", "D", "J", "C", "cdr3_nt", "clonotype",
        "read_count", "consensus_id", "similarity", "complete"
      ),
      too_few = "align_start",
      too_many = "merge"
    ) |>
    dplyr::mutate(
      patient = unname(patient_map[.data$library_id]),
      tissue = "Sural",
      sample = .data$library_id,
      cell_id = paste(.data$library_id, .data$barcode, sep = "_"),
      clonotype = toupper(trimws(.data$clonotype)),
      source = "TRUST4"
    ) |>
    dplyr::filter(
      grepl("^TRAV", .data$V),
      valid_trb_cdr3(.data$clonotype)
    ) |>
    dplyr::select(
      "patient", "tissue", "sample", "cell_id", "clonotype",
      "source", "library_id", "barcode"
    ) |>
    dplyr::distinct(
      .data$patient, .data$tissue, .data$cell_id, .data$clonotype,
      .keep_all = TRUE
    )
}

extract_csf_pbmc_paired_tcr <- function(combined_tcr, patients) {
  selected <- combined_tcr[
    sub(".*_", "", names(combined_tcr)) %in% patients
  ]
  stopifnot(length(selected) == length(patients) * 2L)

  purrr::imap_dfr(selected, function(clones, sample_id) {
    clones |>
      dplyr::filter(
        grepl("^TRAV", .data$TCR1),
        grepl("^TRBV", .data$TCR2)
      ) |>
      dplyr::transmute(
        patient = sub(".*_", "", sample_id),
        tissue = sub("_.*", "", sample_id),
        sample = sample_id,
        cell_id = .data$barcode,
        TRA_CDR3aa = toupper(trimws(.data$cdr3_aa1)),
        TRB_CDR3aa = toupper(trimws(.data$cdr3_aa2)),
        clonotype = paste(.data$TRA_CDR3aa, .data$TRB_CDR3aa, sep = " / "),
        source = "10x VDJ"
      ) |>
      dplyr::filter(
        valid_trb_cdr3(.data$TRA_CDR3aa),
        valid_trb_cdr3(.data$TRB_CDR3aa)
      ) |>
      dplyr::distinct(
        .data$patient, .data$tissue, .data$cell_id, .data$clonotype,
        .keep_all = TRUE
      )
  })
}

extract_sural_paired_tcr <- function(trust4_records, patient_map) {
  alpha <- extract_sural_tra(
    trust4_records, patient_map, primary_only = TRUE
  ) |>
    dplyr::rename(TRA_CDR3aa = "clonotype")
  beta <- extract_sural_trb(
    trust4_records, patient_map, primary_only = TRUE
  ) |>
    dplyr::select(
      "patient", "tissue", "sample", "cell_id", "library_id",
      "barcode", "TRB_CDR3aa"
    )

  alpha |>
    dplyr::inner_join(
      beta,
      by = c(
        "patient", "tissue", "sample", "cell_id", "library_id",
        "barcode"
      ),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      clonotype = paste(
        .data$TRA_CDR3aa, .data$TRB_CDR3aa, sep = " / "
      )
    ) |>
    dplyr::distinct(
      .data$patient, .data$tissue, .data$cell_id, .data$clonotype,
      .keep_all = TRUE
    )
}

build_chain_tissue_tracking <- function(cell_clones, patients) {
  tissue_counts <- cell_clones |>
    dplyr::mutate(
      patient = factor(.data$patient, levels = patients),
      tissue = factor(.data$tissue, levels = c("CSF", "PBMC", "Sural"))
    ) |>
    dplyr::group_by(.data$patient, .data$tissue, .data$clonotype) |>
    dplyr::summarise(
      cell_count = dplyr::n_distinct(.data$cell_id),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$patient, .data$tissue) |>
    dplyr::mutate(frequency = .data$cell_count / sum(.data$cell_count)) |>
    dplyr::ungroup()

  tracking <- tissue_counts |>
    tidyr::pivot_wider(
      names_from = "tissue",
      values_from = c("cell_count", "frequency"),
      names_glue = "{.value}_{tissue}",
      values_fill = 0
    ) |>
    dplyr::mutate(
      n_tissues =
        (.data$cell_count_CSF > 0) +
        (.data$cell_count_PBMC > 0) +
        (.data$cell_count_Sural > 0),
      presence = dplyr::case_when(
        .data$cell_count_CSF > 0 & .data$cell_count_PBMC > 0 &
          .data$cell_count_Sural > 0 ~ "CSF + blood + sural",
        .data$cell_count_CSF > 0 & .data$cell_count_PBMC > 0 ~
          "CSF + blood",
        .data$cell_count_CSF > 0 & .data$cell_count_Sural > 0 ~
          "CSF + sural",
        .data$cell_count_PBMC > 0 & .data$cell_count_Sural > 0 ~
          "blood + sural",
        .data$cell_count_CSF > 0 ~ "CSF only",
        .data$cell_count_PBMC > 0 ~ "blood only",
        TRUE ~ "sural only"
      )
    )

  presence_levels <- c(
    "CSF only", "blood only", "sural only", "CSF + blood",
    "CSF + sural", "blood + sural", "CSF + blood + sural"
  )
  intersections <- tracking |>
    dplyr::count(.data$patient, .data$presence, name = "clonotype_count") |>
    tidyr::complete(
      patient = factor(patients, levels = patients),
      presence = presence_levels,
      fill = list(clonotype_count = 0L)
    ) |>
    dplyr::mutate(
      presence = factor(.data$presence, levels = presence_levels)
    ) |>
    dplyr::arrange(.data$patient, .data$presence)
  sural_shared <- tracking |>
    dplyr::filter(
      .data$cell_count_Sural > 0,
      .data$cell_count_CSF > 0 | .data$cell_count_PBMC > 0
    ) |>
    dplyr::arrange(.data$patient, .data$clonotype)

  list(
    cell_clones = cell_clones,
    tissue_counts = tissue_counts,
    tracking = tracking,
    intersections = intersections,
    sural_shared = sural_shared
  )
}

prepare_sural_tcr_chain_comparisons <- function(
  tcr_contigs, combined_tcr, trust4_records, patient_map
) {
  patient_map <- unlist(patient_map, use.names = TRUE)
  patients <- unname(patient_map)
  alpha <- dplyr::bind_rows(
    extract_csf_pbmc_tra(tcr_contigs, patients),
    extract_sural_tra(trust4_records, patient_map)
  )
  paired <- dplyr::bind_rows(
    extract_csf_pbmc_paired_tcr(combined_tcr, patients),
    extract_sural_paired_tcr(trust4_records, patient_map)
  )

  list(
    alpha = build_chain_tissue_tracking(alpha, patients),
    paired = build_chain_tissue_tracking(paired, patients)
  )
}

write_sural_chain_comparison_workbook <- function(comparisons) {
  path <- file.path(
    sural_tcr_comparison_result_dir(),
    "tcr_alpha_and_paired_chain_tracking.xlsx"
  )
  ensure_parent_dir(path)
  writexl::write_xlsx(
    list(
      alpha_clone_counts = comparisons$alpha$tissue_counts,
      alpha_sural_shared = comparisons$alpha$sural_shared,
      alpha_intersections = comparisons$alpha$intersections,
      paired_clone_counts = comparisons$paired$tissue_counts,
      paired_sural_shared = comparisons$paired$sural_shared,
      paired_intersections = comparisons$paired$intersections
    ),
    path
  )
  path
}

write_sural_chain_intersection_plot <- function(
  comparison, chain_label, file_name
) {
  path <- file.path(sural_tcr_comparison_result_dir(), file_name)
  ensure_parent_dir(path)
  data <- comparison$intersections |>
    dplyr::filter(grepl("\\+", .data$presence))
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$presence,
      y = .data$clonotype_count,
      fill = .data$presence
    )
  ) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$clonotype_count),
      vjust = -0.3,
      size = 3.2
    ) +
    ggplot2::facet_wrap(~patient, nrow = 1L) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_pretty(),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "CSF + blood" = "#6C5CE7",
        "CSF + sural" = "#D1495B",
        "blood + sural" = "#00798C",
        "CSF + blood + sural" = "#E09F3E"
      ),
      guide = "none"
    ) +
    ggplot2::labs(
      title = paste("Exact", chain_label, "clonotypes shared across tissues"),
      subtitle = "Exact CDR3 amino-acid matches within each patient",
      x = NULL,
      y = "Number of clonotypes"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(path, plot, width = 10, height = 4.5)
  path
}

write_sural_chain_tracking_plot <- function(
  comparison, chain_label, file_name
) {
  path <- file.path(sural_tcr_comparison_result_dir(), file_name)
  ensure_parent_dir(path)
  data <- comparison$sural_shared |>
    dplyr::select(
      "patient", "clonotype", dplyr::starts_with("cell_count_")
    ) |>
    tidyr::pivot_longer(
      dplyr::starts_with("cell_count_"),
      names_to = "tissue",
      values_to = "cell_count",
      names_prefix = "cell_count_"
    ) |>
    dplyr::mutate(
      tissue = factor(.data$tissue, levels = c("CSF", "PBMC", "Sural")),
      patient = droplevels(.data$patient)
    )
  if (nrow(data) == 0L) {
    plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = paste0(
          "No exact sural-shared ", chain_label,
          " clonotypes were detected"
        ),
        size = 5
      ) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::labs(
        title = paste("Sural", chain_label, "clonotypes tracked across tissues")
      ) +
      ggplot2::theme_void()
    ggplot2::ggsave(path, plot, width = 9, height = 4)
    return(path)
  }
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$tissue, y = .data$clonotype)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$cell_count,
        fill = .data$tissue,
        alpha = .data$cell_count > 0
      ),
      shape = 21,
      color = "black",
      stroke = 0.25
    ) +
    ggplot2::facet_wrap(~patient, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = c(CSF = "#D1495B", PBMC = "#3264A8", Sural = "#E09F3E"),
      guide = "none"
    ) +
    ggplot2::scale_alpha_manual(
      values = c(`TRUE` = 0.9, `FALSE` = 0.08),
      guide = "none"
    ) +
    ggplot2::scale_size_continuous(range = c(1.5, 7)) +
    ggplot2::labs(
      title = paste("Sural", chain_label, "clonotypes tracked across tissues"),
      x = NULL,
      y = paste(chain_label, "CDR3 amino-acid sequence"),
      size = "Cells"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(path, plot, width = 9, height = 6)
  path
}

prepare_sural_tcr_expansion <- function(beta, chain_comparisons) {
  classify_tracking <- function(tracking, key_column) {
    tracking |>
      dplyr::transmute(
        patient = .data$patient,
        clonotype = .data[[key_column]],
        sharing_status = dplyr::case_when(
          .data$cell_count_Sural == 0 ~ "Not detected in sural",
          .data$cell_count_CSF > 0 & .data$cell_count_PBMC > 0 ~
            "All three",
          .data$cell_count_CSF > 0 ~ "Sural + CSF",
          .data$cell_count_PBMC > 0 ~ "Sural + blood",
          TRUE ~ "Sural only"
        )
      )
  }
  add_status <- function(counts, tracking, chain_type, key_column) {
    status <- classify_tracking(tracking, key_column)
    counts |>
      dplyr::rename(clonotype = dplyr::all_of(key_column)) |>
      dplyr::inner_join(
        status,
        by = c("patient", "clonotype"),
        relationship = "many-to-one"
      ) |>
      dplyr::transmute(
        patient = as.character(.data$patient),
        tissue = as.character(.data$tissue),
        chain_type = chain_type,
        clonotype = .data$clonotype,
        cell_count = .data$cell_count,
        frequency_percent = 100 * .data$frequency,
        sharing_status = .data$sharing_status,
        expanded = .data$cell_count >= 2L
      ) |>
      dplyr::filter(.data$tissue %in% c("CSF", "PBMC"))
  }

  clone_data <- dplyr::bind_rows(
    add_status(
      beta$clone_counts,
      beta$tracking,
      "TRB",
      "TRB_CDR3aa"
    ),
    add_status(
      chain_comparisons$alpha$tissue_counts,
      chain_comparisons$alpha$tracking,
      "TRA",
      "clonotype"
    ),
    add_status(
      chain_comparisons$paired$tissue_counts,
      chain_comparisons$paired$tracking,
      "TRA + TRB",
      "clonotype"
    )
  ) |>
    dplyr::mutate(
      chain_type = factor(.data$chain_type, levels = c("TRA", "TRB", "TRA + TRB")),
      tissue = factor(.data$tissue, levels = c("CSF", "PBMC")),
      sharing_status = factor(
        .data$sharing_status,
        levels = c(
          "Not detected in sural", "Sural + CSF", "Sural + blood",
          "All three"
        )
      )
    )

  summarize_expansion <- function(data, include_patient) {
    groups <- c("chain_type", "tissue", "sharing_status")
    if (include_patient) groups <- c("patient", groups)
    data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
      dplyr::summarise(
        clonotype_count = dplyr::n(),
        expanded_clonotypes = sum(.data$expanded),
        expanded_percent = 100 * mean(.data$expanded),
        median_cell_count = stats::median(.data$cell_count),
        maximum_cell_count = max(.data$cell_count),
        median_frequency_percent = stats::median(.data$frequency_percent),
        maximum_frequency_percent = max(.data$frequency_percent),
        .groups = "drop"
      )
  }

  list(
    clone_data = clone_data,
    pooled_summary = summarize_expansion(clone_data, FALSE),
    patient_summary = summarize_expansion(clone_data, TRUE)
  )
}

write_sural_tcr_expansion_workbook <- function(expansion) {
  path <- file.path(
    sural_tcr_comparison_result_dir(),
    "tcr_10x_expansion_by_sural_sharing.xlsx"
  )
  ensure_parent_dir(path)
  writexl::write_xlsx(
    list(
      pooled_summary = expansion$pooled_summary,
      patient_summary = expansion$patient_summary,
      clonotype_counts = expansion$clone_data
    ),
    path
  )
  path
}

write_sural_tcr_expansion_plot <- function(expansion) {
  path <- file.path(
    sural_tcr_comparison_result_dir(),
    "tcr_10x_expansion_by_sural_sharing.pdf"
  )
  ensure_parent_dir(path)
  data <- expansion$clone_data
  shared <- dplyr::filter(
    data, .data$sharing_status != "Not detected in sural"
  )
  colors <- c(
    "Not detected in sural" = "#A7A9AC",
    "Sural + CSF" = "#D1495B",
    "Sural + blood" = "#00798C",
    "All three" = "#E09F3E"
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$sharing_status,
      y = .data$frequency_percent,
      fill = .data$sharing_status
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.68,
      outlier.shape = NA,
      alpha = 0.7
    ) +
    ggplot2::geom_jitter(
      data = shared,
      ggplot2::aes(shape = .data$patient),
      width = 0.12,
      height = 0,
      size = 2,
      alpha = 0.85,
      color = "black"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$chain_type),
      cols = ggplot2::vars(.data$tissue),
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::scale_y_log10(
      labels = scales::label_number(accuracy = 0.001, suffix = "%")
    ) +
    ggplot2::labs(
      title = "10x TCR expansion by sural-nerve detection status",
      subtitle = paste(
        "Expansion is measured only from within-patient 10x VDJ frequencies;",
        "TRUST4 defines sural detection status only"
      ),
      x = NULL,
      y = "10x clonotype frequency",
      shape = "Patient"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(path, plot, width = 12, height = 8)
  path
}
