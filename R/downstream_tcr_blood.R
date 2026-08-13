tcr_blood_result_dir <- function() {
  file.path(tcr_result_dir(), "blood")
}

#' Put diagnoses in the study's own order rather than alphabetical.
#'
#' These tables carry `diagnosis` as character, which ggplot would otherwise
#' sort alphabetically and put out of step with every other figure in the paper.
#' The names of the palette are already in `diagnosis_order`.
order_by_diagnosis <- function(data, colors) {
  levels <- intersect(names(colors), unique(as.character(data$diagnosis)))
  dplyr::mutate(
    data,
    diagnosis = factor(as.character(.data$diagnosis), levels = levels)
  )
}

#' Clonal expansion in blood, per patient.
#'
#' Answers the part of Reviewer 4's comment 5 that asks whether clones arising
#' from a systemic trigger are detectable in blood at all. Figure 5C cannot
#' address this, because it is computed only over clones already found in both
#' compartments.
#'
#' The depth guard is the whole analysis. A patient with fewer than roughly
#' `top_n` times a few cells reaches a cumulative frequency of 1 by
#' construction, since the top clones then cover the entire repertoire. In this
#' cohort three of six GBS patients fall below any sensible floor, and whether
#' they are included decides the direction of the result: with them GBS is the
#' most clonally expanded group, without them the least. Patients are therefore
#' flagged rather than silently dropped, so the figure and the table both show
#' who is evaluable.
#'
#' @param sc_tcr Seurat object carrying scRepertoire metadata.
#' @param top_n Number of top clones whose frequencies are summed.
#' @param min_cells Blood cells with a recovered TCR required for a patient to
#'   be treated as evaluable.
#' @return List of per-clone counts, per-patient summaries, the repertoire
#'   table, and the excluded patients.
blood_clonal_expansion <- function(sc_tcr, top_n = 10L, min_cells = 50L) {
  required <- c("CTaa", "patient", "tissue", "diagnosis")
  stopifnot(
    inherits(sc_tcr, "Seurat"),
    all(required %in% colnames(sc_tcr[[]])),
    top_n >= 1L,
    min_cells >= 1L
  )

  blood <- sc_tcr[[]] |>
    dplyr::filter(
      !is.na(.data$CTaa),
      as.character(.data$tissue) == "PBMC"
    ) |>
    dplyr::mutate(diagnosis = as.character(.data$diagnosis))
  stopifnot(nrow(blood) > 0L)

  blood_depth <- dplyr::count(blood, .data$patient, name = "blood_cells")
  eligible <- blood_depth$patient[blood_depth$blood_cells >= min_cells]

  clone_counts <- blood |>
    dplyr::count(.data$patient, .data$diagnosis, .data$CTaa, name = "n_cells") |>
    dplyr::group_by(.data$patient) |>
    dplyr::mutate(
      frequency = .data$n_cells / sum(.data$n_cells),
      rank = dplyr::row_number(dplyr::desc(.data$n_cells))
    ) |>
    dplyr::ungroup() |>
    dplyr::left_join(blood_depth, by = "patient") |>
    dplyr::mutate(eligible = .data$patient %in% eligible)

  top_clones <- clone_counts |>
    dplyr::filter(.data$rank <= top_n) |>
    dplyr::group_by(.data$patient, .data$diagnosis) |>
    dplyr::summarise(
      cumulative_frequency = sum(.data$frequency),
      n_clones_counted = dplyr::n(),
      blood_cells = dplyr::first(.data$blood_cells),
      eligible = dplyr::first(.data$eligible),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$diagnosis, dplyr::desc(.data$cumulative_frequency))

  repertoire <- blood |>
    dplyr::group_by(.data$patient, .data$diagnosis) |>
    dplyr::summarise(
      blood_cells = dplyr::n(),
      n_clones = dplyr::n_distinct(.data$CTaa),
      n_expanded = sum(table(.data$CTaa) > 1L),
      fraction_in_expanded =
        sum(table(.data$CTaa)[table(.data$CTaa) > 1L]) / dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(eligible = .data$patient %in% eligible)

  list(
    clone_counts = clone_counts,
    top_clones = top_clones,
    repertoire = repertoire,
    excluded = dplyr::filter(blood_depth, .data$blood_cells < min_cells),
    min_cells = min_cells,
    top_n = top_n
  )
}

#' Per-diagnosis summary computed with and without the depth guard.
#'
#' Both are reported because the difference between them is the finding. Any
#' claim about blood clonal expansion in GBS has to state which set it rests on.
blood_expansion_summary <- function(expansion) {
  stopifnot(is.list(expansion), !is.null(expansion$top_clones))

  summarise_set <- function(data, label) {
    data |>
      dplyr::group_by(.data$diagnosis) |>
      dplyr::summarise(
        patient_set = label,
        n_patients = dplyr::n(),
        median_cumulative = stats::median(.data$cumulative_frequency),
        min_cumulative = min(.data$cumulative_frequency),
        max_cumulative = max(.data$cumulative_frequency),
        min_blood_cells = min(.data$blood_cells),
        .groups = "drop"
      )
  }

  dplyr::bind_rows(
    summarise_set(
      dplyr::filter(expansion$top_clones, .data$eligible), "evaluable only"
    ),
    summarise_set(expansion$top_clones, "all patients")
  ) |>
    dplyr::arrange(.data$patient_set, dplyr::desc(.data$median_cumulative))
}

write_blood_expansion_workbook <- function(expansion, summary, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(
    list(
      per_patient = expansion$top_clones,
      summary = summary,
      repertoire = expansion$repertoire,
      excluded = expansion$excluded
    ),
    path
  )
  path
}

#' Per-patient blood clonal expansion, with sampling depth on the face of it.
#'
#' Every patient is drawn. Those below the depth floor are open symbols, and
#' each point carries its blood cell count, so a reader can see immediately
#' which values are supported and which are a consequence of shallow sampling.
write_blood_expansion_plot <- function(expansion, colors, path,
                                       diagnoses = NULL, seed = 42L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  data <- expansion$top_clones
  if (!is.null(diagnoses)) {
    data <- dplyr::filter(data, .data$diagnosis %in% diagnoses)
  }
  stopifnot(nrow(data) > 0L)
  data <- order_by_diagnosis(data, colors)

  # The box summarises only evaluable patients; drawing it over the shallow
  # ones would reintroduce the artefact the guard exists to remove.
  evaluable <- dplyr::filter(data, .data$eligible)

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$diagnosis, y = .data$cumulative_frequency)
  ) +
    ggplot2::geom_boxplot(
      data = evaluable, ggplot2::aes(fill = .data$diagnosis),
      outlier.shape = NA, alpha = 0.35, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$diagnosis, shape = .data$eligible),
      position = ggplot2::position_jitter(width = 0.12, height = 0, seed = seed),
      size = 2.4
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$blood_cells),
      size = 2.3, colour = "grey30", seed = seed, max.overlaps = 20,
      segment.size = 0.2
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_colour_manual(values = colors, guide = "none") +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      name = paste0(expansion$min_cells, "+ blood cells")
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      x = NULL,
      y = paste0("Cumulative frequency of top ", expansion$top_n, " blood clones"),
      title = "Clonal expansion in blood",
      subtitle = paste(
        "One point per patient, labelled with its blood cell count.",
        "\nOpen symbols are too shallowly sampled to estimate repertoire structure"
      )
    )
  ggplot2::ggsave(path, plot = plot, width = 7, height = 5)
  path
}
