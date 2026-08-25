demographics_result_dir <- function() {
  file.path("results", "targets", "demographics")
}

demographics_boxplot_config <- function() {
  tibble::tribble(
    ~name, ~variable, ~title, ~geom, ~width, ~height,
    "age", "age", "age", "point", 3.5, 3,
    "ncv_tibial_motoric", "ncv_tibial_motoric",
      "motoric NCV tibial nerve (m/s)", "jitter", 3, 3,
    "incat", "incat_at_lumbar_puncture", "INCAT score", "jitter", 3, 3,
    "incat_progress", "incat_progress", "INCAT score progress", "jitter", 5, 5,
    "disease_duration", "disease_duration_in_months",
      "disease duration (months)", "jitter", 3.5, 3,
    "csf_protein", "csf_protein", "CSF protein (mg/L)", "jitter", 3.8, 3
  )
}

demographics_barplot_config <- function() {
  tibble::tribble(
    ~name, ~variable, ~title, ~palette, ~width, ~height,
    "disease", "diagnosis", "diagnosis", "diagnosis", 4.5, 3,
    "sex", "sex", "sex", "sex", 4.5, 3,
    "therapy_status", "therapy", "therapy status", "binary", 4.5, 3,
    "bbb_dysfunction", "bbbd", "BBB dysfunction", "bbb", 5, 5
  )
}

prepare_olink_demographics <- function(quant_file, metadata_file) {
  quant <- readxl::read_xlsx(quant_file)
  metadata <- readxl::read_xlsx(metadata_file)
  required <- c("SampleID", "orbis_id", "age", "sex", "diagnosis")
  stopifnot("SampleID" %in% names(quant), all(required %in% names(metadata)))
  result <- quant |>
    dplyr::select("SampleID") |>
    dplyr::distinct() |>
    dplyr::left_join(metadata, by = "SampleID") |>
    dplyr::select(tidyselect::all_of(required)) |>
    dplyr::distinct(.data$SampleID, .data$orbis_id, .keep_all = TRUE)
  result$diagnosis <- factor(result$diagnosis, levels = c("CTRL", "GBS", "CIDP"))
  stopifnot(!anyNA(result$orbis_id), !anyDuplicated(result$SampleID))
  result
}

prepare_flow_demographics <- function(flow) {
  stopifnot(is.list(flow), "blood" %in% names(flow))
  result <- flow$blood |>
    dplyr::filter(.data$diagnosis %in% c("CTRL", "GBS", "CIDP"))
  result$diagnosis <- factor(result$diagnosis, levels = c("CTRL", "GBS", "CIDP"))
  result
}

demographics_overview <- function(data, patient_column = "orbis_id") {
  stopifnot(all(c(patient_column, "diagnosis", "sex", "age") %in% names(data)))
  data |>
    dplyr::mutate(sex_cat = dplyr::if_else(.data$sex == "male", 1, 0)) |>
    dplyr::group_by(.data$diagnosis) |>
    dplyr::summarise(
      samples = dplyr::n(),
      patients = dplyr::n_distinct(.data[[patient_column]]),
      age = mean(.data$age, na.rm = TRUE),
      female = (1 - mean(.data$sex_cat, na.rm = TRUE)) * 100,
      .groups = "drop"
    )
}

demographics_boxplot <- function(data, variable, title, geom, colors) {
  plot_data <- dplyr::filter(data, !is.na(.data[[variable]]))
  statistics <- scMisc:::compStat(
    x_var = variable, group = "diagnosis", data = plot_data, paired = FALSE
  )
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = diagnosis, y = .data[[variable]], fill = diagnosis)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  if (all(c("comparisons", "annotation") %in% names(statistics)) &&
      length(statistics$comparisons)) {
    plot <- plot + ggsignif::geom_signif(
      comparisons = statistics$comparisons,
      annotations = statistics$annotation,
      textsize = 5,
      step_increase = 0.05,
      vjust = 0.7
    )
  }
  if (geom == "point") plot + ggplot2::geom_point() else
    plot + ggplot2::geom_jitter(height = 0, width = 0.3)
}

demographics_barplot <- function(data, variable, title, colors) {
  ggplot2::ggplot(
    data, ggplot2::aes(x = diagnosis, fill = .data[[variable]])
  ) +
    ggplot2::geom_bar() +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_bw()
}

write_demographics_plot <- function(plot, name, width, height, seed = 123L) {
  path <- file.path(demographics_result_dir(), name)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  withr::with_seed(seed, ggplot2::ggsave(path, plot, width = width, height = height))
  path
}

write_scrna_demographics <- function(lookup, diagnosis_colors, seed = 123L) {
  data <- lookup
  boxes <- demographics_boxplot_config()
  bars <- demographics_barplot_config()
  box_files <- vapply(seq_len(nrow(boxes)), function(index) {
    plot <- demographics_boxplot(
      data, boxes$variable[[index]], boxes$title[[index]], boxes$geom[[index]],
      diagnosis_colors
    )
    write_demographics_plot(
      plot, paste0("boxplot_", boxes$name[[index]], ".pdf"),
      boxes$width[[index]], boxes$height[[index]], seed + index - 1L
    )
  }, character(1))
  palettes <- list(
    diagnosis = pals::cols25(9),
    sex = rev(pals::cols25(2)),
    binary = pals::cols25(2),
    bbb = c(no = "#1F78C8", yes = "#FF0000")
  )
  bar_files <- vapply(seq_len(nrow(bars)), function(index) {
    plot <- demographics_barplot(
      data, bars$variable[[index]], bars$title[[index]],
      palettes[[bars$palette[[index]]]]
    )
    write_demographics_plot(
      plot, paste0("barplot_", bars$name[[index]], ".pdf"),
      bars$width[[index]], bars$height[[index]], seed + 100L + index
    )
  }, character(1))
  c(box_files, bar_files)
}

write_modality_demographics <- function(
  data, modality, patient_column, diagnosis_colors, seed = 123L
) {
  overview <- demographics_overview(data, patient_column)
  table_path <- file.path(
    demographics_result_dir(), paste0("overview_table_", modality, ".xlsx")
  )
  dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(overview, table_path)
  age <- demographics_boxplot(data, "age", "age", "point", diagnosis_colors)
  disease <- demographics_barplot(data, "diagnosis", "diagnosis", pals::cols25(9))
  sex <- demographics_barplot(data, "sex", "sex", rev(pals::cols25(2)))
  c(
    table_path,
    write_demographics_plot(
      age, paste0("boxplot_age_", modality, ".pdf"), 2, 3, seed
    ),
    write_demographics_plot(
      disease, paste0("barplot_disease_", modality, ".pdf"), 3, 3, seed + 1L
    ),
    write_demographics_plot(
      sex, paste0("barplot_sex_", modality, ".pdf"), 3, 3, seed + 2L
    )
  )
}

check_pdf_artifacts <- function(paths) {
  stopifnot(length(paths) > 0L, all(file.exists(paths)), all(file.size(paths) > 0L))
  headers <- vapply(paths, function(path) {
    rawToChar(readBin(path, what = "raw", n = 4L))
  }, character(1))
  stopifnot(all(headers == "%PDF"))
  tibble::tibble(artifact = basename(paths), passed = TRUE)
}
