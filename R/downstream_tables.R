table_result_dir <- function() {
  file.path("results", "targets", "table")
}

write_scrna_overview_table <- function(lookup, path) {
  table <- lookup |>
    dplyr::mutate(sex_cat = dplyr::if_else(.data$sex == "male", 1, 0)) |>
    dplyr::group_by(.data$diagnosis) |>
    dplyr::summarise(
      n = dplyr::n(),
      age = mean(.data$age, na.rm = TRUE),
      female = (1 - mean(.data$sex_cat, na.rm = TRUE)) * 100,
      .groups = "drop"
    )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(table, path)
  path
}

collect_publication_tables <- function(paths, output_dir = table_result_dir()) {
  paths <- unique(paths[grepl("[.]xlsx$", paths)])
  stopifnot(length(paths) > 0L, all(file.exists(paths)), !anyDuplicated(basename(paths)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- file.path(output_dir, basename(paths))
  copy <- normalizePath(paths, mustWork = TRUE) != normalizePath(
    outputs, mustWork = FALSE
  )
  copied <- file.copy(
    paths[copy], outputs[copy], overwrite = TRUE, copy.mode = FALSE
  )
  stopifnot(all(copied), all(file.exists(outputs)), all(file.size(outputs) > 0L))
  outputs
}
