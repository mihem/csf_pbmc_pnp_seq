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

write_table_workbook <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(table, path)
  path
}

write_named_table_workbooks <- function(tables, output_dir = table_result_dir()) {
  stopifnot(length(tables) > 0L, !is.null(names(tables)), !anyDuplicated(names(tables)))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(output_dir, paste0(names(tables), ".xlsx"))
  purrr::walk2(tables, paths, writexl::write_xlsx)
  paths
}
