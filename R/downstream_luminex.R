luminex_numeric_value <- function(value) {
  suppressWarnings(as.numeric(value))
}

luminex_validate_numeric <- function(data, assays) {
  issues <- dplyr::bind_rows(lapply(assays, function(assay) {
    value <- data[[assay]]
    invalid <- !is.na(value) & is.na(suppressWarnings(as.numeric(value)))
    tibble::tibble(
      assay = assay, row = which(invalid), value = value[invalid]
    )
  }))
  if (nrow(issues)) {
    stop(
      "Non-numeric Luminex assay values found: ",
      paste0(issues$assay, "[", issues$row, "]=", issues$value, collapse = ", ")
    )
  }
  invisible(data)
}

read_luminex_input <- function(path) {
  raw <- readxl::read_xlsx(
    path, col_types = "text", na = c("", "NA", "-")
  )
  required <- c(
    "Origine", "IPP", "Age", "Diagnostic_category", "Diagnosis", "IL_1a"
  )
  stopifnot(all(required %in% names(raw)))
  assay_start <- match("IL_1a", names(raw))
  assay_end <- match("VEGALL", names(raw))
  stopifnot(!is.na(assay_start), !is.na(assay_end), assay_start < assay_end)
  assays <- names(raw)[seq.int(assay_start, assay_end)]
  luminex_validate_numeric(raw, assays)
  data <- raw |>
    dplyr::filter(
      .data$Diagnostic_category == "Headache" |
        .data$Diagnosis %in% c("CIDP", "GBS")
    ) |>
    dplyr::transmute(
      patient_id = as.character(.data$IPP),
      age = as.numeric(.data$Age),
      diagnosis = dplyr::if_else(
        .data$Diagnostic_category == "Headache", "CTRL", .data$Diagnosis
      ),
      assay_batch = dplyr::case_when(
        .data$Origine == "Second batch" ~ "Second batch",
        .data$Origine %in% c("First batch", "Test\u00e9 clinique") ~
          "First batch"
      ),
      dplyr::across(tidyselect::all_of(assays), luminex_numeric_value)
    )
  data$diagnosis <- factor(data$diagnosis, levels = c("CTRL", "GBS", "CIDP"))
  data$assay_batch <- factor(
    data$assay_batch, levels = c("First batch", "Second batch")
  )
  patient_diagnoses <- data |>
    dplyr::distinct(.data$patient_id, .data$diagnosis) |>
    dplyr::count(.data$patient_id)
  stopifnot(
    nrow(data) > 0L, !anyNA(data$patient_id), !anyNA(data$assay_batch),
    all(patient_diagnoses$n == 1L), all(table(data$diagnosis) >= 4L)
  )
  list(data = data, assays = assays)
}

luminex_contrasts <- function() {
  tibble::tribble(
    ~comparison, ~group1, ~group2,
    "CIDP_vs_CTRL", "CIDP", "CTRL",
    "GBS_vs_CTRL", "GBS", "CTRL",
    "CIDP_vs_GBS", "CIDP", "GBS"
  )
}
