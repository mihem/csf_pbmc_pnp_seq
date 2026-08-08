luminex_llod_marker_plan <- function() {
  plan <- tibble::tribble(
    ~assay, ~marker_status,
    "IL_1a", "no_test", "IL_1RA", "primary_high_censoring",
    "IL_1b", "descriptive_only", "Il_2", "descriptive_only",
    "Il_4", "no_test", "Il_5", "descriptive_only",
    "IL_6_cobas_", "primary_fully_observed", "Il_7", "exploratory_low_detection",
    "Il_9", "no_test", "IL_10", "descriptive_only",
    "IL_12_p70", "no_test", "IL_13", "descriptive_only",
    "IL_15", "descriptive_only", "IL_16", "primary_left_censored",
    "IL_17A", "exploratory_low_detection", "IL_18", "no_test",
    "IL_21", "descriptive_only", "IL_22", "no_test",
    "IL_23", "no_test", "IL_27", "descriptive_only",
    "IL_29", "descriptive_only", "IL_31", "no_test",
    "IFN_a2", "no_test", "IFN_g", "no_test", "TNF_a", "no_test",
    "CCL1", "primary_high_censoring", "CCL2", "primary_fully_observed",
    "CCL3", "primary_high_censoring", "CCL4", "primary_left_censored",
    "CCL5", "primary_left_censored", "CCL8", "primary_left_censored",
    "CCL11", "descriptive_only", "CCL17", "primary_high_censoring",
    "CCL19", "primary_fully_observed", "CCL20", "primary_left_censored",
    "CCL23", "primary_left_censored", "CXCL1", "primary_left_censored",
    "CXCL8", "primary_fully_observed", "CXCL9", "primary_high_censoring",
    "CXCL10", "primary_fully_observed", "CXCL12", "primary_left_censored",
    "CXCL13", "primary_left_censored", "TNF_b", "no_test",
    "NGF_b", "no_test", "BDNF", "descriptive_only",
    "EGF", "descriptive_only", "FGF_2", "descriptive_only",
    "HGF", "primary_fully_observed", "LIF", "descriptive_only",
    "PDGF_BB", "no_test", "PlGF_1", "descriptive_only",
    "SCF", "primary_left_censored", "BAFF", "exploratory_low_detection",
    "GM_CSF", "descriptive_only", "G_CSF", "descriptive_only",
    "MMP_1", "no_test", "IL2R", "excluded",
    "IL6_non_Cobas", "excluded", "VEFD", "excluded", "VEGALL", "excluded"
  )
  dplyr::mutate(
    plan,
    reason = dplyr::case_when(
      .data$marker_status == "primary_fully_observed" ~
        "No inferred LLOD; unified log-Gaussian model is fully observed",
      .data$assay == "CXCL12" ~
        "Second-batch LLOD only; unified left-censored log-Gaussian model",
      .data$marker_status == "primary_left_censored" ~
        "Batch-specific LLOD; unified left-censored log-Gaussian model",
      .data$marker_status == "primary_high_censoring" ~
        "Highly censored; unified left-censored log-Gaussian model",
      .data$marker_status == "exploratory_low_detection" ~
        "Highly censored with 5-6 detections; exploratory only",
      .data$marker_status == "descriptive_only" ~
        "Only 1-4 detections; descriptive reporting only",
      .data$marker_status == "no_test" ~
        "No detections above the batch-specific LLOD",
      .data$assay %in% c("IL2R", "IL6_non_Cobas") ~
        "Strong batch-dependent availability",
      TRUE ~ "Insufficient observations"
    )
  )
}

infer_luminex_llod <- function(data, assays, minimum_patients = 5L) {
  dplyr::bind_rows(lapply(assays, function(assay) {
    dplyr::bind_rows(lapply(levels(data$assay_batch), function(assay_batch) {
      batch_data <- dplyr::filter(data, .data$assay_batch == .env$assay_batch)
      patient_values <- tibble::tibble(
        patient_id = batch_data$patient_id, value = batch_data[[assay]]
      ) |>
        dplyr::filter(is.finite(.data$value)) |>
        dplyr::distinct(.data$patient_id, .data$value)
      counts <- table(patient_values$value)
      candidates <- as.numeric(names(counts)[counts >= minimum_patients])
      llod <- if (length(candidates)) min(candidates) else NA_real_
      values <- batch_data[[assay]]
      observed <- is.finite(values)
      below <- observed & !is.na(llod) & values < llod
      at <- observed & !is.na(llod) & values == llod
      tibble::tibble(
        assay = assay, assay_batch = assay_batch, llod = llod,
        llod_basis = paste0(
          "Lowest value observed in >=", minimum_patients,
          " distinct patients within batch"
        ),
        observations = sum(observed),
        patients = dplyr::n_distinct(batch_data$patient_id[observed]),
        patients_at_llod = dplyr::n_distinct(batch_data$patient_id[at]),
        values_below_llod = sum(below), values_at_llod = sum(at),
        censored_values = sum(below | at),
        censored_percent = if (sum(observed)) {
          100 * sum(below | at) / sum(observed)
        } else NA_real_
      )
    }))
  }))
}

luminex_llod_assay_data <- function(data, assay, llod_qc) {
  data |>
    dplyr::transmute(
      patient_id, diagnosis, age, assay_batch, value = .data[[assay]]
    ) |>
    dplyr::left_join(
      dplyr::filter(llod_qc, .data$assay == .env$assay) |>
        dplyr::select("assay_batch", "llod"),
      by = "assay_batch"
    ) |>
    dplyr::mutate(
      censored = is.finite(.data$value) & !is.na(.data$llod) &
        .data$value <= .data$llod,
      detected = dplyr::if_else(
        is.finite(.data$value) & !is.na(.data$llod),
        .data$value > .data$llod, NA
      ),
      analysis_value = dplyr::if_else(.data$censored, .data$llod, .data$value)
    )
}

luminex_llod_assay_qc <- function(data, assays, llod_qc) {
  dplyr::bind_rows(lapply(assays, function(assay) {
    assay_data <- luminex_llod_assay_data(data, assay, llod_qc)
    observed <- is.finite(assay_data$value)
    tibble::tibble(
      assay = assay, observations = sum(observed),
      batches_with_llod = sum(!is.na(
        dplyr::filter(llod_qc, .data$assay == .env$assay)$llod
      )),
      censored_values = sum(assay_data$censored[observed]),
      censored_percent = if (sum(observed)) {
        100 * sum(assay_data$censored[observed]) / sum(observed)
      } else NA_real_,
      detected_count = sum(assay_data$detected, na.rm = TRUE)
    )
  }))
}

luminex_llod_detection_rates <- function(data, assays, llod_qc) {
  dplyr::bind_rows(lapply(assays, function(assay) {
    luminex_llod_assay_data(data, assay, llod_qc) |>
      dplyr::filter(!is.na(.data$detected)) |>
      dplyr::group_by(.data$diagnosis) |>
      dplyr::summarise(
        observations = dplyr::n(), detected_percent = 100 * mean(.data$detected),
        detected_count = sum(.data$detected),
        censored_count = sum(!.data$detected), .groups = "drop"
      ) |>
      dplyr::mutate(assay = assay, .before = 1L)
  }))
}

analyze_luminex_llod <- function(input) {
  marker_plan <- luminex_llod_marker_plan()
  stopifnot(setequal(input$assays, marker_plan$assay))
  llod_qc <- infer_luminex_llod(input$data, input$assays)
  list(
    data_wide = input$data, assays = input$assays,
    marker_plan = marker_plan, llod_by_batch = llod_qc,
    assay_qc = luminex_llod_assay_qc(
      input$data, input$assays, llod_qc
    ) |>
      dplyr::left_join(marker_plan, by = "assay"),
    detection_rates = luminex_llod_detection_rates(
      input$data, input$assays, llod_qc
    )
  )
}
