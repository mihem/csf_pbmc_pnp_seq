proteomics_result_dir <- function() {
  file.path("results", "targets", "proteomics")
}

proteomics_settings <- function(config) {
  fdr <- as.numeric(config$proteomics$fdr)
  stopifnot(length(fdr) == 1L, is.finite(fdr), fdr > 0, fdr < 1)
  list(fdr = fdr)
}

proteomics_contrast <- function(model) {
  suppressWarnings(
    emmeans::emmeans(model, "diagnosis") |>
      emmeans::contrast(list(difference = c(-1, 1)), adjust = "none") |>
      broom::tidy()
  )
}

proteomics_fit_olink_marker <- function(data, assay, group1, group2) {
  model_data <- data |>
    dplyr::filter(
      .data$diagnosis %in% c(.env$group1, .env$group2),
      is.finite(.data[[assay]]), !is.na(.data$age), !is.na(.data$sex)
    ) |>
    dplyr::mutate(
      diagnosis = factor(.data$diagnosis, levels = c(group2, group1)),
      outcome = log2(.data[[assay]] + 1)
    )
  counts <- table(model_data$diagnosis)
  if (length(counts) < 2L || any(counts < 4L)) return(NULL)
  if (dplyr::n_distinct(model_data$orbis_id) < nrow(model_data)) {
    model <- suppressWarnings(suppressMessages(lme4::lmer(
      outcome ~ diagnosis + age + sex + (1 | orbis_id), data = model_data
    )))
    model_type <- "linear_mixed_log2_quantified"
    singular <- lme4::isSingular(model)
  } else {
    model <- stats::lm(outcome ~ diagnosis + age + sex, data = model_data)
    model_type <- "linear_log2_quantified"
    singular <- NA
  }
  contrast <- proteomics_contrast(model)
  tibble::tibble(
    assay = assay, model = model_type,
    effect_scale = "adjusted log2 concentration difference",
    effect = contrast$estimate[[1L]], p_value = contrast$p.value[[1L]],
    observations = nrow(model_data),
    patients = dplyr::n_distinct(model_data$orbis_id), singular = singular
  )
}

proteomics_fit_luminex_marker <- function(
  result, assay, marker_status, group1, group2
) {
  model_data <- luminex_llod_assay_data(
    result$data_wide, assay, result$llod_by_batch
  ) |>
    dplyr::filter(
      .data$diagnosis %in% c(.env$group1, .env$group2),
      is.finite(.data$value), .data$value >= 0, !is.na(.data$age)
    ) |>
    dplyr::mutate(
      diagnosis = factor(.data$diagnosis, levels = c(group2, group1))
    )
  counts <- table(model_data$diagnosis)
  if (length(counts) < 2L || any(counts < 4L) ||
      length(unique(model_data$value)) < 2L) return(NULL)
  model <- survival::survreg(
    survival::Surv(
      log2(analysis_value + 1), !censored, type = "left"
    ) ~ diagnosis + age + assay_batch,
    data = model_data, dist = "gaussian", robust = TRUE,
    cluster = patient_id
  )
  contrast <- proteomics_contrast(model)
  tibble::tibble(
    assay = assay, marker_status = marker_status,
    model = "unified_censored_log2_gaussian",
    effect_scale = "log2 concentration difference",
    effect = contrast$estimate[[1L]],
    p_value = contrast$p.value[[1L]], observations = nrow(model_data),
    patients = dplyr::n_distinct(model_data$patient_id), singular = NA
  )
}

proteomics_olink_lod_qc <- function(input, metadata, group1, group2) {
  input |>
    dplyr::left_join(
      dplyr::select(metadata, "SampleID", "diagnosis"), by = "SampleID"
    ) |>
    dplyr::filter(.data$diagnosis %in% c(.env$group1, .env$group2)) |>
    dplyr::mutate(
      PlateLOD = as.numeric(.data$PlateLOD), LLOQ = as.numeric(.data$LLOQ),
      below_lod = is.finite(.data$Quantified_value) &
        is.finite(.data$PlateLOD) & .data$Quantified_value < .data$PlateLOD,
      below_lloq = is.finite(.data$Quantified_value) &
        is.finite(.data$LLOQ) & .data$Quantified_value < .data$LLOQ
    ) |>
    dplyr::group_by(.data$Assay, .data$diagnosis) |>
    dplyr::summarise(
      samples = dplyr::n(), quantified = sum(is.finite(.data$Quantified_value)),
      below_lod = sum(.data$below_lod), below_lloq = sum(.data$below_lloq),
      below_lod_percent = 100 * .data$below_lod / .data$samples,
      below_lloq_percent = 100 * .data$below_lloq / .data$samples,
      .groups = "drop"
    ) |>
    dplyr::rename(assay = "Assay")
}

proteomics_luminex_marker_inventory <- function(result) {
  rates <- result$detection_rates |>
    dplyr::select(
      "assay", "diagnosis", "detected_count", "detected_percent"
    ) |>
    tidyr::pivot_wider(
      names_from = "diagnosis",
      values_from = c("detected_count", "detected_percent"),
      names_glue = "{diagnosis}_{.value}"
    )
  llods <- result$llod_by_batch |>
    dplyr::select("assay", "assay_batch", "llod") |>
    tidyr::pivot_wider(
      names_from = "assay_batch", values_from = "llod",
      names_glue = "llod_{assay_batch}"
    ) |>
    janitor::clean_names()
  result$assay_qc |>
    dplyr::mutate(
      used_in_primary_analysis = .data$marker_status %in% c(
        "primary_fully_observed", "primary_left_censored",
        "primary_high_censoring"
      ),
      included_in_primary_bh = .data$used_in_primary_analysis,
      primary_model = dplyr::if_else(
        .data$used_in_primary_analysis,
        "unified_censored_log2_gaussian", NA_character_
      ),
      reporting_status = dplyr::case_when(
        .data$used_in_primary_analysis ~ "primary_inference",
        .data$marker_status == "exploratory_low_detection" ~
          "exploratory_low_detection",
        .data$marker_status == "descriptive_only" ~ "descriptive_only",
        .data$marker_status == "no_test" ~ "not_testable",
        TRUE ~ "excluded"
      ),
      why_not_primary = dplyr::if_else(
        .data$used_in_primary_analysis, NA_character_, .data$reason
      )
    ) |>
    dplyr::left_join(rates, by = "assay") |>
    dplyr::left_join(llods, by = "assay") |>
    dplyr::arrange(
      dplyr::desc(.data$used_in_primary_analysis), .data$reporting_status,
      .data$assay
    )
}

analyze_proteomics_olink <- function(
  result, input, metadata, comparison, group1, group2, fdr = 0.1
) {
  statistics <- dplyr::bind_rows(lapply(result$assays, function(assay) {
    proteomics_fit_olink_marker(result$data_wide, assay, group1, group2)
  })) |>
    dplyr::mutate(
      comparison = comparison, group1 = group1, group2 = group2,
      bh_adjusted_p_value = stats::p.adjust(.data$p_value, method = "BH"),
      significant = .data$bh_adjusted_p_value < fdr
    ) |>
    dplyr::arrange(.data$p_value)
  list(
    statistics = statistics,
    lod_qc = proteomics_olink_lod_qc(input, metadata, group1, group2),
    platform = "Olink quantified", comparison = comparison,
    group1 = group1, group2 = group2, fdr = fdr
  )
}

analyze_proteomics_luminex <- function(
  result, comparison, group1, group2, fdr = 0.1
) {
  plan <- dplyr::filter(
    result$marker_plan,
    .data$marker_status %in% c(
      "primary_fully_observed", "primary_left_censored",
      "primary_high_censoring"
    )
  )
  statistics <- dplyr::bind_rows(lapply(seq_len(nrow(plan)), function(index) {
    proteomics_fit_luminex_marker(
      result, plan$assay[[index]], plan$marker_status[[index]], group1, group2
    )
  })) |>
    dplyr::mutate(
      comparison = comparison, group1 = group1, group2 = group2,
      bh_adjusted_p_value = stats::p.adjust(.data$p_value, method = "BH"),
      significant = .data$bh_adjusted_p_value < fdr
    ) |>
    dplyr::arrange(.data$p_value)
  list(
    statistics = statistics,
    marker_inventory = proteomics_luminex_marker_inventory(result),
    platform = "Luminex", comparison = comparison, group1 = group1,
    group2 = group2, fdr = fdr
  )
}

proteomics_volcano <- function(result) {
  data <- dplyr::mutate(
    result$statistics, neg_log10_q = -log10(.data$bh_adjusted_p_value)
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = effect, y = neg_log10_q, label = assay, color = significant
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(
      yintercept = -log10(result$fdr), color = "blue", linetype = "dashed"
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "red") +
    ggrepel::geom_text_repel(
      data = dplyr::filter(data, .data$significant), seed = 900L
    ) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
    ggplot2::labs(
      title = paste(result$platform, result$comparison),
      subtitle = "Separate two-group adjusted model; BH within comparison",
      x = paste0(
        "Adjusted log2 concentration difference: ", result$group1,
        " - ", result$group2
      ),
      y = expression(-Log[10] ~ "BH-adjusted p value"),
      color = "BH significant"
    ) +
    ggplot2::theme_classic()
}

write_proteomics_artifacts <- function(result, platform) {
  root <- proteomics_result_dir()
  prefix <- paste0(platform, "_", result$comparison, "_bh")
  paths <- c(
    object = file.path(root, paste0(prefix, ".qs")),
    statistics = file.path(root, paste0(prefix, ".xlsx")),
    volcano = file.path(root, paste0(prefix, ".pdf"))
  )
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  qs::qsave(result, paths[["object"]])
  sheets <- list(statistics = result$statistics)
  if (!is.null(result$lod_qc)) sheets$lod_qc <- result$lod_qc
  if (!is.null(result$marker_inventory)) {
    sheets$marker_inventory <- result$marker_inventory
  }
  writexl::write_xlsx(sheets, paths[["statistics"]])
  ggplot2::ggsave(paths[["volcano"]], proteomics_volcano(result), width = 7, height = 5)
  unname(paths)
}

compare_proteomics_markers <- function(olink, luminex, markers) {
  olink$statistics |>
    dplyr::filter(.data$assay %in% .env$markers) |>
    dplyr::select(
      "assay", "comparison", olink_effect = "effect",
      olink_p_value = "p_value", olink_q_value = "bh_adjusted_p_value",
      olink_significant = "significant"
    ) |>
    dplyr::left_join(
      luminex$statistics |>
        dplyr::filter(.data$assay %in% .env$markers) |>
        dplyr::select(
          "assay", "comparison", luminex_model = "model",
          luminex_effect = "effect", luminex_p_value = "p_value",
          luminex_q_value = "bh_adjusted_p_value",
          luminex_significant = "significant"
        ),
      by = c("assay", "comparison")
    ) |>
    dplyr::mutate(
      direction_agrees = sign(.data$olink_effect) == sign(.data$luminex_effect),
      significant_in_both = .data$olink_significant & .data$luminex_significant
    )
}

write_proteomics_comparison <- function(data) {
  path <- file.path(proteomics_result_dir(), "shared_marker_comparison.xlsx")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(data, path)
  path
}

write_proteomics_luminex_inventory <- function(result) {
  path <- file.path(proteomics_result_dir(), "luminex_marker_inventory.xlsx")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(
    proteomics_luminex_marker_inventory(result), path
  )
  path
}
