aligned_bh_result_dir <- function() {
  file.path("results", "targets", "proteomics")
}

aligned_bh_settings <- function(config) {
  fdr <- as.numeric(config$aligned_proteomics$fdr)
  stopifnot(length(fdr) == 1L, is.finite(fdr), fdr > 0, fdr < 1)
  list(fdr = fdr)
}

aligned_bh_contrast <- function(model) {
  suppressWarnings(
    emmeans::emmeans(model, "diagnosis") |>
      emmeans::contrast(list(difference = c(-1, 1)), adjust = "none") |>
      broom::tidy()
  )
}

aligned_bh_fit_olink_marker <- function(data, assay, group1, group2) {
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
  contrast <- aligned_bh_contrast(model)
  tibble::tibble(
    assay = assay, model = model_type,
    effect_scale = "adjusted log2 concentration difference",
    effect = contrast$estimate[[1L]], p_value = contrast$p.value[[1L]],
    observations = nrow(model_data),
    patients = dplyr::n_distinct(model_data$orbis_id), singular = singular
  )
}

aligned_bh_fit_luminex_marker <- function(
  result, assay, analysis, group1, group2, detection_sensitivity = FALSE
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
  if (detection_sensitivity) {
    model_data <- dplyr::filter(model_data, !is.na(.data$detected))
  }
  counts <- table(model_data$diagnosis)
  outcome <- if (detection_sensitivity) model_data$detected else model_data$value
  if (length(counts) < 2L || any(counts < 4L) ||
      length(unique(outcome)) < 2L) return(NULL)
  if (!detection_sensitivity) {
    model <- survival::survreg(
      survival::Surv(
        log2(analysis_value + 1), !censored, type = "left"
      ) ~ diagnosis + age + assay_batch,
      data = model_data, dist = "gaussian", robust = TRUE,
      cluster = patient_id
    )
    model_type <- "unified_censored_log2_gaussian"
    effect_scale <- "log2 concentration difference"
  } else {
    model <- stats::glm(
      detected ~ diagnosis + age + assay_batch,
      data = model_data, family = stats::binomial("logit"),
      method = brglm2::brglmFit, type = "AS_mean"
    )
    model_type <- "bias_reduced_logistic"
    effect_scale <- "log odds ratio"
  }
  contrast <- aligned_bh_contrast(model)
  tibble::tibble(
    assay = assay, analysis = analysis, model = model_type,
    effect_scale = effect_scale, effect = contrast$estimate[[1L]],
    p_value = contrast$p.value[[1L]], observations = nrow(model_data),
    patients = dplyr::n_distinct(model_data$patient_id), singular = NA
  )
}

aligned_bh_olink_lod_qc <- function(input, metadata, group1, group2) {
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

aligned_bh_luminex_marker_inventory <- function(result) {
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
      used_in_primary_analysis = .data$analysis %in% c(
        "primary_linear", "primary_censored", "primary_detection"
      ),
      included_in_primary_bh = .data$used_in_primary_analysis,
      primary_model = dplyr::if_else(
        .data$used_in_primary_analysis,
        "unified_censored_log2_gaussian", NA_character_
      ),
      reporting_status = dplyr::case_when(
        .data$used_in_primary_analysis ~ "primary_inference",
        .data$analysis == "exploratory_detection" ~ "exploratory_detection",
        .data$analysis == "descriptive_only" ~ "descriptive_only",
        .data$analysis == "no_test" ~ "not_testable",
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

analyze_aligned_bh_olink <- function(
  result, input, metadata, comparison, group1, group2, fdr = 0.1
) {
  statistics <- dplyr::bind_rows(lapply(result$assays, function(assay) {
    aligned_bh_fit_olink_marker(result$data_wide, assay, group1, group2)
  })) |>
    dplyr::mutate(
      comparison = comparison, group1 = group1, group2 = group2,
      bh_adjusted_p_value = stats::p.adjust(.data$p_value, method = "BH"),
      significant = .data$bh_adjusted_p_value < fdr
    ) |>
    dplyr::arrange(.data$p_value)
  list(
    statistics = statistics,
    lod_qc = aligned_bh_olink_lod_qc(input, metadata, group1, group2),
    platform = "Olink quantified", comparison = comparison,
    group1 = group1, group2 = group2, fdr = fdr
  )
}

analyze_aligned_bh_luminex <- function(
  result, comparison, group1, group2, fdr = 0.1
) {
  plan <- dplyr::filter(
    result$marker_plan,
    .data$analysis %in% c(
      "primary_linear", "primary_censored", "primary_detection"
    )
  )
  statistics <- dplyr::bind_rows(lapply(seq_len(nrow(plan)), function(index) {
    aligned_bh_fit_luminex_marker(
      result, plan$assay[[index]], plan$analysis[[index]], group1, group2
    )
  })) |>
    dplyr::mutate(
      comparison = comparison, group1 = group1, group2 = group2,
      bh_adjusted_p_value = stats::p.adjust(.data$p_value, method = "BH"),
      significant = .data$bh_adjusted_p_value < fdr
    ) |>
    dplyr::arrange(.data$p_value)
  sensitivity_plan <- dplyr::filter(plan, .data$analysis == "primary_detection")
  detection <- dplyr::bind_rows(lapply(
    seq_len(nrow(sensitivity_plan)), function(index) {
      aligned_bh_fit_luminex_marker(
        result, sensitivity_plan$assay[[index]],
        sensitivity_plan$analysis[[index]], group1, group2,
        detection_sensitivity = TRUE
      )
    }
  ))
  sensitivity_family <- dplyr::bind_rows(
    dplyr::filter(statistics, !.data$assay %in% sensitivity_plan$assay) |>
      dplyr::transmute(assay, p_value, source = "primary"),
    dplyr::transmute(
      detection, assay, p_value, source = "detection_sensitivity"
    )
  ) |>
    dplyr::mutate(
      sensitivity_bh_adjusted_p_value = stats::p.adjust(
        .data$p_value, method = "BH"
      )
    ) |>
    dplyr::filter(.data$source == "detection_sensitivity") |>
    dplyr::select("assay", "sensitivity_bh_adjusted_p_value")
  sensitivity <- dplyr::filter(
    statistics, .data$assay %in% sensitivity_plan$assay
  ) |>
    dplyr::select(
      "assay", "comparison", "group1", "group2",
      primary_effect = "effect", primary_p_value = "p_value",
      primary_bh_adjusted_p_value = "bh_adjusted_p_value"
    ) |>
    dplyr::left_join(
      dplyr::select(
        detection, "assay", detection_log_odds_ratio = "effect",
        detection_p_value = "p_value"
      ),
      by = "assay"
    ) |>
    dplyr::left_join(sensitivity_family, by = "assay")
  list(
    statistics = statistics, detection_sensitivity = sensitivity,
    marker_inventory = aligned_bh_luminex_marker_inventory(result),
    platform = "Luminex", comparison = comparison, group1 = group1,
    group2 = group2, fdr = fdr
  )
}

aligned_bh_volcano <- function(result) {
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

write_aligned_bh_artifacts <- function(result, platform) {
  root <- aligned_bh_result_dir()
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
  if (!is.null(result$detection_sensitivity)) {
    sheets$detection_sensitivity <- result$detection_sensitivity
  }
  if (!is.null(result$marker_inventory)) {
    sheets$marker_inventory <- result$marker_inventory
  }
  writexl::write_xlsx(sheets, paths[["statistics"]])
  ggplot2::ggsave(paths[["volcano"]], aligned_bh_volcano(result), width = 7, height = 5)
  unname(paths)
}

compare_aligned_bh_markers <- function(olink, luminex, markers) {
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

write_aligned_bh_comparison <- function(data) {
  path <- file.path(aligned_bh_result_dir(), "shared_marker_comparison.xlsx")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(data, path)
  path
}

write_aligned_bh_luminex_inventory <- function(result) {
  path <- file.path(aligned_bh_result_dir(), "luminex_marker_inventory.xlsx")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(
    aligned_bh_luminex_marker_inventory(result), path
  )
  path
}
