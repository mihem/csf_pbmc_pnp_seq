luminex_llod_result_dir <- function() {
  file.path(luminex_result_dir(), "llod")
}

luminex_llod_marker_plan <- function() {
  plan <- tibble::tribble(
    ~assay, ~analysis,
    "IL_1a", "no_test",
    "IL_1RA", "primary_detection",
    "IL_1b", "descriptive_only",
    "Il_2", "descriptive_only",
    "Il_4", "no_test",
    "Il_5", "descriptive_only",
    "IL_6_cobas_", "primary_linear",
    "Il_7", "exploratory_detection",
    "Il_9", "no_test",
    "IL_10", "descriptive_only",
    "IL_12_p70", "no_test",
    "IL_13", "descriptive_only",
    "IL_15", "descriptive_only",
    "IL_16", "primary_censored",
    "IL_17A", "exploratory_detection",
    "IL_18", "no_test",
    "IL_21", "descriptive_only",
    "IL_22", "no_test",
    "IL_23", "no_test",
    "IL_27", "descriptive_only",
    "IL_29", "descriptive_only",
    "IL_31", "no_test",
    "IFN_a2", "no_test",
    "IFN_g", "no_test",
    "TNF_a", "no_test",
    "CCL1", "primary_detection",
    "CCL2", "primary_linear",
    "CCL3", "primary_detection",
    "CCL4", "primary_censored",
    "CCL5", "primary_censored",
    "CCL8", "primary_censored",
    "CCL11", "descriptive_only",
    "CCL17", "primary_detection",
    "CCL19", "primary_linear",
    "CCL20", "primary_censored",
    "CCL23", "primary_censored",
    "CXCL1", "primary_censored",
    "CXCL8", "primary_linear",
    "CXCL9", "primary_detection",
    "CXCL10", "primary_linear",
    "CXCL12", "primary_censored",
    "CXCL13", "primary_censored",
    "TNF_b", "no_test",
    "NGF_b", "no_test",
    "BDNF", "descriptive_only",
    "EGF", "descriptive_only",
    "FGF_2", "descriptive_only",
    "HGF", "primary_linear",
    "LIF", "descriptive_only",
    "PDGF_BB", "no_test",
    "PlGF_1", "descriptive_only",
    "SCF", "primary_censored",
    "BAFF", "exploratory_detection",
    "GM_CSF", "descriptive_only",
    "G_CSF", "descriptive_only",
    "MMP_1", "no_test",
    "IL2R", "excluded",
    "IL6_non_Cobas", "excluded",
    "VEFD", "excluded",
    "VEGALL", "excluded"
  )
  dplyr::mutate(
    plan,
    reason = dplyr::case_when(
      .data$analysis == "primary_linear" ~
        "No inferred LLOD; unified log-Gaussian model is fully observed",
      .data$assay == "CXCL12" ~
        "Second-batch LLOD only; batch-specific left-censored model",
      .data$analysis == "primary_censored" ~
        "Batch-specific LLOD; unified left-censored log-Gaussian model",
      .data$analysis == "primary_detection" ~
        "Highly censored; unified censored primary with detection sensitivity",
      .data$analysis == "exploratory_detection" ~
        "Highly censored with 5-6 detections; exploratory bias-reduced logistic model",
      .data$analysis == "descriptive_only" ~
        "Only 1-4 detections; descriptive reporting only",
      .data$analysis == "no_test" ~
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
      batch_data <- dplyr::filter(
        data, .data$assay_batch == .env$assay_batch
      )
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
        assay = assay,
        assay_batch = assay_batch,
        llod = llod,
        llod_basis = paste0(
          "Lowest value observed in >=", minimum_patients,
          " distinct patients within batch"
        ),
        observations = sum(observed),
        patients = dplyr::n_distinct(batch_data$patient_id[observed]),
        patients_at_llod = dplyr::n_distinct(batch_data$patient_id[at]),
        values_below_llod = sum(below),
        values_at_llod = sum(at),
        censored_values = sum(below | at),
        censored_percent = if (sum(observed)) {
          100 * sum(below | at) / sum(observed)
        } else {
          NA_real_
        }
      )
    }))
  }))
}

luminex_llod_assay_data <- function(data, assay, llod_qc) {
  data |>
    dplyr::transmute(
      patient_id, diagnosis, age, assay_batch,
      value = .data[[assay]]
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
        .data$value > .data$llod,
        NA
      ),
      analysis_value = dplyr::if_else(
        .data$censored, .data$llod, .data$value
      )
    )
}

luminex_llod_assay_qc <- function(data, assays, llod_qc) {
  dplyr::bind_rows(lapply(assays, function(assay) {
    assay_data <- luminex_llod_assay_data(data, assay, llod_qc)
    observed <- is.finite(assay_data$value)
    tibble::tibble(
      assay = assay,
      observations = sum(observed),
      batches_with_llod = sum(!is.na(
        dplyr::filter(llod_qc, .data$assay == .env$assay)$llod
      )),
      censored_values = sum(assay_data$censored[observed]),
      censored_percent = if (sum(observed)) {
        100 * sum(assay_data$censored[observed]) / sum(observed)
      } else {
        NA_real_
      },
      detected_count = sum(assay_data$detected, na.rm = TRUE)
    )
  }))
}

luminex_values_below_llod <- function(data, assays, llod_qc) {
  dplyr::bind_rows(lapply(assays, function(assay) {
    luminex_llod_assay_data(data, assay, llod_qc) |>
      dplyr::filter(
        is.finite(.data$value), !is.na(.data$llod), .data$value < .data$llod
      ) |>
      dplyr::transmute(
        assay = assay, patient_id, diagnosis, age, assay_batch, value, llod
      )
  }))
}

luminex_llod_emmeans <- function(model) {
  suppressWarnings(
    emmeans::emmeans(model, "diagnosis") |>
      emmeans::contrast(luminex_contrast_methods(), adjust = "none") |>
      broom::tidy()
  )
}

luminex_llod_fit_models <- function(
  data, llod_qc, marker_plan, analyses, adjust_bh = FALSE,
  force_censored = FALSE
) {
  selected <- dplyr::filter(marker_plan, .data$analysis %in% .env$analyses)
  results <- list()
  exclusions <- list()
  for (index in seq_len(nrow(selected))) {
    assay <- selected$assay[[index]]
    analysis <- selected$analysis[[index]]
    model_data <- luminex_llod_assay_data(data, assay, llod_qc) |>
      dplyr::filter(
        is.finite(.data$value), .data$value >= 0, !is.na(.data$age)
      )
    if (analysis %in% c("primary_detection", "exploratory_detection")) {
      model_data <- dplyr::filter(model_data, !is.na(.data$detected))
    }
    counts <- table(model_data$diagnosis)
    outcome <- if (analysis %in% c(
      "primary_detection", "exploratory_detection"
    )) model_data$detected else model_data$value
    if (length(counts) < 3L || any(counts < 4L) ||
        length(unique(outcome)) < 2L) {
      exclusions[[assay]] <- tibble::tibble(
        assay = assay, analysis = analysis,
        reason = "Insufficient variable observations by group"
      )
      next
    }
    fitted <- tryCatch({
      if (force_censored) {
        raw_model <- survival::survreg(
          survival::Surv(
            log2(analysis_value + 1), !censored, type = "left"
          ) ~ diagnosis,
          data = model_data, dist = "gaussian", robust = TRUE,
          cluster = patient_id
        )
        adjusted_model <- survival::survreg(
          survival::Surv(
            log2(analysis_value + 1), !censored, type = "left"
          ) ~ diagnosis + age + assay_batch,
          data = model_data, dist = "gaussian", robust = TRUE,
          cluster = patient_id
        )
        model_type <- "unified_censored_log2_gaussian"
        effect_scale <- "log2 concentration difference"
      } else if (analysis == "primary_linear") {
        raw_model <- stats::lm(
          log2(value + 1) ~ diagnosis, data = model_data
        )
        adjusted_model <- stats::lm(
          log2(value + 1) ~ diagnosis + age + assay_batch,
          data = model_data
        )
        model_type <- "linear_log2"
        effect_scale <- "log2 concentration difference"
      } else if (analysis == "primary_censored") {
        raw_model <- survival::survreg(
          survival::Surv(
            log2(analysis_value + 1), !censored, type = "left"
          ) ~ diagnosis,
          data = model_data, dist = "gaussian", robust = TRUE,
          cluster = patient_id
        )
        adjusted_model <- survival::survreg(
          survival::Surv(
            log2(analysis_value + 1), !censored, type = "left"
          ) ~ diagnosis + age + assay_batch,
          data = model_data, dist = "gaussian", robust = TRUE,
          cluster = patient_id
        )
        model_type <- "left_censored_log2"
        effect_scale <- "log2 concentration difference"
      } else {
        raw_model <- stats::glm(
          detected ~ diagnosis, data = model_data,
          family = stats::binomial("logit"), method = brglm2::brglmFit,
          type = "AS_mean"
        )
        adjusted_model <- stats::glm(
          detected ~ diagnosis + age + assay_batch, data = model_data,
          family = stats::binomial("logit"), method = brglm2::brglmFit,
          type = "AS_mean"
        )
        model_type <- "bias_reduced_logistic"
        effect_scale <- "log odds ratio"
      }
      raw <- luminex_llod_emmeans(raw_model) |>
        dplyr::select(comparison = "contrast", raw_p_value = "p.value")
      adjusted <- luminex_llod_emmeans(adjusted_model) |>
        dplyr::rename(
          comparison = "contrast", effect = "estimate",
          adjusted_p_value = "p.value"
        )
      dplyr::left_join(adjusted, raw, by = "comparison") |>
        dplyr::left_join(
          dplyr::select(luminex_contrasts(), "comparison", "group1", "group2"),
          by = "comparison"
        ) |>
        dplyr::mutate(
          assay = assay,
          analysis = analysis,
          model = model_type,
          effect_scale = effect_scale,
          odds_ratio = dplyr::if_else(
            .data$model == "bias_reduced_logistic", exp(.data$effect), NA_real_
          ),
          .before = 1L
        )
    }, error = function(error) {
      exclusions[[assay]] <<- tibble::tibble(
        assay = assay, analysis = analysis, reason = conditionMessage(error)
      )
      NULL
    })
    if (!is.null(fitted)) results[[assay]] <- fitted
  }
  statistics <- dplyr::bind_rows(results)
  if (adjust_bh && nrow(statistics)) {
    statistics <- statistics |>
      dplyr::group_by(.data$comparison) |>
      dplyr::mutate(
        bh_adjusted_p_value = stats::p.adjust(
          .data$adjusted_p_value, method = "BH"
        )
      ) |>
      dplyr::ungroup()
  }
  list(
    statistics = dplyr::arrange(statistics, .data$adjusted_p_value),
    exclusions = dplyr::bind_rows(exclusions)
  )
}

luminex_llod_detection_rates <- function(data, assays, llod_qc) {
  dplyr::bind_rows(lapply(assays, function(assay) {
    luminex_llod_assay_data(data, assay, llod_qc) |>
      dplyr::filter(!is.na(.data$detected)) |>
      dplyr::group_by(.data$diagnosis) |>
      dplyr::summarise(
        observations = dplyr::n(),
        detected_percent = 100 * mean(.data$detected),
        detected_count = sum(.data$detected),
        censored_count = sum(!.data$detected),
        .groups = "drop"
      ) |>
      dplyr::mutate(assay = assay, .before = 1L)
  }))
}

luminex_llod_detection_sensitivity <- function(primary, detection) {
  sensitivity_assays <- unique(detection$assay)
  sensitivity_family <- dplyr::bind_rows(
    dplyr::filter(primary, !.data$assay %in% .env$sensitivity_assays) |>
      dplyr::transmute(
        assay, comparison, adjusted_p_value, source = "primary"
      ),
    dplyr::transmute(
      detection, assay, comparison, adjusted_p_value,
      source = "detection_sensitivity"
    )
  ) |>
    dplyr::group_by(.data$comparison) |>
    dplyr::mutate(
      sensitivity_bh_adjusted_p_value = stats::p.adjust(
        .data$adjusted_p_value, method = "BH"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$source == "detection_sensitivity") |>
    dplyr::select(
      "assay", "comparison", "sensitivity_bh_adjusted_p_value"
    )
  dplyr::filter(primary, .data$assay %in% .env$sensitivity_assays) |>
    dplyr::select(
      "assay", "comparison", "group1", "group2",
      primary_effect = "effect", primary_p_value = "adjusted_p_value",
      primary_bh_adjusted_p_value = "bh_adjusted_p_value"
    ) |>
    dplyr::left_join(
      dplyr::select(
        detection, "assay", "comparison",
        detection_log_odds_ratio = "effect", "odds_ratio",
        detection_p_value = "adjusted_p_value"
      ),
      by = c("assay", "comparison")
    ) |>
    dplyr::left_join(
      sensitivity_family, by = c("assay", "comparison")
    ) |>
    dplyr::arrange(.data$comparison, .data$detection_p_value)
}

analyze_luminex_llod <- function(input) {
  marker_plan <- luminex_llod_marker_plan()
  stopifnot(setequal(input$assays, marker_plan$assay))
  llod_qc <- infer_luminex_llod(input$data, input$assays)
  assay_qc <- luminex_llod_assay_qc(input$data, input$assays, llod_qc) |>
    dplyr::left_join(marker_plan, by = "assay")
  primary <- luminex_llod_fit_models(
    input$data, llod_qc, marker_plan,
    c("primary_linear", "primary_censored", "primary_detection"),
    adjust_bh = TRUE, force_censored = TRUE
  )
  detection_sensitivity <- luminex_llod_fit_models(
    input$data, llod_qc, marker_plan, "primary_detection"
  )
  exploratory <- luminex_llod_fit_models(
    input$data, llod_qc, marker_plan, "exploratory_detection"
  )
  sensitivity <- luminex_llod_detection_sensitivity(
    primary$statistics, detection_sensitivity$statistics
  )
  list(
    data_wide = input$data,
    assays = input$assays,
    marker_plan = marker_plan,
    llod_by_batch = llod_qc,
    assay_qc = assay_qc,
    below_llod_values = luminex_values_below_llod(
      input$data, input$assays, llod_qc
    ),
    primary_statistics = primary$statistics,
    detection_sensitivity = sensitivity,
    exploratory_statistics = exploratory$statistics,
    model_exclusions = dplyr::bind_rows(
      primary$exclusions, detection_sensitivity$exclusions,
      exploratory$exclusions
    ),
    detection_rates = luminex_llod_detection_rates(
      input$data, input$assays, llod_qc
    )
  )
}

luminex_llod_boxplot <- function(result) {
  colors <- stats::setNames(
    pals::cols25(nlevels(result$data_wide$diagnosis)),
    levels(result$data_wide$diagnosis)
  )
  plots <- lapply(gtools::mixedsort(result$assays), function(assay) {
    qc <- dplyr::filter(result$assay_qc, .data$assay == .env$assay)
    plot_data <- luminex_llod_assay_data(
      result$data_wide, assay, result$llod_by_batch
    ) |>
      dplyr::mutate(
        measurement = dplyr::case_when(
          !is.finite(.data$value) ~ NA_character_,
          is.na(.data$llod) ~ "No inferred LLOD",
          .data$censored ~ "<= batch LLOD",
          TRUE ~ "Detected"
        ),
        plot_value = log2(.data$value + 1)
      )
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = diagnosis, y = plot_value, fill = diagnosis)
    ) +
      ggplot2::geom_boxplot(na.rm = TRUE, outlier.shape = NA) +
      ggplot2::geom_jitter(
        ggplot2::aes(shape = measurement), width = 0.2, height = 0,
        na.rm = TRUE
      ) +
      ggplot2::facet_wrap(~assay_batch) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_shape_manual(values = c(
        "<= batch LLOD" = 25, "Detected" = 21, "No inferred LLOD" = 1
      )) +
      ggplot2::labs(
        title = paste0(
          assay, " (", round(qc$censored_percent[[1L]], 1), "% censored; ",
          qc$analysis[[1L]], ")"
        ),
        x = NULL, y = "log2(pg/ml + 1)", shape = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")
  })
  patchwork::wrap_plots(plots)
}

write_luminex_llod_artifacts <- function(result) {
  root <- luminex_llod_result_dir()
  paths <- c(
    object = file.path(root, "luminex_llod.qs"),
    statistics = file.path(root, "luminex_llod_statistics.xlsx"),
    boxplots = file.path(root, "luminex_llod_boxplots.pdf")
  )
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  qs::qsave(result, paths[["object"]])
  writexl::write_xlsx(
    list(
      primary_statistics = result$primary_statistics,
      detection_sensitivity = result$detection_sensitivity,
      exploratory_detection = result$exploratory_statistics,
      detection_rates = result$detection_rates,
      llod_by_batch = result$llod_by_batch,
      assay_qc = result$assay_qc,
      below_llod_values = result$below_llod_values,
      model_exclusions = result$model_exclusions,
      marker_plan = result$marker_plan
    ),
    paths[["statistics"]]
  )
  ggplot2::ggsave(
    paths[["boxplots"]], luminex_llod_boxplot(result),
    width = 14, height = 36, limitsize = FALSE
  )
  unname(paths)
}

write_luminex_llod_volcano <- function(
  result, comparison, group1, group2, seed
) {
  data <- result$primary_statistics |>
    dplyr::filter(.data$comparison == .env$comparison) |>
    dplyr::mutate(
      neg_log10_q = -log10(.data$bh_adjusted_p_value),
      significant = .data$bh_adjusted_p_value < 0.1
    )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = effect, y = neg_log10_q, label = assay, color = significant
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(
      yintercept = -log10(0.1), color = "blue", linetype = "dashed"
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "red") +
    ggrepel::geom_text_repel(
      data = dplyr::filter(data, .data$significant), seed = seed
    ) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
    ggplot2::labs(
      x = paste0(
        "Adjusted underlying log2 concentration difference: ",
        group1, " - ", group2
      ),
      y = expression(-Log[10] ~ "BH-adjusted p value"),
      color = "BH < 0.1"
    ) +
    ggplot2::theme_classic()
  path <- file.path(
    luminex_llod_result_dir(), paste0("volcano_", comparison, "_llod.pdf")
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = 6, height = 5)
  path
}
