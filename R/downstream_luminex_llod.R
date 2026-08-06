luminex_llod_result_dir <- function() {
  file.path(luminex_result_dir(), "llod")
}

infer_luminex_llod <- function(data, assays, minimum_patients = 5L) {
  dplyr::bind_rows(lapply(assays, function(assay) {
    patient_values <- tibble::tibble(
      patient_id = data$patient_id, value = data[[assay]]
    ) |>
      dplyr::filter(is.finite(.data$value)) |>
      dplyr::distinct(.data$patient_id, .data$value)
    counts <- table(patient_values$value)
    candidates <- as.numeric(names(counts)[counts >= minimum_patients])
    llod <- if (length(candidates)) min(candidates) else NA_real_
    values <- data[[assay]]
    observed <- is.finite(values)
    below <- if (is.na(llod)) rep(FALSE, length(values)) else
      observed & values < llod
    at <- if (is.na(llod)) rep(FALSE, length(values)) else
      observed & values == llod
    censored_percent <- if (sum(observed)) {
      100 * sum(below | at) / sum(observed)
    } else {
      NA_real_
    }
    strategy <- dplyr::case_when(
      is.na(llod) ~ "uncensored_concentration",
      sum(below) > 0L ~ "llod_conflict",
      llod <= 0 ~ "detection_only",
      censored_percent > 80 ~ "detection_only",
      censored_percent >= 20 ~ "censored_concentration_and_detection",
      TRUE ~ "censored_concentration"
    )
    tibble::tibble(
      assay = assay,
      llod = llod,
      llod_basis = paste0(
        "Lowest value observed in >=", minimum_patients, " distinct patients"
      ),
      observations = sum(observed),
      patients = dplyr::n_distinct(data$patient_id[observed]),
      patients_at_llod = if (is.na(llod)) 0L else
        dplyr::n_distinct(data$patient_id[which(at)]),
      values_below_llod = sum(below),
      values_at_llod = sum(at),
      censored_values = sum(below | at),
      censored_percent = censored_percent,
      strategy = strategy
    )
  }))
}

luminex_values_below_llod <- function(data, llod_qc) {
  dplyr::bind_rows(lapply(seq_len(nrow(llod_qc)), function(index) {
    assay_name <- llod_qc$assay[[index]]
    llod <- llod_qc$llod[[index]]
    if (is.na(llod)) return(NULL)
    tibble::tibble(
      assay = assay_name,
      patient_id = data$patient_id,
      diagnosis = data$diagnosis,
      age = data$age,
      value = data[[assay_name]],
      llod = llod
    ) |>
      dplyr::filter(!is.na(.data$value), .data$value < .data$llod)
  }))
}

luminex_llod_emmeans <- function(model) {
  suppressWarnings(
    emmeans::emmeans(model, "diagnosis") |>
      emmeans::contrast(luminex_contrast_methods(), adjust = "none") |>
      broom::tidy()
  )
}

fit_luminex_llod_concentration <- function(data, llod_qc) {
  results <- list()
  exclusions <- list()
  for (index in seq_len(nrow(llod_qc))) {
    assay <- llod_qc$assay[[index]]
    llod <- llod_qc$llod[[index]]
    strategy <- llod_qc$strategy[[index]]
    if (strategy %in% c("detection_only", "llod_conflict")) {
      exclusions[[assay]] <- tibble::tibble(
        assay = assay,
        reason = if (strategy == "llod_conflict") {
          "Observed values below the inferred LLOD"
        } else {
          "More than 80% censored or non-positive LLOD"
        }
      )
      next
    }
    model_data <- data |>
      dplyr::transmute(
        patient_id, diagnosis, age, value = .data[[assay]],
        censored = !is.na(.env$llod) & .data[[assay]] <= .env$llod
      ) |>
      dplyr::filter(
        is.finite(.data$value), .data$value > 0, !is.na(.data$age)
      ) |>
      dplyr::mutate(
        analysis_value = dplyr::if_else(
          .data$censored, .env$llod, .data$value
        )
      )
    counts <- table(model_data$diagnosis)
    if (length(counts) < 3L || any(counts < 4L) ||
        length(unique(model_data$value)) < 2L) {
      exclusions[[assay]] <- tibble::tibble(
        assay = assay, reason = "Insufficient variable observations by group"
      )
      next
    }
    fitted <- tryCatch({
      if (is.na(llod)) {
        raw_model <- stats::lm(log2(value) ~ diagnosis, data = model_data)
        adjusted_model <- stats::lm(
          log2(value) ~ diagnosis + age, data = model_data
        )
        model_type <- "linear_log2"
      } else {
        raw_model <- survival::survreg(
          survival::Surv(
            log2(analysis_value), !censored, type = "left"
          ) ~ diagnosis,
          data = model_data, dist = "gaussian", robust = TRUE,
          cluster = patient_id
        )
        adjusted_model <- survival::survreg(
          survival::Surv(
            log2(analysis_value), !censored, type = "left"
          ) ~ diagnosis + age,
          data = model_data, dist = "gaussian", robust = TRUE,
          cluster = patient_id
        )
        model_type <- "left_censored_log2"
      }
      raw <- luminex_llod_emmeans(raw_model) |>
        dplyr::select(comparison = "contrast", raw_p_value = "p.value")
      adjusted <- luminex_llod_emmeans(adjusted_model) |>
        dplyr::rename(
          comparison = "contrast",
          log2_difference = "estimate",
          age_adjusted_p_value = "p.value"
        )
      dplyr::left_join(adjusted, raw, by = "comparison") |>
        dplyr::left_join(
          dplyr::select(luminex_contrasts(), "comparison", "group1", "group2"),
          by = "comparison"
        ) |>
        dplyr::mutate(
          assay = assay, model = model_type, llod = llod,
          censored_percent = llod_qc$censored_percent[[index]], .before = 1L
        )
    }, error = function(error) {
      exclusions[[assay]] <<- tibble::tibble(
        assay = assay, reason = conditionMessage(error)
      )
      NULL
    })
    if (!is.null(fitted)) results[[assay]] <- fitted
  }
  statistics <- dplyr::bind_rows(results) |>
    dplyr::mutate(
      bh_adjusted_p_value = stats::p.adjust(
        .data$age_adjusted_p_value, method = "BH"
      )
    ) |>
    dplyr::arrange(.data$age_adjusted_p_value)
  list(statistics = statistics, exclusions = dplyr::bind_rows(exclusions))
}

luminex_llod_detection_rates <- function(data, llod_qc) {
  dplyr::bind_rows(lapply(seq_len(nrow(llod_qc)), function(index) {
    assay_name <- llod_qc$assay[[index]]
    llod <- llod_qc$llod[[index]]
    if (is.na(llod)) return(NULL)
    tibble::tibble(
      assay = assay_name, diagnosis = data$diagnosis,
      detected = dplyr::if_else(
        is.na(data[[assay_name]]), NA, data[[assay_name]] > llod
      )
    ) |>
      dplyr::filter(!is.na(.data$detected)) |>
      dplyr::group_by(.data$assay, .data$diagnosis) |>
      dplyr::summarise(
        observations = dplyr::n(),
        detected_percent = 100 * mean(.data$detected),
        detected_count = sum(.data$detected),
        censored_count = sum(!.data$detected),
        .groups = "drop"
      )
  }))
}

fit_luminex_llod_detection <- function(data, llod_qc) {
  eligible <- dplyr::filter(
    llod_qc,
    .data$strategy %in% c(
      "censored_concentration_and_detection", "detection_only"
    )
  )
  comparisons <- list(
    CIDP_vs_CTRL = c("CIDP", "CTRL"),
    GBS_vs_CTRL = c("GBS", "CTRL"),
    CIDP_vs_GBS = c("CIDP", "GBS")
  )
  results <- dplyr::bind_rows(lapply(seq_len(nrow(eligible)), function(index) {
    assay_name <- eligible$assay[[index]]
    llod <- eligible$llod[[index]]
    dplyr::bind_rows(lapply(names(comparisons), function(comparison) {
      groups <- comparisons[[comparison]]
      test_data <- tibble::tibble(
        diagnosis = data$diagnosis, value = data[[assay_name]]
      ) |>
        dplyr::filter(
          .data$diagnosis %in% .env$groups, !is.na(.data$value)
        ) |>
        dplyr::mutate(detected = .data$value > .env$llod)
      contingency <- table(
        droplevels(test_data$diagnosis), test_data$detected
      )
      if (nrow(contingency) < 2L || ncol(contingency) < 2L) return(NULL)
      test <- stats::fisher.test(contingency)
      tibble::tibble(
        assay = assay_name, comparison = comparison, llod = llod,
        detection_odds_ratio = unname(test$estimate),
        detection_p_value = test$p.value
      )
    }))
  }))
  dplyr::mutate(
    results,
    bh_adjusted_detection_p_value = stats::p.adjust(
      .data$detection_p_value, method = "BH"
    )
  ) |>
    dplyr::arrange(.data$detection_p_value)
}

analyze_luminex_llod <- function(input) {
  llod_qc <- infer_luminex_llod(input$data, input$assays)
  concentration <- fit_luminex_llod_concentration(input$data, llod_qc)
  list(
    data_wide = input$data,
    assays = input$assays,
    llod_qc = llod_qc,
    below_llod_values = luminex_values_below_llod(input$data, llod_qc),
    concentration_statistics = concentration$statistics,
    concentration_exclusions = concentration$exclusions,
    detection_rates = luminex_llod_detection_rates(input$data, llod_qc),
    detection_statistics = fit_luminex_llod_detection(input$data, llod_qc)
  )
}

luminex_llod_boxplot <- function(result) {
  colors <- stats::setNames(
    pals::cols25(nlevels(result$data_wide$diagnosis)),
    levels(result$data_wide$diagnosis)
  )
  plots <- lapply(gtools::mixedsort(result$assays), function(assay) {
    qc <- dplyr::filter(result$llod_qc, .data$assay == .env$assay)
    llod <- qc$llod[[1L]]
    conflict <- qc$strategy[[1L]] == "llod_conflict"
    plot_data <- result$data_wide |>
      dplyr::transmute(
        diagnosis, value = .data[[assay]],
        measurement = dplyr::case_when(
          .env$conflict & !is.na(.data$value) & .data$value <= .env$llod ~
            "LLOD conflict",
          !is.na(.env$llod) & .data$value <= .env$llod ~ "<= LLOD",
          TRUE ~ "Detected"
        )
      )
    ggplot2::ggplot(
      plot_data, ggplot2::aes(x = diagnosis, y = value, fill = diagnosis)
    ) +
      ggplot2::geom_boxplot(na.rm = TRUE, outlier.shape = NA) +
      ggplot2::geom_jitter(
        ggplot2::aes(shape = measurement), width = 0.2, height = 0,
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_shape_manual(values = c(
        "<= LLOD" = 25, Detected = 21, "LLOD conflict" = 24
      )) +
      ggplot2::labs(
        title = paste0(
          assay, " (", round(qc$censored_percent[[1L]], 1), "% censored)"
        ),
        x = NULL, y = "pg/ml", shape = NULL
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
      concentration_statistics = result$concentration_statistics,
      detection_statistics = result$detection_statistics,
      detection_rates = result$detection_rates,
      llod_qc = result$llod_qc,
      below_llod_values = result$below_llod_values,
      concentration_exclusions = result$concentration_exclusions
    ),
    paths[["statistics"]]
  )
  ggplot2::ggsave(
    paths[["boxplots"]], luminex_llod_boxplot(result),
    width = 12, height = 32, limitsize = FALSE
  )
  unname(paths)
}

write_luminex_llod_volcano <- function(
  result, comparison, group1, group2, seed
) {
  data <- result$concentration_statistics |>
    dplyr::filter(.data$comparison == .env$comparison) |>
    dplyr::mutate(neg_log10_q = -log10(.data$bh_adjusted_p_value))
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = log2_difference, y = neg_log10_q, label = assay,
      color = bh_adjusted_p_value < 0.1
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(
      yintercept = -log10(0.1), color = "blue", linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-0.5, 0, 0.5), color = "red",
      linetype = c("dashed", "solid", "dashed")
    ) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(
        data, .data$bh_adjusted_p_value < 0.1,
        abs(.data$log2_difference) >= 0.5
      ),
      seed = seed
    ) +
    ggplot2::scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue")) +
    ggplot2::labs(
      x = paste0("Estimated log2 difference (", group1, " - ", group2, ")"),
      y = expression(-Log[10] ~ "BH-adjusted p value"), color = "BH < 0.1"
    ) +
    ggplot2::theme_classic()
  path <- file.path(
    luminex_llod_result_dir(), paste0("volcano_", comparison, "_llod.pdf")
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = 4, height = 4)
  path
}
