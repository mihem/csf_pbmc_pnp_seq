luminex_result_dir <- function() {
  file.path("results", "targets", "luminex")
}

luminex_numeric_value <- function(value) {
  numeric <- as.numeric(value)
  is_excel_date <- !is.na(numeric) & numeric >= 40000 & numeric <= 60000
  date <- as.Date(numeric[is_excel_date], origin = "1899-12-30")
  numeric[is_excel_date] <- as.integer(format(date, "%d")) +
    as.integer(format(date, "%m")) / 100
  numeric
}

read_luminex_input <- function(path) {
  raw <- readxl::read_xlsx(path, col_types = "text", na = c("", "NA"))
  required <- c("ID_seq", "age", "diagnosis", "IL_1a")
  stopifnot(all(required %in% names(raw)))

  assay_start <- match("IL_1a", names(raw))
  assays <- names(raw)[seq.int(assay_start, ncol(raw))]
  date_corrections <- tibble::tibble(
    assay = assays,
    corrected_values = vapply(raw[assays], function(value) {
      numeric <- suppressWarnings(as.numeric(value))
      sum(!is.na(numeric) & numeric >= 40000 & numeric <= 60000)
    }, integer(1))
  ) |>
    dplyr::filter(.data$corrected_values > 0L)
  data <- raw |>
    dplyr::filter(!is.na(.data$ID_seq), .data$diagnosis %in% c("CIDP", "GBS")) |>
    dplyr::transmute(
      patient_id = as.character(.data$ID_seq),
      age = as.numeric(.data$age),
      diagnosis = factor(.data$diagnosis, levels = c("CIDP", "GBS")),
      dplyr::across(tidyselect::all_of(assays), luminex_numeric_value)
    )

  patient_diagnoses <- data |>
    dplyr::distinct(.data$patient_id, .data$diagnosis) |>
    dplyr::count(.data$patient_id)
  stopifnot(
    nrow(data) > 0L,
    !anyNA(data$patient_id),
    all(patient_diagnoses$n == 1L),
    all(table(data$diagnosis) >= 4L),
    length(assays) > 0L
  )
  list(data = data, assays = assays, date_corrections = date_corrections)
}

fit_luminex_assays <- function(data, assays) {
  results <- lapply(assays, function(assay) {
    model_data <- data |>
      dplyr::filter(
        !is.na(.data[[assay]]), !is.na(.data$diagnosis), !is.na(.data$age)
      )
    counts <- table(model_data$diagnosis)
    if (length(unique(model_data[[assay]])) < 2L || any(counts < 4L)) {
      return(NULL)
    }

    raw_model <- stats::lm(
      stats::reformulate("diagnosis", response = assay),
      data = model_data
    )
    adjusted_model <- stats::lm(
      stats::reformulate(c("diagnosis", "age"), response = assay),
      data = model_data
    )
    raw_p_value <- emmeans::emmeans(raw_model, "diagnosis", adjust = "none") |>
      graphics::pairs(adjust = "none") |>
      broom::tidy() |>
      dplyr::pull("p.value")
    result <- emmeans::emmeans(adjusted_model, "diagnosis", adjust = "none") |>
      graphics::pairs(adjust = "none") |>
      broom::tidy()
    separated <- tidyr::separate_wider_delim(
      result, "contrast", delim = " - ", names = c("group1", "group2")
    )
    dplyr::mutate(
      separated,
      var = assay,
      raw_p_value = raw_p_value,
      age_adjusted_p_value = .data$p.value,
      .before = 1L
    ) |>
      dplyr::select(-"p.value")
  })

  dplyr::bind_rows(results) |>
    dplyr::mutate(
      bh_adjusted_p_value = stats::p.adjust(
        .data$age_adjusted_p_value, method = "BH"
      ),
      bh_signif = as.character(stats::symnum(
        .data$bh_adjusted_p_value,
        corr = FALSE,
        na = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.1, 1),
        symbols = c("***", "**", "*", " ")
      ))
    ) |>
    dplyr::arrange(.data$age_adjusted_p_value)
}

analyze_luminex_data <- function(input, seed = 400L) {
  stats <- withr::with_seed(
    seed, fit_luminex_assays(input$data, input$assays)
  )
  tested <- unique(stats$var)
  excluded <- tibble::tibble(
    assay = setdiff(input$assays, tested),
    reason = "Invariant or insufficient complete observations"
  )
  list(
    data_wide = input$data,
    stats = stats,
    assays = tested,
    excluded_assays = excluded,
    date_corrections = input$date_corrections,
    group = "diagnosis",
    unit = "pg/ml"
  )
}

luminex_boxplot <- function(result) {
  colors <- stats::setNames(
    pals::cols25(nlevels(result$data_wide[[result$group]])),
    levels(result$data_wide[[result$group]])
  )
  plots <- lapply(result$assays, function(assay) {
    annotations <- dplyr::filter(
      result$stats,
      .data$var == .env$assay,
      .data$bh_adjusted_p_value < 0.1
    )
    plot <- ggplot2::ggplot(
      result$data_wide,
      ggplot2::aes(
        x = .data[[result$group]], y = .data[[assay]],
        fill = .data[[result$group]]
      )
    ) +
      ggplot2::geom_boxplot(na.rm = TRUE) +
      ggplot2::geom_jitter(width = 0.2, na.rm = TRUE) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(
        title = paste0(assay, " (", result$unit, ")"), x = NULL, y = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
    if (nrow(annotations)) {
      plot <- plot + ggsignif::geom_signif(
        comparisons = Map(c, annotations$group1, annotations$group2),
        annotations = annotations$bh_signif,
        step_increase = 0.05,
        vjust = 0.7
      )
    }
    plot
  })
  patchwork::wrap_plots(plots)
}

write_luminex_artifacts <- function(result) {
  root <- luminex_result_dir()
  paths <- c(
    object = file.path(root, "luminex.qs"),
    stats = file.path(root, "luminex_stats_age_adjusted.xlsx"),
    boxplot = file.path(root, "luminex_boxplots_age_adjusted.pdf")
  )
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  qs::qsave(result, paths[["object"]])
  writexl::write_xlsx(
    list(
      statistics = result$stats,
      excluded_assays = result$excluded_assays,
      excel_date_corrections = result$date_corrections
    ),
    paths[["stats"]]
  )
  ggplot2::ggsave(
    paths[["boxplot"]], luminex_boxplot(result), width = 10, height = 24,
    limitsize = FALSE
  )
  unname(paths)
}

luminex_volcano_data <- function(
  result, seed = 500L, fdr = 0.1, permutations = 10000L
) {
  tested <- dplyr::filter(
    result$stats, is.finite(.data$age_adjusted_p_value)
  ) |>
    dplyr::transmute(var = .data$var, p.value = .data$age_adjusted_p_value)
  means <- result$data_wide |>
    dplyr::group_by(.data[[result$group]]) |>
    dplyr::summarise(
      dplyr::across(
        tidyselect::all_of(tested$var), \(value) mean(value, na.rm = TRUE)
      ),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(-tidyselect::all_of(result$group), names_to = "var") |>
    tidyr::pivot_wider(
      names_from = tidyselect::all_of(result$group), values_from = "value"
    ) |>
    dplyr::mutate(log2_ratio = log2(.data$CIDP / .data$GBS))

  matrix <- t(as.matrix(result$data_wide[tested$var]))
  complete <- stats::complete.cases(matrix)
  stopifnot(any(complete), permutations > 0L)
  design <- ifelse(result$data_wide[[result$group]] == "CIDP", 1L, 2L)
  threshold <- withr::with_seed(seed, permFDP::permFDP.adjust.threshold(
    pVals = tested$p.value[complete],
    threshold = fdr,
    myDesign = design,
    intOnly = matrix[complete, , drop = FALSE],
    nPerms = permutations
  ))
  dplyr::left_join(means, tested, by = "var") |>
    dplyr::mutate(
      p.adj.threshold = threshold,
      significant = .data$p.value < threshold,
      neg_log10_p = -log10(.data$p.value)
    )
}

write_luminex_volcano <- function(data, seed = 500L) {
  threshold <- -log10(data$p.adj.threshold[[1L]])
  plot <- ggplot2::ggplot(
    data, ggplot2::aes(x = log2_ratio, y = neg_log10_p, label = var)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(
      yintercept = threshold, color = "blue", linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-0.5, 0, 0.5), color = "red",
      linetype = c("dashed", "solid", "dashed")
    ) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(
        data, .data$neg_log10_p >= threshold, abs(.data$log2_ratio) >= 0.5
      ),
      seed = seed
    ) +
    ggplot2::labs(
      x = expression(Log[2] ~ "fold change (CIDP / GBS)"),
      y = expression(-Log[10] ~ "p value")
    ) +
    ggplot2::theme_classic()
  path <- file.path(luminex_result_dir(), "volcano_CIDP_vs_GBS.pdf")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = 3, height = 3)
  path
}
