luminex_result_dir <- function() {
  file.path("results", "targets", "luminex")
}

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
  raw <- readxl::read_xlsx(path, col_types = "text", na = c("", "NA"))
  required <- c("ID_seq", "age", "diagnosis", "IL_1a")
  stopifnot(all(required %in% names(raw)))

  assay_start <- match("IL_1a", names(raw))
  assays <- names(raw)[seq.int(assay_start, ncol(raw))]
  luminex_validate_numeric(raw, assays)
  data <- raw |>
    dplyr::filter(!is.na(.data$ID_seq), .data$diagnosis %in% c("CIDP", "GBS")) |>
    dplyr::transmute(
      patient_id = as.character(.data$ID_seq),
      age = as.numeric(.data$age),
      diagnosis = as.character(.data$diagnosis),
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
  list(data = data, assays = assays)
}

read_luminex_controls <- function(path) {
  raw <- readxl::read_xlsx(
    path, col_types = "text", na = c("", "NA", "-")
  )
  required <- c("biobank_id", "age", "Group...8", "IL_1a")
  stopifnot(all(required %in% names(raw)))
  assays <- names(raw)[seq.int(match("IL_1a", names(raw)), ncol(raw))]
  luminex_validate_numeric(raw, assays)
  data <- raw |>
    dplyr::filter(.data$Group...8 == "control") |>
    dplyr::transmute(
      patient_id = as.character(.data$biobank_id),
      age = as.numeric(.data$age),
      diagnosis = "CTRL",
      dplyr::across(tidyselect::all_of(assays), luminex_numeric_value)
    )
  stopifnot(nrow(data) > 0L, !anyNA(data[c("patient_id", "age")]))
  list(data = data, assays = assays)
}

combine_luminex_inputs <- function(disease, controls) {
  stopifnot(setequal(disease$assays, controls$assays))
  disease_data <- dplyr::mutate(
    disease$data, patient_id = paste0("disease_", .data$patient_id)
  )
  control_data <- dplyr::mutate(
    controls$data, patient_id = paste0("control_", .data$patient_id)
  )
  data <- dplyr::bind_rows(disease_data, control_data)
  data$diagnosis <- factor(data$diagnosis, levels = c("CTRL", "GBS", "CIDP"))
  stopifnot(all(table(data$diagnosis) >= 4L))
  list(data = data, assays = disease$assays)
}

luminex_contrasts <- function() {
  tibble::tribble(
    ~comparison, ~group1, ~group2, ~seed,
    "CIDP_vs_CTRL", "CIDP", "CTRL", 500L,
    "GBS_vs_CTRL", "GBS", "CTRL", 600L,
    "CIDP_vs_GBS", "CIDP", "GBS", 700L
  )
}

luminex_contrast_methods <- function() {
  list(
    CIDP_vs_CTRL = c(-1, 0, 1),
    GBS_vs_CTRL = c(-1, 1, 0),
    CIDP_vs_GBS = c(0, -1, 1)
  )
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
    raw <- emmeans::emmeans(raw_model, "diagnosis") |>
      emmeans::contrast(luminex_contrast_methods(), adjust = "none") |>
      broom::tidy() |>
      dplyr::select(comparison = "contrast", raw_p_value = "p.value")
    adjusted <- emmeans::emmeans(adjusted_model, "diagnosis") |>
      emmeans::contrast(luminex_contrast_methods(), adjust = "none") |>
      broom::tidy() |>
      dplyr::rename(
        comparison = "contrast", age_adjusted_p_value = "p.value"
      )
    dplyr::left_join(adjusted, raw, by = "comparison") |>
      dplyr::left_join(
        dplyr::select(luminex_contrasts(), "comparison", "group1", "group2"),
        by = "comparison"
      ) |>
      dplyr::mutate(var = assay, .before = 1L)
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
    group = "diagnosis",
    unit = "pg/ml"
  )
}

luminex_boxplot <- function(result) {
  colors <- stats::setNames(
    pals::cols25(nlevels(result$data_wide[[result$group]])),
    levels(result$data_wide[[result$group]])
  )
  plots <- lapply(gtools::mixedsort(result$assays), function(assay) {
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
      ggplot2::geom_jitter(width = 0.2, height = 0, na.rm = TRUE) +
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
      excluded_assays = result$excluded_assays
    ),
    paths[["stats"]]
  )
  ggplot2::ggsave(
    paths[["boxplot"]], luminex_boxplot(result), width = 11, height = 24,
    limitsize = FALSE
  )
  unname(paths)
}

luminex_volcano_data <- function(
  result, comparison, group1, group2, seed = 500L, fdr = 0.1,
  permutations = 10000L
) {
  tested <- dplyr::filter(
    result$stats,
    .data$comparison == .env$comparison,
    is.finite(.data$age_adjusted_p_value)
  ) |>
    dplyr::transmute(var = .data$var, p.value = .data$age_adjusted_p_value)
  data <- result$data_wide |>
    dplyr::filter(.data[[result$group]] %in% c(group1, group2)) |>
    droplevels()
  means <- data |>
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
    dplyr::mutate(log2_ratio = log2(.data[[group1]] / .data[[group2]]))

  matrix <- t(as.matrix(data[tested$var]))
  complete <- stats::complete.cases(matrix)
  pooled_variance <- vapply(tested$var, function(assay) {
    values <- data[[assay]]
    sum(vapply(c(group1, group2), function(group) {
      stats::var(values[data[[result$group]] == group], na.rm = TRUE)
    }, numeric(1)), na.rm = TRUE)
  }, numeric(1))
  permutation_rows <- complete & is.finite(pooled_variance) & pooled_variance > 0
  stopifnot(any(permutation_rows), permutations > 0L)
  design <- ifelse(data[[result$group]] == group1, 1L, 2L)
  threshold <- withr::with_seed(seed, permFDP::permFDP.adjust.threshold(
    pVals = tested$p.value[permutation_rows],
    threshold = fdr,
    myDesign = design,
    intOnly = matrix[permutation_rows, , drop = FALSE],
    nPerms = permutations
  ))
  dplyr::left_join(means, tested, by = "var") |>
    dplyr::mutate(
      p.adj.threshold = threshold,
      significant = .data$p.value < threshold,
      neg_log10_p = -log10(.data$p.value)
    )
}

write_luminex_volcano <- function(
  data, comparison, group1, group2, seed = 500L
) {
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
      x = paste0("Log2 fold change (", group1, " / ", group2, ")"),
      y = expression(-Log[10] ~ "p value")
    ) +
    ggplot2::theme_classic()
  path <- file.path(
    luminex_result_dir(), paste0("volcano_", comparison, ".pdf")
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = 3, height = 3)
  path
}
