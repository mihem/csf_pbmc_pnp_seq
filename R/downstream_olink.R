olink_result_dir <- function() {
  file.path("results", "targets", "olink")
}

olink_assays <- function() {
  c(
    "CCL2", "CCL3", "CXCL8", "GZMA", "KLRD1", "SIRT2", "TNFSF14",
    "TREM2", "VEGFA", "IFNG", "CXCL12", "CXCL10", "TNFRSF9", "CCL7",
    "TNFSF12", "IL1B", "IL18"
  )
}

read_olink_input <- function(path, value_column, assays = olink_assays()) {
  result <- readxl::read_xlsx(path)
  stopifnot(all(c("SampleID", "Assay", value_column) %in% names(result)))
  result[[value_column]] <- as.numeric(result[[value_column]])
  result <- dplyr::filter(result, .data$Assay %in% .env$assays)
  stopifnot(nrow(result) > 0L, all(unique(result$Assay) %in% assays))
  result
}

read_olink_metadata <- function(path) {
  result <- readxl::read_xlsx(path)
  required <- c("SampleID", "orbis_id", "diagnosis", "group2", "age", "sex")
  stopifnot(all(required %in% names(result)), !anyDuplicated(result$SampleID))
  result
}

olink_group_configuration <- function(group) {
  switch(
    group,
    diagnosis = list(levels = c("CTRL", "GBS", "CIDP"), suffix = ""),
    group2 = list(levels = c("CTRL", "IN"), suffix = "_meta"),
    stop("Unknown Olink group: ", group)
  )
}

fit_olink_assays <- function(data, assays, group) {
  results <- lapply(assays, function(assay) {
    counts <- data |>
      dplyr::filter(!is.na(.data[[assay]])) |>
      dplyr::count(.data[[group]]) |>
      dplyr::filter(.data$n >= 4L)
    if (nrow(counts) < nlevels(data[[group]])) {
      return(NULL)
    }
    model <- lme4::lmer(
      stats::reformulate(
        c(group, "sex", "age", "(1 | orbis_id)"), response = assay
      ),
      data = data
    )
    result <- emmeans::emmeans(model, group, adjust = "none") |>
      graphics::pairs(adjust = "none") |>
      broom::tidy()
    separated <- tidyr::separate_wider_delim(
      result, "contrast", delim = " - ", names = c("group1", "group2")
    )
    dplyr::mutate(separated, var = assay, .before = 1L)
  })
  dplyr::bind_rows(results) |>
    dplyr::mutate(
      p.adj = stats::p.adjust(.data$p.value, method = "BH"),
      p.adj.signif = as.character(stats::symnum(
        .data$p.adj,
        corr = FALSE,
        na = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.1, 1),
        symbols = c("***", "**", "*", " ")
      ))
    )
}

analyze_olink_data <- function(
  data, metadata, value_column, group, unit, seed = 123L
) {
  configuration <- olink_group_configuration(group)
  joined <- dplyr::left_join(data, metadata, by = "SampleID")
  stopifnot(!anyNA(joined$orbis_id), !anyNA(joined[[group]]))
  wide <- tidyr::pivot_wider(
    joined,
    id_cols = tidyselect::all_of(c(
      "SampleID", "orbis_id", group, "age", "sex"
    )),
    names_from = "Assay",
    values_from = tidyselect::all_of(value_column)
  )
  wide[[group]] <- factor(wide[[group]], levels = configuration$levels)
  assays <- sort(intersect(unique(joined$Assay), names(wide)))
  stopifnot(length(assays) > 0L)
  stats <- withr::with_seed(seed, fit_olink_assays(wide, assays, group))
  list(
    data_wide = wide,
    stats = stats,
    assays = assays,
    group = group,
    unit = unit,
    suffix = configuration$suffix
  )
}

olink_boxplot <- function(result) {
  colors <- stats::setNames(pals::cols25(nlevels(result$data_wide[[result$group]])),
                            levels(result$data_wide[[result$group]]))
  plots <- lapply(result$assays, function(assay) {
    annotations <- dplyr::filter(
      result$stats, .data$var == .env$assay, .data$p.adj < 0.1
    )
    plot <- ggplot2::ggplot(
      result$data_wide,
      ggplot2::aes(
        x = .data[[result$group]], y = .data[[assay]],
        fill = .data[[result$group]]
      )
    ) +
      ggplot2::geom_boxplot() +
      ggplot2::geom_jitter(width = 0.2) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(title = paste0(assay, " (", result$unit, ")"), x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
    if (nrow(annotations)) {
      plot <- plot + ggsignif::geom_signif(
        comparisons = Map(c, annotations$group1, annotations$group2),
        annotations = annotations$p.adj.signif,
        step_increase = 0.05,
        vjust = 0.7
      )
    }
    plot
  })
  patchwork::wrap_plots(plots)
}

olink_artifact_paths <- function(type, suffix = "") {
  stopifnot(type %in% c("quant", "npx"))
  root <- olink_result_dir()
  if (type == "quant") {
    c(
      object = file.path(root, paste0("olink_quant", suffix, ".qs")),
      stats = file.path(root, paste0("olink_quant_stats_lme4", suffix, ".xlsx")),
      boxplot = file.path(root, paste0("olink_quant_boxplot_lme4", suffix, ".pdf"))
    )
  } else {
    c(
      object = file.path(root, paste0("olink_npx_wide", suffix, ".qs")),
      stats = file.path(root, paste0("olink_npx_stats_all", suffix, ".xlsx")),
      boxplot = file.path(root, paste0("olink_npx_boxplots", suffix, ".pdf"))
    )
  }
}

write_olink_artifacts <- function(result, type) {
  paths <- olink_artifact_paths(type, result$suffix)
  dir.create(olink_result_dir(), recursive = TRUE, showWarnings = FALSE)
  qs::qsave(if (type == "quant") result else result$data_wide, paths[["object"]])
  writexl::write_xlsx(result$stats, paths[["stats"]])
  ggplot2::ggsave(paths[["boxplot"]], olink_boxplot(result), width = 8, height = 12)
  unname(paths)
}

olink_contrasts <- function() {
  tibble::tribble(
    ~comparison, ~group1, ~group2,
    "CIDP_vs_CTRL", "CIDP", "CTRL",
    "GBS_vs_CTRL", "GBS", "CTRL",
    "CIDP_vs_GBS", "CIDP", "GBS"
  )
}

olink_volcano_data <- function(
  result, group1, group2, seed = 123L, fdr = 0.1, permutations = 10000L
) {
  data <- dplyr::filter(
    result$data_wide, .data[[result$group]] %in% c(group1, group2)
  ) |>
    droplevels()
  tested <- lapply(result$assays, function(assay) {
    counts <- data |>
      dplyr::filter(!is.na(.data[[assay]])) |>
      dplyr::count(.data[[result$group]]) |>
      dplyr::filter(.data$n >= 4L)
    if (nrow(counts) < 2L) return(NULL)
    model <- lme4::lmer(
      stats::reformulate(
        c(result$group, "sex", "age", "(1 | orbis_id)"), response = assay
      ),
      data = data
    )
    contrast <- emmeans::emmeans(model, result$group, adjust = "none") |>
      graphics::pairs(adjust = "none") |>
      broom::tidy()
    tibble::tibble(var = assay, p.value = contrast$p.value[[1L]])
  }) |>
    dplyr::bind_rows()
  means <- data |>
    dplyr::group_by(.data[[result$group]]) |>
    dplyr::summarise(dplyr::across(
      tidyselect::all_of(tested$var), \(value) mean(value, na.rm = TRUE)
    ),
                     .groups = "drop") |>
    tidyr::pivot_longer(-tidyselect::all_of(result$group), names_to = "var") |>
    tidyr::pivot_wider(names_from = tidyselect::all_of(result$group), values_from = "value") |>
    dplyr::mutate(log2_ratio = log2(.data[[group1]] / .data[[group2]]))
  matrix <- t(as.matrix(data[tested$var]))
  complete <- stats::complete.cases(matrix)
  stopifnot(any(complete), permutations > 0L)
  design <- ifelse(data[[result$group]] == levels(data[[result$group]])[[1L]], 1L, 2L)
  threshold <- withr::with_seed(seed, permFDP::permFDP.adjust.threshold(
    pVals = tested$p.value[complete], threshold = fdr,
    myDesign = design, intOnly = matrix[complete, , drop = FALSE],
    nPerms = permutations
  ))
  dplyr::left_join(means, tested, by = "var") |>
    dplyr::mutate(
      p.adj.threshold = threshold,
      significant = .data$p.value < threshold,
      neg_log10_p = -log10(.data$p.value)
    )
}

write_olink_volcano <- function(data, comparison, type, seed = 123L) {
  threshold <- -log10(data$p.adj.threshold[[1L]])
  plot <- ggplot2::ggplot(
    data, ggplot2::aes(x = log2_ratio, y = neg_log10_p, label = var)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_hline(yintercept = threshold, color = "blue", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-0.5, 0, 0.5),
                        color = "red", linetype = c("dashed", "solid", "dashed")) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(data, .data$neg_log10_p >= threshold,
                           abs(.data$log2_ratio) >= 0.5),
      seed = seed
    ) +
    ggplot2::labs(x = expression(Log[2] ~ "fold change"),
                  y = expression(-Log[10] ~ "p value")) +
    ggplot2::theme_classic()
  path <- file.path(olink_result_dir(), paste0("volcano_", comparison, "_", type, ".pdf"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = 3, height = 3)
  path
}

validate_olink_stats <- function(current_files, legacy_files) {
  current_stats <- current_files[grepl("stats.*xlsx$", current_files)]
  stopifnot(length(current_stats) == length(legacy_files))
  checks <- lapply(legacy_files, function(legacy) {
    current <- current_stats[basename(current_stats) == basename(legacy)]
    stopifnot(length(current) == 1L)
    observed <- readxl::read_xlsx(current)
    expected <- readxl::read_xlsx(legacy)
    stopifnot(
      identical(names(observed), names(expected)),
      nrow(observed) == nrow(expected),
      isTRUE(all.equal(observed, expected, check.attributes = FALSE, tolerance = 1e-8))
    )
    tibble::tibble(artifact = basename(legacy), passed = TRUE)
  })
  dplyr::bind_rows(checks)
}
