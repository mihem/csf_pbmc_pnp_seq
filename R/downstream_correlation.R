correlation_result_dir <- function() {
  file.path("results", "targets", "correlation")
}

correlation_scores <- function() {
  c(
    "incat_at_lumbar_puncture", "incat_follow_up", "incat_progress",
    "onls_at_lumbar_puncture", "onls_follow_up", "onls_progress",
    "mrc_sum_score_60_at_lumbar_puncture", "mrc_sum_score_60_follow_up",
    "mrc_sum_score_60_progress", "csf_protein", "dml_ulnar_motoric",
    "cmap_ulnar_motoric", "ncv_ulnar_motoric",
    "min_f_latency_ulnar_motoric", "cmap_tibial_motoric",
    "ncv_tibial_motoric"
  )
}

correlation_jobs <- function() {
  tidyr::crossing(
    score = correlation_scores(),
    tissue = c("CSF", "PBMC")
  ) |>
    dplyr::mutate(job = paste(.data$score, tolower(.data$tissue), sep = "_"))
}

write_severity_subgroup_plot <- function(lookup, seed = 123L) {
  required <- c("patient", "diagnosis", "incat_at_lumbar_puncture")
  stopifnot(all(required %in% names(lookup)))
  high_nk <- c("P16", "P17", "P18", "P21", "P23")
  data <- dplyr::mutate(
    lookup,
    subgroup = dplyr::if_else(
      .data$patient %in% .env$high_nk, "NKhi", as.character(.data$diagnosis)
    )
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = subgroup, y = incat_at_lumbar_puncture, fill = subgroup)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")
  path <- file.path(correlation_result_dir(), "severity_subgroup.pdf")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  withr::with_seed(seed, ggplot2::ggsave(path, plot, width = 6, height = 4))
  path
}

run_abundance_clinical_correlation <- function(object, lookup, score, tissue) {
  required_meta <- c("cluster", "patient", "tissue")
  required_lookup <- c("patient", "sex", "age", score)
  stopifnot(
    inherits(object, "Seurat"),
    all(required_meta %in% names(object[[]])),
    all(required_lookup %in% names(lookup)),
    !anyDuplicated(lookup$patient),
    tissue %in% c("CSF", "PBMC")
  )
  metadata <- object[[]] |>
    dplyr::select(-tidyselect::any_of(score)) |>
    dplyr::left_join(
      dplyr::select(lookup, "patient", tidyselect::all_of(score)),
      by = "patient"
    ) |>
    dplyr::filter(.data$tissue == .env$tissue, !is.na(.data[[score]]))
  stopifnot(dplyr::n_distinct(metadata$patient) >= 5L)
  props <- speckle::getTransformedProps(
    clusters = metadata$cluster,
    sample = metadata$patient,
    transform = "logit"
  )$TransformedProps
  model_data <- tibble::tibble(patient = colnames(props)) |>
    dplyr::left_join(lookup, by = "patient") |>
    dplyr::distinct(.data$patient, .keep_all = TRUE)
  stopifnot(
    nrow(model_data) == ncol(props),
    !anyNA(model_data[[score]]), !anyNA(model_data$sex), !anyNA(model_data$age)
  )
  design <- stats::model.matrix(
    stats::as.formula(paste0("~0 + ", score, " + sex + age")),
    data = model_data
  )
  fit <- limma::lmFit(props, design) |>
    limma::eBayes()
  coefficient <- which(colnames(design) == score)
  stopifnot(length(coefficient) == 1L)
  degrees <- fit$df.residual[[1L]]
  statistics <- tibble::tibble(
    cell_type = rownames(fit$coefficients),
    limma_coef = fit$coefficients[, coefficient],
    limma_t_stat = fit$t[, coefficient],
    limma_p_value = fit$p.value[, coefficient],
    limma_sigma = fit$sigma,
    adjusted_pearson_r = fit$t[, coefficient] /
      sqrt(fit$t[, coefficient]^2 + degrees)
  ) |>
    dplyr::mutate(
      limma_p_adj = stats::p.adjust(.data$limma_p_value, method = "BH")
    )
  list(
    statistics = statistics,
    fit = fit,
    proportions = props,
    metadata = model_data,
    score = score,
    tissue = tissue
  )
}

correlation_plot <- function(result, correlation_threshold = 0.3,
                             significance_threshold = 0.05) {
  selected <- result$statistics |>
    dplyr::filter(
      abs(.data$adjusted_pearson_r) > correlation_threshold,
      .data$limma_p_adj < significance_threshold
    ) |>
    dplyr::arrange(dplyr::desc(abs(.data$adjusted_pearson_r)))
  if (!nrow(selected)) return(NULL)
  fitted <- result$fit$design %*% t(result$fit$coefficients)
  residuals <- result$proportions - t(fitted)
  plots <- lapply(selected$cell_type, function(cell_type) {
    stats <- dplyr::filter(selected, .data$cell_type == .env$cell_type)
    partial <- residuals[cell_type, ] +
      result$fit$coefficients[cell_type, result$score] *
      result$metadata[[result$score]]
    data <- tibble::tibble(
      score = result$metadata[[result$score]], partial_residual = partial
    ) |>
      tidyr::drop_na()
    ggplot2::ggplot(
      data, ggplot2::aes(x = score, y = partial_residual)
    ) +
      ggplot2::geom_point(alpha = 0.7, size = 1.5) +
      ggplot2::geom_smooth(method = "lm", se = TRUE, color = "blue") +
      ggplot2::labs(
        title = gsub("_", " ", cell_type),
        subtitle = sprintf("r = %.3f, adj. p = %.3g",
                           stats$adjusted_pearson_r, stats$limma_p_adj),
        x = tools::toTitleCase(gsub("_", " ", result$score)),
        y = "Proportion\n(adj. for age & sex)"
      ) +
      ggplot2::theme_classic()
  })
  patchwork::wrap_plots(plots, ncol = 3L)
}

write_correlation_result <- function(result) {
  root <- correlation_result_dir()
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0("correlation_", result$score, "_", result$tissue)
  csv <- file.path(root, paste0(stem, ".csv"))
  readr::write_csv(result$statistics, csv)
  plot <- correlation_plot(result)
  if (is.null(plot)) return(csv)
  pdf <- file.path(root, paste0(stem, ".pdf"))
  n <- sum(
    abs(result$statistics$adjusted_pearson_r) > 0.3 &
      result$statistics$limma_p_adj < 0.05,
    na.rm = TRUE
  )
  ggplot2::ggsave(pdf, plot, width = 12, height = 4 * ceiling(n / 3))
  c(csv, pdf)
}

run_rcna_association <- function(object, lookup, score, seed = 123L) {
  stopifnot(
    inherits(object, "Seurat"),
    all(c("patient", "sex", "age", "sample") %in% names(object[[]])),
    all(c("patient", score) %in% names(lookup)),
    !anyDuplicated(lookup$patient),
    "RNA_nn" %in% names(object@graphs),
    "umap.stacas.ss.all" %in% names(object@reductions)
  )
  object[[score]] <- lookup[[score]][match(object$patient, lookup$patient)]
  keep <- !is.na(object[[score, drop = TRUE]])
  stopifnot(sum(keep) > 0L)
  filtered <- object[, keep]
  filtered$sex_numeric <- as.numeric(filtered$sex == "male")
  # rcna 0.1.0 calls glue() internally without importing the package.
  withr::local_package("glue")
  tryCatch(
    withr::with_seed(seed, {
      rcna_data <- rcna:::create_object.Seurat(
        filtered,
        samplem_key = "sample",
        samplem_vars = c(score, "sex_numeric", "age"),
        graph_use = "RNA_nn"
      )
      association <- rcna::association(
        data = rcna_data,
        y = rcna_data$samplem[[score]],
        batches = NULL,
        covs = c("sex_numeric", "age"),
        return_nam = FALSE,
        verbose = FALSE
      )
      filtered$cna_ncorrs <- association$ncorrs[
        colnames(filtered), , drop = TRUE
      ]
      filtered
    }),
    error = function(condition) {
      structure(
        list(score = score, error = conditionMessage(condition)),
        class = "rcna_failure"
      )
    }
  )
}

write_rcna_plot <- function(object, score) {
  if (inherits(object, "rcna_failure")) {
    path <- file.path(correlation_result_dir(), paste0("rcna_", score, ".txt"))
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines(
      c(
        paste("score:", score),
        "status: failed",
        paste("error:", object$error)
      ),
      path
    )
    return(path)
  }
  plot <- Seurat::FeaturePlot(
    object,
    features = "cna_ncorrs",
    reduction = "umap.stacas.ss.all",
    pt.size = 0.1,
    order = FALSE,
    coord.fixed = TRUE,
    raster = FALSE,
    alpha = 0.2
  ) +
    viridis::scale_color_viridis(option = "inferno") +
    ggplot2::labs(
      title = gsub("_", " ", score), color = "Correlation",
      x = "UMAP1", y = "UMAP2"
    ) +
    ggplot2::theme_classic()
  path <- file.path(correlation_result_dir(), paste0("rcna_", score, ".pdf"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = 8, height = 6)
  path
}

validate_correlation_result <- function(current_files, legacy_file,
                                        tolerance = 1e-8) {
  current <- current_files[basename(current_files) == basename(legacy_file)]
  stopifnot(length(current) == 1L)
  observed <- readr::read_csv(current, show_col_types = FALSE)
  expected <- readr::read_csv(legacy_file, show_col_types = FALSE)
  stopifnot(
    identical(names(observed), names(expected)),
    identical(as.character(observed$cell_type), as.character(expected$cell_type)),
    all(vapply(setdiff(names(observed), "cell_type"), function(column) {
      isTRUE(all.equal(observed[[column]], expected[[column]],
                       check.attributes = FALSE, tolerance = tolerance))
    }, logical(1)))
  )
  tibble::tibble(artifact = basename(legacy_file), passed = TRUE)
}
