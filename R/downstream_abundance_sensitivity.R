prepare_treatment_naive_abundance <- function(object, lookup) {
  required_lookup <- c("patient", "diagnosis", "therapy", "immunomodulators")
  stopifnot(
    inherits(object, "Seurat"),
    all(required_lookup %in% names(lookup)),
    all(c("patient", "sample", "tissue", "diagnosis") %in% colnames(object[[]]))
  )

  excluded <- lookup |>
    dplyr::filter(
      as.character(.data$diagnosis) %in% c("CIDP", "GBS"),
      .data$therapy == "yes"
    ) |>
    dplyr::select(tidyselect::all_of(required_lookup))
  excluded_counts <- table(as.character(excluded$diagnosis))
  stopifnot(
    !anyDuplicated(excluded$patient),
    identical(as.integer(excluded_counts[c("CIDP", "GBS")]), c(2L, 2L)),
    all(excluded$patient %in% object$patient)
  )

  filtered_lookup <- dplyr::filter(lookup, !.data$patient %in% excluded$patient)
  filtered_object <- object[, !object$patient %in% excluded$patient]
  cohort <- filtered_object[[]] |>
    dplyr::filter(as.character(.data$diagnosis) %in% c("CIDP", "GBS", "CTRL")) |>
    dplyr::group_by(.data$diagnosis, .data$tissue) |>
    dplyr::summarise(
      patients = dplyr::n_distinct(.data$patient),
      samples = dplyr::n_distinct(.data$sample),
      cells = dplyr::n(),
      .groups = "drop"
    )

  list(
    object = filtered_object,
    lookup = filtered_lookup,
    excluded = excluded,
    cohort = cohort
  )
}

treatment_naive_abundance_comparisons <- function(comparisons) {
  result <- comparisons |>
    dplyr::filter(
      .data$group_column == "diagnosis",
      (.data$condition1 %in% c("CIDP", "GBS") & .data$condition2 == "CTRL") |
        (.data$condition1 == "CIDP" & .data$condition2 == "GBS")
    )
  stopifnot(
    nrow(result) == 6L,
    setequal(result$tissue, c("CSF", "PBMC")),
    setequal(result$condition1, c("CIDP", "GBS")),
    setequal(result$condition2, c("GBS", "CTRL"))
  )
  result
}

run_treatment_naive_abundance_propeller <- function(
  object,
  lookup,
  comparisons,
  seed = 123L,
  n_perms = 1000L
) {
  stopifnot(
    all(comparisons$group_column %in% colnames(object[[]])),
    all(c("patient", "sex", "age", "diagnosis") %in% names(lookup)),
    !anyDuplicated(comparisons$comparison)
  )

  result <- lapply(seq_len(nrow(comparisons)), function(index) {
    comparison <- comparisons[index, ]
    tissue_object <- object[, object$tissue == comparison$tissue]
    tissue_patients <- unique(tissue_object$patient)
    comparison_lookup <- lookup |>
      dplyr::filter(
        .data$patient %in% .env$tissue_patients,
        .data[[comparison$group_column]] %in% c(
          comparison$condition1, comparison$condition2
        )
      ) |>
      droplevels()
    covariates <- c(
      comparison$group_column,
      if (dplyr::n_distinct(comparison_lookup$sex) > 1L) "sex",
      if (dplyr::n_distinct(comparison_lookup$age) > 1L) "age"
    )
    formula <- stats::as.formula(paste("~0 +", paste(covariates, collapse = " + ")))
    model_formula <- paste(deparse(formula), collapse = "")
    comparison_result <- scMisc::propellerCalc(
      seu_obj1 = tissue_object,
      condition1 = comparison$condition1,
      condition2 = comparison$condition2,
      cluster_col = "cluster",
      meta_col = comparison$group_column,
      lookup = comparison_lookup,
      sample_col = "patient",
      formula = formula,
      min_cells = 9,
      adjustment_method = "BH"
    )
    threshold <- deterministic_propeller_threshold(
      comparison_result,
      tissue_object,
      comparison_lookup,
      comparison$group_column,
      comparison$condition1,
      comparison$condition2,
      seed = seed + index - 1L,
      n_perms = n_perms
    )
    dplyr::mutate(
      comparison_result,
      model_formula = .env$model_formula,
      permFDP_sig = .data$P.Value < threshold,
      permFDP_threshold = threshold
    )
  })
  names(result) <- substr(comparisons$comparison, 1L, 31L)
  result
}

compare_propeller_sensitivity <- function(primary, sensitivity) {
  stopifnot(all(names(sensitivity) %in% names(primary)))
  result <- lapply(names(sensitivity), function(comparison) {
    full <- primary[[comparison]] |>
      dplyr::select(
        "cluster",
        full_log2ratio = "log2ratio",
        full_p_value = "P.Value",
        full_fdr = "FDR",
        full_perm_fdp_sig = "permFDP_sig"
      )
    untreated <- sensitivity[[comparison]] |>
      dplyr::select(
        "cluster",
        untreated_log2ratio = "log2ratio",
        untreated_p_value = "P.Value",
        untreated_fdr = "FDR",
        untreated_perm_fdp_sig = "permFDP_sig"
      )
    dplyr::full_join(full, untreated, by = "cluster") |>
      dplyr::mutate(
        log2ratio_change = .data$untreated_log2ratio - .data$full_log2ratio,
        direction_concordant = sign(.data$full_log2ratio) ==
          sign(.data$untreated_log2ratio)
      )
  })
  names(result) <- names(sensitivity)
  result
}
