flow_result_dir <- function() {
  file.path("results", "targets", "flow")
}

flow_diagnosis_order <- function() {
  c("CTRL", "CIAP", "CIDP", "GBS", "MAG", "MFS", "PNC", "CAN", "PPN")
}

flow_diagnosis_colors <- function() {
  stats::setNames(
    pals::cols25(length(flow_diagnosis_order())),
    flow_diagnosis_order()
  )
}

read_flow_seed_lookup <- function(path) {
  lookup <- readxl::read_excel(path) |>
    janitor::clean_names() |>
    dplyr::select(patient, pseudonym, cohort, sex, age, diagnosis)
  stopifnot(
    !anyDuplicated(lookup$patient),
    !anyNA(lookup$patient),
    all(c("sex", "age", "diagnosis") %in% names(lookup))
  )
  lookup
}

read_flow_frontiers_lookup <- function(path) {
  readxl::read_excel(path) |>
    dplyr::select(
      name, id, date...6, birth, dx, age, sex, protein, disruption:lac
    ) |>
    dplyr::mutate(
      measure_date = suppressWarnings(lubridate::dmy(date...6)),
      birth_date = suppressWarnings(lubridate::dmy(birth))
    ) |>
    dplyr::filter(dx %in% c("CIDP", "IIH", "GBS")) |>
    dplyr::mutate(
      dx = dplyr::if_else(dx == "IIH", "CTRL", dx),
      sex = dplyr::if_else(sex == "M", "male", "female")
    ) |>
    tidyr::separate_wider_delim(
      name, delim = ", ", names = c("last_name", "first_name")
    ) |>
    dplyr::mutate(
      first_name_letter = substring(first_name, 1L, 1L),
      last_name_letter = substring(last_name, 1L, 1L),
      id = as.numeric(id)
    )
}

prepare_flow_frontiers <- function(flow_raw, frontiers_lookup) {
  flow_frontiers <- flow_raw |>
    dplyr::filter(cohort == "frontiers") |>
    dplyr::select(-measure_date) |>
    tidyr::separate_wider_regex(
      flow_file_name,
      patterns = c(
        last_name_letter = "^[A-Za-z]",
        first_name_letter = "[A-Za-z]",
        birth_date = "\\d{8}",
        extra = ".*",
        measure_date = "20\\d{2}-\\d{2}-\\d{2}",
        extra2 = ".*"
      ),
      cols_remove = FALSE
    ) |>
    dplyr::mutate(
      birth_date1 = suppressWarnings(lubridate::dmy(birth_date)),
      birth_date2 = suppressWarnings(lubridate::mdy(birth_date)),
      birth_date = dplyr::coalesce(birth_date1, birth_date2),
      measure_date = lubridate::ymd(measure_date)
    ) |>
    dplyr::select(-extra, -extra2, -birth_date1, -birth_date2)

  matched_1 <- flow_frontiers |>
    dplyr::inner_join(
      frontiers_lookup,
      dplyr::join_by(
        first_name_letter, last_name_letter, birth_date, measure_date
      ),
      suffix = c("_flow", "_lookup")
    ) |>
    dplyr::mutate(match_level = "all_four")
  unmatched_1 <- flow_frontiers |>
    dplyr::anti_join(
      frontiers_lookup,
      dplyr::join_by(
        first_name_letter, last_name_letter, birth_date, measure_date
      )
    )

  matched_2 <- unmatched_1 |>
    dplyr::inner_join(
      frontiers_lookup,
      dplyr::join_by(first_name_letter, birth_date, measure_date),
      suffix = c("_flow", "_lookup")
    ) |>
    dplyr::mutate(match_level = "three_no_last")
  unmatched_2 <- unmatched_1 |>
    dplyr::anti_join(
      frontiers_lookup,
      dplyr::join_by(first_name_letter, birth_date, measure_date)
    )

  matched_3 <- unmatched_2 |>
    dplyr::inner_join(
      dplyr::select(frontiers_lookup, -measure_date, -last_name_letter),
      dplyr::join_by(first_name_letter, birth_date),
      suffix = c("_flow", "_lookup")
    ) |>
    dplyr::mutate(match_level = "two_name_birth")
  unmatched_3 <- unmatched_2 |>
    dplyr::anti_join(
      frontiers_lookup,
      dplyr::join_by(first_name_letter, birth_date)
    )

  matched_4 <- unmatched_3 |>
    dplyr::inner_join(
      dplyr::select(frontiers_lookup, -birth_date, -last_name_letter),
      dplyr::join_by(first_name_letter, measure_date),
      suffix = c("_flow", "_lookup")
    ) |>
    dplyr::mutate(match_level = "two_name_date")
  unmatched_4 <- unmatched_3 |>
    dplyr::anti_join(
      frontiers_lookup,
      dplyr::join_by(first_name_letter, measure_date)
    )

  matched <- dplyr::bind_rows(matched_1, matched_2, matched_3, matched_4) |>
    dplyr::rename(diagnosis = dx)
  summary <- tibble::tibble(
    match_level = c(
      "all_four", "three_no_last", "two_name_birth", "two_name_date",
      "unmatched"
    ),
    records = c(
      nrow(matched_1), nrow(matched_2), nrow(matched_3), nrow(matched_4),
      nrow(unmatched_4)
    )
  )
  stopifnot(nrow(matched) > 0L, sum(summary$records) == nrow(flow_frontiers))
  list(matched = matched, summary = summary)
}

prepare_flow_data <- function(flow_file, frontiers_file, seed_lookup_file) {
  flow_raw <- readxl::read_excel(flow_file)
  frontiers_raw <- readxl::read_excel(frontiers_file)
  frontiers_lookup <- read_flow_frontiers_lookup(frontiers_file)
  frontiers <- prepare_flow_frontiers(flow_raw, frontiers_lookup)
  seed_lookup <- read_flow_seed_lookup(seed_lookup_file)

  flow_seed <- flow_raw |>
    dplyr::filter(!is.na(patient)) |>
    dplyr::group_by(patient, tissue) |>
    dplyr::filter(measure_date == min(measure_date)) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::inner_join(seed_lookup, dplyr::join_by(patient))

  flow_final <- dplyr::bind_rows(flow_seed, frontiers$matched) |>
    dplyr::mutate(
      diagnosis = factor(diagnosis, levels = flow_diagnosis_order()),
      id = as.character(id),
      patient_id = dplyr::coalesce(patient, id)
    ) |>
    dplyr::distinct(patient_id, tissue, .keep_all = TRUE)
  stopifnot(
    !anyNA(flow_final$diagnosis),
    !anyNA(flow_final$age),
    !anyNA(flow_final$sex),
    !anyNA(flow_final$tissue),
    !anyDuplicated(flow_final[c("patient_id", "tissue")])
  )
  flow <- split(flow_final, flow_final$tissue)
  stopifnot(identical(names(flow), c("blood", "CSF")))

  list(
    flow = flow,
    gating_comparison = make_flow_gating_comparison(
      flow_raw, frontiers_raw, frontiers$matched
    ),
    match_summary = frontiers$summary
  )
}

make_flow_gating_comparison <- function(flow_raw, frontiers_raw, flow_matched) {
  flow_louisa <- flow_matched |>
    dplyr::filter(tissue == "CSF") |>
    dplyr::inner_join(dplyr::select(flow_raw, flow_file_name, Gran:intMono))
  flow_andi <- frontiers_raw |>
    dplyr::mutate(
      measure_date = suppressWarnings(lubridate::dmy(date...6)),
      birth_date = suppressWarnings(lubridate::dmy(birth))
    ) |>
    dplyr::filter(dx %in% c("CIDP", "IIH", "GBS")) |>
    dplyr::mutate(
      dx = dplyr::if_else(dx == "IIH", "CTRL", dx),
      sex = dplyr::if_else(sex == "M", "male", "female")
    ) |>
    tidyr::separate_wider_delim(
      name, delim = ", ", names = c("last_name", "first_name")
    ) |>
    dplyr::mutate(
      first_name_letter = substring(first_name, 1L, 1L),
      last_name_letter = substring(last_name, 1L, 1L)
    ) |>
    dplyr::select(
      measure_date, birth_date, first_name_letter, last_name_letter,
      id, dx, age, sex, protein, disruption:cd4cd8ratio
    ) |>
    dplyr::mutate(tissue = "CSF") |>
    dplyr::select(
      -dplyr::contains("_abs"), -cd45cells_pct, -cd4cd8cells_pct,
      -tcellshladr_pct, -nkcellshladr_pct, -cd4cd8ratio
    ) |>
    dplyr::rename(
      patient = id, diagnosis = dx, Lymph = lymphos_pct, Mono = monos_pct,
      T = tcells_pct, CD4 = cd4cells_pct, CD8 = cd8cells_pct,
      B = bcells_pct, NK = nkcells_pct, NKT = nktcells_pct,
      actCD4 = cd4cellshladr_pct, actCD8 = cd8cellshladr_pct,
      Plasma = plasmacells_pct, cMono = monosclassical_pct,
      ncMono = monosatypical_pct, brightNK = nkcellsbright_pct,
      dimNK = nkcellsdim_pct
    )

  flow_vars <- names(dplyr::select(flow_louisa, Gran:intMono))
  flow_vars <- intersect(flow_vars, names(flow_andi))
  comparison <- flow_louisa |>
    dplyr::select(
      flow_file_name, diagnosis, first_name_letter, last_name_letter,
      birth_date, measure_date, dplyr::all_of(flow_vars)
    ) |>
    dplyr::inner_join(
      dplyr::select(
        flow_andi, first_name_letter, last_name_letter, birth_date,
        measure_date, dplyr::all_of(flow_vars)
      ),
      by = c(
        "first_name_letter", "last_name_letter", "birth_date", "measure_date"
      ),
      suffix = c("_louisa", "_andi")
    )
  for (variable in flow_vars) {
    comparison[[paste0(variable, "_diff")]] <-
      comparison[[paste0(variable, "_louisa")]] -
      comparison[[paste0(variable, "_andi")]]
    comparison[[paste0(variable, "_rel_diff")]] <-
      comparison[[paste0(variable, "_diff")]] /
      comparison[[paste0(variable, "_andi")]]
  }
  id_columns <- c(
    "flow_file_name", "diagnosis", "first_name_letter", "last_name_letter",
    "birth_date", "measure_date"
  )
  value_columns <- unlist(lapply(flow_vars, function(variable) {
    paste0(
      variable,
      c("_louisa", "_andi", "_diff", "_rel_diff")
    )
  }), use.names = FALSE)
  dplyr::select(comparison, dplyr::all_of(c(id_columns, value_columns)))
}

flow_variables <- function(flow) {
  variables <- names(dplyr::select(flow$CSF, Lymph:intMono))
  stopifnot(length(variables) == 17L, all(variables %in% names(flow$blood)))
  variables
}

flow_comparisons <- function() {
  tidyr::crossing(
    tissue = c("CSF", "blood"),
    comparison = c("CIDP_GBS", "CIDP_CTRL", "GBS_CTRL")
  ) |>
    tidyr::separate_wider_delim(
      comparison, delim = "_", names = c("group1", "group2"),
      cols_remove = FALSE
    ) |>
    dplyr::mutate(name = paste(comparison, tissue, sep = "_"))
}

flow_fold_changes <- function(data, variables, group1, group2) {
  data |>
    dplyr::group_by(diagnosis) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(variables), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(dplyr::all_of(variables), names_to = "variable") |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = value) |>
    dplyr::mutate(log2_ratio = log2(.data[[group1]] / .data[[group2]]))
}

flow_regression_statistics <- function(
  data, variables, seed, fdr_threshold = 0.1, n_permutations = 1000L
) {
  statistics <- dplyr::bind_rows(lapply(variables, function(variable) {
    model <- stats::lm(
      stats::as.formula(paste0(variable, "~diagnosis + sex + age")),
      data = data
    )
    tibble::tibble(
      variable = variable,
      p.value = dplyr::pull(dplyr::slice(broom::tidy(model), 2L), p.value)
    )
  }))
  quantities <- t(as.data.frame(dplyr::select(data, dplyr::all_of(variables))))
  complete <- stats::complete.cases(quantities)
  quantities <- quantities[complete, , drop = FALSE]
  statistics <- statistics[complete, , drop = FALSE]
  stopifnot(nrow(quantities) > 0L, n_permutations > 0L)
  levels_seen <- unique(data$diagnosis)
  design <- ifelse(data$diagnosis == levels_seen[[1L]], 1L, 2L)
  threshold <- withr::with_seed(seed, {
    permFDP::permFDP.adjust.threshold(
      pVals = statistics$p.value,
      threshold = fdr_threshold,
      myDesign = design,
      intOnly = quantities,
      nPerms = n_permutations
    )
  })
  dplyr::mutate(
    statistics,
    p.adj.threshold = threshold,
    significant = p.value < threshold,
    neg_log10_p = -log10(p.value)
  )
}

run_flow_volcanoes <- function(
  flow, comparisons, seed = 123L, n_permutations = 1000L
) {
  variables <- flow_variables(flow)
  results <- lapply(seq_len(nrow(comparisons)), function(index) {
    comparison <- comparisons[index, ]
    data <- flow[[comparison$tissue]] |>
      dplyr::filter(diagnosis %in% c(comparison$group1, comparison$group2)) |>
      droplevels()
    fold_changes <- flow_fold_changes(
      data, variables, comparison$group1, comparison$group2
    )
    statistics <- flow_regression_statistics(
      data,
      fold_changes$variable,
      seed = seed + index - 1L,
      n_permutations = n_permutations
    )
    dplyr::left_join(fold_changes, statistics, by = "variable")
  })
  names(results) <- comparisons$name
  results
}

write_flow_object <- function(flow, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  qs::qsave(flow, path)
  path
}

write_flow_table <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(table, path)
  path
}

write_flow_volcano_tables <- function(results, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(results, path)
  path
}

flow_boxplot <- function(data, variable, colors, seed) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = diagnosis, y = .data[[variable]], fill = diagnosis)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(
      position = ggplot2::position_jitter(width = 0.2, height = 0, seed = seed),
      alpha = 0.5,
      size = 0.5
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = NULL, y = NULL, title = variable) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        size = 12, angle = 90, vjust = 0.5, hjust = 1
      ),
      plot.title = ggplot2::element_text(size = 20)
    )
}

write_flow_boxplots <- function(flow, output_dir, seed = 123L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  variables <- flow_variables(flow)
  colors <- flow_diagnosis_colors()
  paths <- file.path(output_dir, c("csf_con_plots_dx.pdf", "blood_con_plots_dx.pdf"))
  tissues <- c("CSF", "blood")
  for (index in seq_along(tissues)) {
    plots <- lapply(seq_along(variables), function(variable_index) {
      flow_boxplot(
        flow[[tissues[[index]]]], variables[[variable_index]], colors,
        seed + 100L * index + variable_index
      )
    })
    patch <- patchwork::wrap_plots(plots, ncol = 4L)
    ggplot2::ggsave(paths[[index]], patch, width = 10, height = 15)
  }
  paths
}

write_flow_volcano_plots <- function(results, output_dir, seed = 123L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  variables <- unique(unlist(lapply(results, `[[`, "variable"), use.names = FALSE))
  colors <- stats::setNames(pals::cols25(length(variables)), variables)
  paths <- file.path(output_dir, paste0("volcano_", names(results), ".pdf"))
  for (index in seq_along(results)) {
    data <- results[[index]]
    threshold <- unique(data$p.adj.threshold)
    stopifnot(length(threshold) == 1L, is.finite(threshold), threshold > 0)
    labels <- ifelse(
      data$neg_log10_p >= -log10(threshold) & abs(data$log2_ratio) >= 0.5,
      data$variable,
      NA_character_
    )
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = log2_ratio, y = neg_log10_p, color = variable)
    ) +
      ggplot2::geom_point(size = 3) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = labels), na.rm = TRUE, seed = seed + index - 1L
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(threshold), color = "blue", linetype = "dashed"
      ) +
      ggplot2::geom_vline(xintercept = 0, color = "red") +
      ggplot2::geom_vline(
        xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed"
      ) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(
        x = expression(Log[2] ~ "fold change"),
        y = expression(-Log[10] ~ "p value")
      ) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(paths[[index]], plot, width = 3, height = 3)
  }
  paths
}

write_flow_batch_plots <- function(flow, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  variables <- flow_variables(flow)
  paths <- file.path(output_dir, paste0(variables, "_csf.pdf"))
  for (index in seq_along(variables)) {
    variable <- variables[[index]]
    data <- flow$CSF |>
      dplyr::mutate(
        batch = ifelse(grepl("P\\d{2}", patient), "scRNAseq", "Frontiers")
      ) |>
      dplyr::filter(diagnosis %in% c("CTRL", "CIDP", "GBS"))
    plot <- ggplot2::ggplot(
      data, ggplot2::aes(x = batch, y = .data[[variable]])
    ) +
      ggplot2::geom_boxplot() +
      ggplot2::facet_wrap(~diagnosis) +
      ggplot2::labs(y = variable, x = "Batch")
    ggplot2::ggsave(paths[[index]], plot, width = 5, height = 3)
  }
  paths
}

validate_flow_object <- function(current, legacy_path, tolerance = 1e-10) {
  legacy <- qs::qread(legacy_path)
  stopifnot(
    identical(names(current), names(legacy)),
    identical(vapply(current, names, character(ncol(current[[1L]]))),
              vapply(legacy, names, character(ncol(legacy[[1L]]))))
  )
  maximum_difference <- 0
  for (tissue in names(current)) {
    stopifnot(identical(dim(current[[tissue]]), dim(legacy[[tissue]])))
    for (column in names(current[[tissue]])) {
      observed <- current[[tissue]][[column]]
      expected <- legacy[[tissue]][[column]]
      stopifnot(identical(is.na(observed), is.na(expected)))
      if (is.numeric(observed) && is.numeric(expected)) {
        difference <- abs(observed - expected)
        difference <- difference[!is.na(difference)]
        if (length(difference)) {
          maximum_difference <- max(maximum_difference, difference)
        }
      } else {
        stopifnot(identical(as.character(observed), as.character(expected)))
        if (is.factor(observed) || is.factor(expected)) {
          stopifnot(identical(levels(observed), levels(expected)))
        }
      }
    }
  }
  stopifnot(maximum_difference <= tolerance)
  tibble::tibble(
    artifact = basename(legacy_path),
    check = "exact_non_numeric_tolerant_numeric",
    max_abs_difference = maximum_difference,
    tolerance = tolerance,
    passed = TRUE
  )
}

validate_flow_gating_workbook <- function(
  current_path, legacy_path, tolerance = 1e-10
) {
  current <- readxl::read_excel(current_path)
  legacy <- readxl::read_excel(legacy_path)
  stopifnot(identical(names(current), names(legacy)), identical(dim(current), dim(legacy)))
  identifier_columns <- names(current)[seq_len(6L)]
  stopifnot(all(vapply(identifier_columns, function(column) {
    identical(as.character(current[[column]]), as.character(legacy[[column]]))
  }, logical(1))))

  # Nine Louisa-gating channels in the workbook predate flowbasic_v7. The
  # reference channels and unchanged Louisa channels remain valid baselines.
  stable_louisa <- c("Lymph", "B", "Plasma", "Mono", "cMono", "ncMono")
  compared_columns <- c(
    grep("_andi$", names(current), value = TRUE),
    unlist(lapply(stable_louisa, function(variable) {
      paste0(variable, c("_louisa", "_diff", "_rel_diff"))
    }), use.names = FALSE)
  )
  differences <- unlist(lapply(compared_columns, function(column) {
    observed <- suppressWarnings(as.numeric(current[[column]]))
    expected <- suppressWarnings(as.numeric(legacy[[column]]))
    stopifnot(
      identical(is.na(observed), is.na(expected)),
      identical(is.infinite(observed), is.infinite(expected)),
      identical(sign(observed[is.infinite(observed)]),
                sign(expected[is.infinite(expected)]))
    )
    difference <- abs(observed - expected)
    difference[is.finite(difference)]
  }), use.names = FALSE)
  maximum_difference <- if (length(differences)) max(differences) else 0
  stopifnot(is.finite(maximum_difference), maximum_difference <= tolerance)
  compared_count <- length(identifier_columns) + length(compared_columns)
  excluded_count <- ncol(current) - compared_count
  tibble::tibble(
    artifact = basename(legacy_path),
    check = "exact_identifiers_tolerant_stable_numeric",
    compared_columns = compared_count,
    excluded_stale_columns = excluded_count,
    max_abs_difference = maximum_difference,
    tolerance = tolerance,
    passed = TRUE
  )
}
