abundance_result_dir <- function() {
  file.path("results", "targets", "abundance")
}

prepare_abundance_object <- function(object, lookup) {
  required_meta <- c("cluster", "sample", "patient", "tissue", "diagnosis", "group")
  required_lookup <- c("patient", "group2", "sex", "age")
  stopifnot(
    inherits(object, "Seurat"),
    all(required_meta %in% colnames(object[[]])),
    all(required_lookup %in% names(lookup)),
    !anyDuplicated(lookup$patient),
    all(object$patient %in% lookup$patient),
    all(c("CSF", "PBMC") %in% unique(as.character(object$tissue)))
  )

  object$group2 <- lookup$group2[match(object$patient, lookup$patient)]
  stopifnot(!anyNA(object$group2))
  object
}

make_abundance_table <- function(object, row_var, col_var) {
  stopifnot(all(c(row_var, col_var) %in% colnames(object[[]])))
  absolute <- table(object[[row_var, drop = TRUE]], object[[col_var, drop = TRUE]]) |>
    as.data.frame.matrix() |>
    tibble::rownames_to_column("cell")
  percentage <- absolute |>
    dplyr::mutate(dplyr::across(
      tidyselect::where(is.numeric),
      ~ round(.x / sum(.x) * 100, 2)
    ))
  list(absolute = absolute, percentage = percentage)
}

make_abundance_tables <- function(object) {
  list(
    sample = make_abundance_table(object, "cluster", "sample"),
    tissue_diagnosis = make_abundance_table(
      object, "cluster", "tissue_diagnosis"
    ),
    tissue_group = make_abundance_table(object, "cluster", "tissue_group")
  )
}

write_abundance_table <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(table, path)
  path
}

abundance_comparisons <- function() {
  make_pairs <- function(conditions, group_column) {
    pairs <- unlist(lapply(seq.int(2L, length(conditions)), function(second) {
      lapply(seq_len(second - 1L), function(first) {
        conditions[c(first, second)]
      })
    }), recursive = FALSE)
    dplyr::bind_rows(lapply(pairs, function(pair) {
      tibble::tibble(
        tissue = c("CSF", "PBMC"),
        condition1 = pair[[1L]],
        condition2 = pair[[2L]],
        group_column = group_column
      )
    }))
  }

  dplyr::bind_rows(
    make_pairs(c("PNP", "CTRL"), "group"),
    make_pairs(c("IN", "NIN", "CTRL"), "group2"),
    make_pairs(c("CIDP", "GBS", "CIAP", "CTRL"), "diagnosis")
  ) |>
    dplyr::mutate(
      comparison = paste(tissue, condition1, "vs", condition2, sep = "_")
    )
}

run_abundance_propeller <- function(
  object,
  lookup,
  comparisons,
  seed = 123L,
  n_perms = 1000L
) {
  stopifnot(
    all(comparisons$group_column %in% colnames(object[[]])),
    all(c("patient", "sex", "age", "group", "group2", "diagnosis") %in%
      names(lookup)),
    !anyDuplicated(comparisons$comparison)
  )

  result <- lapply(seq_len(nrow(comparisons)), function(index) {
    comparison <- comparisons[index, ]
    tissue_object <- object[, object$tissue == comparison$tissue]
    comparison_lookup <- lookup |>
      dplyr::filter(.data[[comparison$group_column]] %in% c(
        comparison$condition1, comparison$condition2
      )) |>
      droplevels()
    formula <- stats::as.formula(paste0(
      "~0 + ", comparison$group_column, " + sex + age"
    ))
    result <- scMisc::propellerCalc(
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
      result,
      tissue_object,
      lookup,
      comparison$group_column,
      comparison$condition1,
      comparison$condition2,
      seed = seed + index - 1L,
      n_perms = n_perms
    )
    dplyr::mutate(
      result,
      permFDP_sig = .data$P.Value < threshold,
      permFDP_threshold = threshold
    )
  })
  names(result) <- substr(comparisons$comparison, 1L, 31L)
  result
}

deterministic_propeller_threshold <- function(
  result,
  object,
  lookup,
  group_column,
  condition1,
  condition2,
  seed,
  n_perms,
  fdr = 0.1
) {
  metadata <- object[[]] |>
    dplyr::filter(.data[[group_column]] %in% c(condition1, condition2))
  counts <- table(metadata$cluster, metadata$patient)
  proportions <- sweep(counts, 2L, colSums(counts), "/")
  proportions <- proportions[result$cluster, , drop = FALSE]
  sample_group <- lookup[[group_column]][match(colnames(proportions), lookup$patient)]
  stopifnot(!anyNA(sample_group), n_perms > 0L)

  nc <- sum(sample_group == condition1)
  nt <- sum(sample_group == condition2)
  c_into_c <- round(nc * nc / (nc + nt))
  t_into_c <- nc - c_into_c
  control_design <- c(rep.int(1L, c_into_c), rep.int(2L, t_into_c))
  test_design <- c(rep.int(1L, t_into_c), rep.int(2L, nt - t_into_c))

  permutation_p <- withr::with_seed(seed, replicate(n_perms, {
    design <- c(sample(control_design), sample(test_design))
    first <- proportions[, design == 1L, drop = FALSE]
    second <- proportions[, design == 2L, drop = FALSE]
    pooled_variance <- (
      (ncol(first) - 1L) * apply(first, 1L, stats::var) +
        (ncol(second) - 1L) * apply(second, 1L, stats::var)
    ) / (ncol(first) + ncol(second) - 2L)
    statistic <- (rowMeans(first) - rowMeans(second)) /
      sqrt(pooled_variance * (1 / ncol(first) + 1 / ncol(second)))
    2 * stats::pt(
      -abs(statistic), df = ncol(first) + ncol(second) - 2L
    )
  }))
  permutation_p <- matrix(
    permutation_p, nrow = nrow(proportions), ncol = n_perms
  )
  observed <- sort(result$P.Value)
  false_hits <- vapply(observed, function(threshold) {
    mean(colSums(permutation_p <= threshold, na.rm = TRUE))
  }, numeric(1))
  estimated_fdp <- false_hits / seq_along(observed)
  accepted <- which(estimated_fdp <= fdr)
  if (!length(accepted)) {
    return(observed[[1L]] / 2)
  }
  best <- max(accepted)
  if (best == length(observed)) {
    return(if (observed[[best]] + 0.05 <= 1) {
      observed[[best]] + 0.05
    } else {
      (observed[[best]] + 1) / 2
    })
  }
  mean(observed[c(best, best + 1L)])
}

write_propeller_results <- function(results, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(results, path)
  path
}

write_propeller_plots <- function(results, colors, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(output_dir, paste0("propeller_", names(results), ".pdf"))
  for (index in seq_along(results)) {
    data <- results[[index]]
    threshold <- unique(data$permFDP_threshold)[[1L]]
    stopifnot(is.finite(threshold), threshold > 0)
    y <- -log10(data$P.Value)
    labels <- ifelse(
      y >= -log10(threshold) & abs(data$log2ratio) >= 0.5,
      data$cluster,
      NA_character_
    )
    plot <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = log2ratio, y = -log10(P.Value), color = cluster)
    ) +
      ggplot2::geom_point(size = 3) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = labels),
        nudge_y = 0.07,
        max.overlaps = 20,
        na.rm = TRUE,
        seed = 123L + index - 1L
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(threshold), color = "blue", linetype = "dashed"
      ) +
      ggplot2::geom_vline(xintercept = 0, color = "red") +
      ggplot2::geom_vline(
        xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed"
      ) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::labs(x = expression(Log[2] ~ "fold change"),
                    y = expression(-Log[10] ~ "p value")) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(paths[[index]], plot, width = 3, height = 3)
  }
  paths
}

stacked_abundance_data <- function(object, tissue, x_axis) {
  metadata <- object[[]] |>
    dplyr::filter(.data$tissue == .env$tissue)
  result <- table(metadata$cluster, metadata[[x_axis]]) |>
    as.data.frame.matrix() |>
    tibble::rownames_to_column("cluster") |>
    tidyr::pivot_longer(-cluster, names_to = "group", values_to = "count")
  result$cluster <- factor(result$cluster, levels = object@misc$cluster_order)
  result$group <- factor(result$group, levels = unique(metadata[[x_axis]]))
  dplyr::filter(result, count != 0)
}

write_stacked_abundance_plot <- function(object, tissue, x_axis, path, width) {
  data <- stacked_abundance_data(object, tissue, x_axis)
  plot <- ggplot2::ggplot(
    data, ggplot2::aes(x = group, y = count, fill = cluster)
  ) +
    ggplot2::geom_col(
      color = "black", linewidth = 0.1, position = "fill"
    ) +
    ggplot2::scale_fill_manual(values = object@misc$cluster_col) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::guides(fill = ggplot2::guide_legend(title = NULL)) +
    ggplot2::labs(x = NULL, y = "Proportion of cells") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.3)
    )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, plot, width = width, height = 10)
  path
}

write_abundance_boxplot <- function(object, tissue, path) {
  tissue_object <- object[, object$tissue == tissue]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  generated <- file.path(
    dirname(path), "boxplot_cluster_tissue_object_diagnosis.pdf"
  )
  unlink(generated)
  suppressMessages(suppressWarnings(scMisc::abBoxPlot(
      object = tissue_object,
      cluster_idents = "cluster",
      sample = "sample",
      cluster_order = object@misc$cluster_order,
      group_by = "diagnosis",
      group_order = object@misc$diagnosis_order,
      color = object@misc$diagnosis_col,
      number_of_tests = choose(9, 2),
      width = 14,
      dir_output = dirname(path)
    )))
  stopifnot(file.exists(generated), file.rename(generated, path))
  path
}

abundance_pca_input <- function(object, tissues, diagnoses = NULL) {
  metadata <- object[[]] |>
    dplyr::filter(.data$tissue %in% .env$tissues)
  if (!is.null(diagnoses)) {
    metadata <- dplyr::filter(metadata, .data$diagnosis %in% .env$diagnoses)
  }
  stopifnot(nrow(metadata) > 0L)

  transformed <- lapply(tissues, function(tissue) {
    tissue_metadata <- dplyr::filter(metadata, .data$tissue == .env$tissue)
    props <- speckle::getTransformedProps(
      clusters = tissue_metadata$cluster,
      sample = tissue_metadata$sample,
      transform = "logit"
    )$TransformedProps
    matrix <- t(as.data.frame.matrix(props))
    sample_lookup <- dplyr::distinct(
      tissue_metadata, .data$sample, .data$patient, .data$diagnosis
    )
    patient <- sample_lookup$patient[match(rownames(matrix), sample_lookup$sample)]
    stopifnot(!anyNA(patient), !anyDuplicated(patient))
    rownames(matrix) <- patient
    colnames(matrix) <- paste0(colnames(matrix), "_", tissue)
    matrix
  })

  if (length(transformed) == 1L) {
    matrix <- transformed[[1L]]
  } else {
    common_patients <- Reduce(intersect, lapply(transformed, rownames))
    stopifnot(length(common_patients) >= 3L)
    matrix <- do.call(cbind, lapply(transformed, function(x) {
      x[common_patients, , drop = FALSE]
    }))
  }
  patient_lookup <- dplyr::distinct(metadata, .data$patient, .data$diagnosis)
  condition <- patient_lookup$diagnosis[match(rownames(matrix), patient_lookup$patient)]
  stopifnot(!anyNA(condition), all(is.finite(matrix)))
  list(matrix = matrix, condition = condition)
}

run_abundance_pca <- function(input) {
  ncp <- min(30L, nrow(input$matrix) - 1L, ncol(input$matrix))
  stopifnot(ncp >= 2L)
  result <- FactoMineR::PCA(
    input$matrix, scale.unit = FALSE, ncp = ncp, graph = FALSE
  )
  list(result = result, condition = input$condition)
}

make_abundance_pcas <- function(object) {
  main_diagnoses <- c("CIDP", "GBS", "CIAP", "CTRL")
  inputs <- list(
    csf_main = abundance_pca_input(object, "CSF", main_diagnoses),
    pbmc_main = abundance_pca_input(object, "PBMC", main_diagnoses),
    combined_main = abundance_pca_input(
      object, c("CSF", "PBMC"), main_diagnoses
    ),
    combined_all = abundance_pca_input(object, c("CSF", "PBMC"))
  )
  lapply(inputs, run_abundance_pca)
}

write_abundance_pca_tables <- function(pcas, path) {
  sheets <- list()
  for (name in names(pcas)) {
    pca <- pcas[[name]]$result
    sheets[[paste0(name, "_scores")]] <- pca$ind$coord |>
      as.data.frame() |>
      tibble::rownames_to_column("patient") |>
      dplyr::mutate(condition = as.character(pcas[[name]]$condition), .after = patient)
    sheets[[paste0(name, "_loadings")]] <- pca$var$coord |>
      as.data.frame() |>
      tibble::rownames_to_column("cluster")
    sheets[[paste0(name, "_variance")]] <- pca$eig |>
      as.data.frame() |>
      tibble::rownames_to_column("component")
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(sheets, path)
  path
}

write_abundance_pca_plots <- function(pcas, colors, output_dir, seed = 123L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(output_dir, paste0("pca_", names(pcas), ".pdf"))
  for (index in seq_along(pcas)) {
    pca <- pcas[[index]]$result
    condition <- pcas[[index]]$condition
    plot <- withr::with_seed(seed + index - 1L, {
      eigen <- factoextra::fviz_eig(
        pca, addlabels = TRUE, ylim = c(0, 50), ncp = min(7L, nrow(pca$eig))
      )
      variables <- factoextra::fviz_pca_var(
        pca,
        col.var = "contrib",
        gradient.cols = viridis::viridis(100),
        repel = TRUE,
        select.var = list(contrib = min(20L, nrow(pca$var$coord)))
      ) + ggplot2::labs(title = NULL) + ggplot2::theme_classic()
      individuals <- factoextra::fviz_pca_ind(
        pca, pointsize = 5, pointshape = 21, fill = "#E7B800",
        col.ind = "black", axes.linetype = "solid"
      ) + ggplot2::labs(title = NULL, x = "PC1", y = "PC2") +
        ggplot2::theme_classic()
      group_data <- pca$ind$coord |>
        as.data.frame() |>
        dplyr::mutate(condition = .env$condition)
      groups <- ggplot2::ggplot(
        group_data,
        ggplot2::aes(x = Dim.1, y = Dim.2, fill = condition)
      ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        ggplot2::geom_point(shape = 21, size = 5, color = "black") +
        ggplot2::scale_fill_manual(values = colors, name = "group") +
        ggplot2::labs(title = NULL, x = "PC1", y = "PC2") +
        ggplot2::theme_classic()
      patchwork::wrap_plots(eigen, variables, individuals, groups, ncol = 4)
    })
    ggplot2::ggsave(paths[[index]], plot, width = 20, height = 5)
  }
  paths
}

legacy_abundance_table_paths <- function() {
  paths <- file.path(
    "results", "abundance",
    c(
      "abundance_tbl_sc_merge_sample.xlsx",
      "abundance_tbl_sc_merge_tissue_diagnosis.xlsx",
      "abundance_tbl_sc_merge_tissue_group.xlsx"
    )
  )
  stopifnot(all(file.exists(paths)))
  paths
}

validate_abundance_tables <- function(tables, legacy_paths) {
  names(legacy_paths) <- c("sample", "tissue_diagnosis", "tissue_group")
  checks <- lapply(names(legacy_paths), function(name) {
    for (sheet in c("absolute", "percentage")) {
      legacy <- readxl::read_excel(legacy_paths[[name]], sheet = sheet)
      current <- tables[[name]][[sheet]]
      stopifnot(
        identical(names(current), names(legacy)),
        nrow(current) == nrow(legacy),
        identical(as.character(current$cell), as.character(legacy$cell)),
        isTRUE(all.equal(
          as.matrix(current[-1L]), as.matrix(legacy[-1L]),
          check.attributes = FALSE, tolerance = 1e-12
        ))
      )
    }
    tibble::tibble(artifact = basename(legacy_paths[[name]]), passed = TRUE)
  })
  dplyr::bind_rows(checks)
}

legacy_propeller_path <- function() {
  path <- file.path("results", "abundance", "propeller_results_all.xlsx")
  stopifnot(file.exists(path))
  path
}

validate_propeller_results <- function(results, legacy_path) {
  stopifnot(identical(names(results), readxl::excel_sheets(legacy_path)))
  checks <- lapply(names(results), function(sheet) {
    legacy <- readxl::read_excel(legacy_path, sheet = sheet)
    current <- results[[sheet]]
    deterministic_columns <- setdiff(
      names(current), c("permFDP_sig", "permFDP_threshold")
    )
    stopifnot(
      identical(deterministic_columns, names(legacy)[seq_along(deterministic_columns)]),
      identical(as.character(current$cluster), as.character(legacy$cluster)),
      all(vapply(deterministic_columns[-1L], function(column) {
        isTRUE(all.equal(
          as.numeric(current[[column]]), as.numeric(legacy[[column]]),
          check.attributes = FALSE, tolerance = 1e-10
        ))
      }, logical(1)))
    )
    tibble::tibble(comparison = sheet, passed = TRUE)
  })
  dplyr::bind_rows(checks)
}
