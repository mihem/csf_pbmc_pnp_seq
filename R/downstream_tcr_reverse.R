reverse_phenotype_result_dir <- function() {
  file.path(tcr_result_dir(), "reverse_phenotype")
}

reverse_phenotype_gene_panel <- function() {
  c(
    "GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "IFNG",
    "TOX", "TIGIT", "ENTPD1", "PDCD1", "LAG3", "HAVCR2",
    "KLRG1", "CX3CR1", "TBX21",
    "CXCR6", "ITGA1", "ITGAE", "ZNF683", "CD69",
    "CCR4", "CCR5", "SELL", "TCF7", "IL7R"
  )
}

#' Select expanded clones enriched in the target diagnoses relative to the
#' non-inflammatory background.
#'
#' Unit-of-analysis caveat, stated here because it must also be stated in the
#' reply: a clone is almost always confined to a single patient, so `p_value`
#' ranks clones by how atypical their target-pool expansion is, and is not an
#' inferential claim about the diagnosis group. `n_patients` is returned so any
#' clone resting on one patient is visible. Treat the FDR as a selection
#' threshold, not as evidence of a group difference.
#'
#' @param sc_tcr Seurat object carrying scRepertoire metadata.
#' @param targets Diagnoses forming the disease pool.
#' @param background Diagnoses forming the non-inflammatory background pool.
#' @param min_cells Minimum cells in the target pool for a clone to be tested.
#' @param fdr Benjamini-Hochberg threshold marking the `selected` column.
#' @return Tibble of every tested clone, ranked, with a `selected` flag.
select_disease_enriched_clones <- function(sc_tcr,
                                           targets = c("GBS", "CIDP"),
                                           background = c("CTRL", "CIAP"),
                                           min_cells = 3L,
                                           fdr = 0.1) {
  required <- c("CTaa", "patient", "diagnosis", "tissue", "cluster")
  stopifnot(
    inherits(sc_tcr, "Seurat"),
    all(required %in% colnames(sc_tcr[[]])),
    length(intersect(targets, background)) == 0L,
    min_cells >= 1L
  )

  metadata <- sc_tcr[[]] |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::filter(!is.na(.data$CTaa)) |>
    dplyr::mutate(diagnosis = as.character(.data$diagnosis))

  stopifnot(
    any(metadata$diagnosis %in% targets),
    any(metadata$diagnosis %in% background)
  )

  # Repertoire depth of each pool sets what an unremarkable clone looks like.
  n_target_pool <- sum(metadata$diagnosis %in% targets)
  n_background_pool <- sum(metadata$diagnosis %in% background)

  clones <- metadata |>
    dplyr::group_by(.data$CTaa) |>
    dplyr::summarise(
      n_target = sum(.data$diagnosis %in% targets),
      n_background = sum(.data$diagnosis %in% background),
      n_cells = dplyr::n(),
      n_patients = dplyr::n_distinct(.data$patient[.data$diagnosis %in% targets]),
      n_background_patients = dplyr::n_distinct(
        .data$patient[.data$diagnosis %in% background]
      ),
      patients = paste(
        sort(unique(.data$patient[.data$diagnosis %in% targets])),
        collapse = ","
      ),
      diagnoses = paste(sort(unique(.data$diagnosis)), collapse = ","),
      tissues = paste(sort(unique(as.character(.data$tissue))), collapse = ","),
      top_cluster = names(sort(table(as.character(.data$cluster)), decreasing = TRUE))[1L],
      .groups = "drop"
    ) |>
    dplyr::filter(.data$n_target >= min_cells)

  stopifnot(nrow(clones) > 0L)

  tested <- clones |>
    dplyr::mutate(
      p_value = purrr::map2_dbl(
        .data$n_target, .data$n_background,
        function(in_target, in_background) {
          stats::fisher.test(
            matrix(
              c(
                in_target, in_background,
                n_target_pool - in_target, n_background_pool - in_background
              ),
              nrow = 2L
            ),
            alternative = "greater"
          )$p.value
        }
      ),
      odds_ratio = (.data$n_target / pmax(n_target_pool - .data$n_target, 1L)) /
        pmax(
          .data$n_background / pmax(n_background_pool - .data$n_background, 1L),
          .Machine$double.eps
        ),
      target_pool_size = n_target_pool,
      background_pool_size = n_background_pool
    ) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p_value, "BH")) |>
    dplyr::mutate(selected = .data$p_adj < fdr) |>
    dplyr::arrange(.data$p_adj, dplyr::desc(.data$n_target))

  stopifnot(!anyNA(tested$p_value))
  tested
}

#' Select clones belonging to GBS/CIDP-enriched GLIPH2 specificity groups.
#'
#' The per-clone test in `select_disease_enriched_clones()` cannot answer
#' Reviewer 4's comment 4 on this cohort: clonotypes are almost entirely
#' patient-private, so the background count is zero for nearly every clone and
#' the ranking collapses onto clone size. Every clone it selects rests on a
#' single patient, which is the criticism the reviewer made in the first place.
#'
#' Moving the unit up to the GLIPH2 specificity group fixes that, because a
#' group collects convergent CDR3b sequences from different people. Selection
#' then rests on a quantity that several patients contribute to, and
#' `min_patients` makes the requirement explicit.
#'
#' @param sc_tcr Seurat object carrying scRepertoire metadata.
#' @param cluster_membership CDR3b-to-cluster map; `prepare_tcr_comparison()`
#'   returns this as `cluster_membership`.
#' @param diagnosis_enrichment Per-cluster, per-diagnosis enrichment table from
#'   the same call.
#' @param targets Diagnoses whose enriched groups are selected.
#' @param fdr Threshold on the enrichment FDR.
#' @param min_patients Minimum distinct study patients contributing cells to a
#'   group for it to be kept.
#' @return Tibble of one row per clonotype, carrying `CTaa` so the downstream
#'   reverse-phenotype functions take it unchanged, plus the group it came from.
select_gliph_enriched_clones <- function(sc_tcr,
                                         cluster_membership,
                                         diagnosis_enrichment,
                                         targets = c("GBS", "CIDP"),
                                         fdr = 0.1,
                                         min_patients = 2L) {
  stopifnot(
    inherits(sc_tcr, "Seurat"),
    is.data.frame(cluster_membership),
    all(c("cluster", "CDR3b") %in% colnames(cluster_membership)),
    is.data.frame(diagnosis_enrichment),
    all(c("cluster", "diagnosis") %in% colnames(diagnosis_enrichment))
  )

  # The pipeline names these p_adj/odds_ratio; the cached standalone run used
  # p.adj. Accept either rather than forcing a particular provenance.
  enrichment <- diagnosis_enrichment
  if (!"p_adj" %in% names(enrichment) && "p.adj" %in% names(enrichment)) {
    enrichment$p_adj <- enrichment$p.adj
  }
  stopifnot(all(c("p_adj", "odds_ratio") %in% names(enrichment)))

  enriched_clusters <- enrichment |>
    dplyr::filter(
      .data$diagnosis %in% targets,
      .data$p_adj < fdr,
      .data$odds_ratio > 1
    ) |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      enriched_for = paste(sort(unique(.data$diagnosis)), collapse = ","),
      min_p_adj = min(.data$p_adj),
      max_odds_ratio = max(.data$odds_ratio),
      .groups = "drop"
    )
  stopifnot(nrow(enriched_clusters) > 0L)

  metadata <- sc_tcr[[]] |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::filter(!is.na(.data$CTaa)) |>
    dplyr::mutate(
      CDR3b = sub(".*_", "", .data$CTaa),
      diagnosis = as.character(.data$diagnosis)
    )

  # `cluster` means cell type in the Seurat metadata and specificity group in
  # the GLIPH tables; rename before joining so neither becomes cluster.x/.y.
  members <- cluster_membership |>
    dplyr::distinct(.data$cluster, .data$CDR3b) |>
    dplyr::semi_join(enriched_clusters, by = "cluster") |>
    dplyr::rename(gliph_cluster = "cluster")

  # A CDR3b can belong to more than one specificity group; collapse to one row
  # per clonotype so no clone is counted twice downstream.
  clones <- metadata |>
    dplyr::inner_join(members, by = "CDR3b", relationship = "many-to-many") |>
    dplyr::group_by(.data$gliph_cluster) |>
    dplyr::mutate(group_patients = dplyr::n_distinct(.data$patient)) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$group_patients >= min_patients) |>
    dplyr::group_by(.data$CTaa) |>
    dplyr::summarise(
      CDR3b = dplyr::first(.data$CDR3b),
      gliph_clusters = paste(sort(unique(.data$gliph_cluster)), collapse = ","),
      group_patients = max(.data$group_patients),
      n_cells = dplyr::n_distinct(.data$cell_id),
      n_patients = dplyr::n_distinct(.data$patient),
      patients = paste(sort(unique(.data$patient)), collapse = ","),
      diagnoses = paste(sort(unique(.data$diagnosis)), collapse = ","),
      top_cluster = names(sort(table(as.character(.data$cluster)), decreasing = TRUE))[1L],
      .groups = "drop"
    ) |>
    dplyr::mutate(selected = TRUE) |>
    dplyr::arrange(dplyr::desc(.data$n_cells))

  stopifnot(nrow(clones) > 0L)
  clones
}

#' Every clonotype in the object, flagged by whether it was selected.
#'
#' `match_clones_by_size()` draws its comparison set from the unselected rows of
#' whatever it is given. Handing it the selected clones alone leaves nothing to
#' match against, because every row there is selected by construction, and the
#' composition then has no comparison arm at all. This builds the full universe
#' so the matching has a pool to draw from.
reverse_phenotype_clone_universe <- function(sc_tcr, clone_set) {
  stopifnot(
    inherits(sc_tcr, "Seurat"),
    is.data.frame(clone_set),
    "CTaa" %in% colnames(clone_set)
  )
  sc_tcr[[]] |>
    dplyr::filter(!is.na(.data$CTaa)) |>
    dplyr::count(.data$CTaa, name = "n_cells") |>
    dplyr::mutate(
      selected = .data$CTaa %in% clone_set$CTaa,
      n_target = .data$n_cells
    )
}

#' Size-matched comparison set for the enriched clones.
#'
#' Panel B asks whether disease-enriched clones sit in different clusters from
#' clones of the same size that are not enriched. Matching on cell count keeps
#' the comparison from being a restatement of "big clones are effector cells".
match_clones_by_size <- function(clone_set, seed = 42L) {
  stopifnot(
    is.data.frame(clone_set),
    all(c("CTaa", "n_target", "selected") %in% colnames(clone_set))
  )
  enriched <- dplyr::filter(clone_set, .data$selected)
  pool <- dplyr::filter(clone_set, !.data$selected)
  stopifnot(nrow(enriched) > 0L)

  if (nrow(pool) == 0L) {
    return(dplyr::mutate(enriched, clone_group = "disease_enriched"))
  }

  matched <- withr::with_seed(seed, {
    purrr::map_dfr(enriched$n_target, function(size) {
      distance <- abs(pool$n_target - size)
      candidates <- which(distance == min(distance))
      pool[sample(candidates, 1L), , drop = FALSE]
    })
  })

  dplyr::bind_rows(
    dplyr::mutate(enriched, clone_group = "disease_enriched"),
    dplyr::mutate(dplyr::distinct(matched), clone_group = "size_matched")
  )
}

#' Cluster composition of enriched versus size-matched clones.
reverse_phenotype_cluster_composition <- function(sc_tcr, matched_clones) {
  stopifnot(
    inherits(sc_tcr, "Seurat"),
    all(c("CTaa", "clone_group") %in% colnames(matched_clones))
  )
  lookup <- dplyr::distinct(matched_clones, .data$CTaa, .data$clone_group)

  sc_tcr[[]] |>
    dplyr::filter(!is.na(.data$CTaa)) |>
    dplyr::inner_join(lookup, by = "CTaa", relationship = "many-to-many") |>
    dplyr::count(.data$clone_group, .data$cluster, name = "n_cells") |>
    dplyr::group_by(.data$clone_group) |>
    dplyr::mutate(fraction = .data$n_cells / sum(.data$n_cells)) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$clone_group, dplyr::desc(.data$fraction))
}

#' Per-cluster enrichment estimated over repeated size-matched draws.
#'
#' A single size-matched comparison set is one arbitrary draw. Small clusters
#' can land at zero matched cells purely by that draw, which sends the exact
#' odds ratio to infinity and makes the estimate look undefined when it is only
#' unstable. CD8TEM_3 is exactly that case: 13 cells among disease-enriched
#' clones and 0 in one particular matched draw.
#'
#' Repeating the match many times and reporting the distribution of the odds
#' ratio across draws removes the artefact. The point estimate becomes the
#' median across draws and the interval reports how much the matching itself
#' moves it, which is the uncertainty a reader of the figure actually cares
#' about.
#'
#' @param sc_tcr Seurat object carrying scRepertoire metadata.
#' @param clone_set Selected clones, from `select_gliph_enriched_clones()`.
#' @param all_clones One row per clonotype with `n_target` and `selected`.
#' @param n_resamples Number of independent size-matched draws.
#' @param seed Base seed; each draw uses `seed + index - 1`.
#' @return Tibble of one row per cluster with the median odds ratio, the 2.5
#'   and 97.5 percentiles across draws, and the median Fisher p-value.
reverse_phenotype_cluster_enrichment_resampled <- function(
    sc_tcr,
    all_clones,
    n_resamples = 200L,
    seed = 42L,
    pool = list(CD8TEM_pooled = c("CD8TEM_1", "CD8TEM_2", "CD8TEM_3"))) {
  stopifnot(
    inherits(sc_tcr, "Seurat"),
    is.data.frame(all_clones),
    all(c("CTaa", "n_target", "selected") %in% colnames(all_clones)),
    n_resamples >= 1L
  )

  draws <- purrr::map_dfr(seq_len(n_resamples), function(index) {
    matched <- match_clones_by_size(all_clones, seed = seed + index - 1L)
    composition <- reverse_phenotype_cluster_composition(sc_tcr, matched)
    reverse_phenotype_cluster_enrichment(composition, pool = pool) |>
      dplyr::mutate(draw = index)
  })

  draws |>
    dplyr::group_by(.data$cluster, .data$pooled) |>
    dplyr::summarise(
      n_draws = dplyr::n(),
      n_enriched = stats::median(.data$n_enriched),
      n_matched = stats::median(.data$n_matched),
      draws_with_zero_matched = sum(.data$n_matched == 0L),
      odds_ratio = stats::median(.data$odds_ratio_corrected),
      ci_low = stats::quantile(
        .data$odds_ratio_corrected, 0.025, names = FALSE, na.rm = TRUE
      ),
      ci_high = stats::quantile(
        .data$odds_ratio_corrected, 0.975, names = FALSE, na.rm = TRUE
      ),
      p_value = stats::median(.data$p_value),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      odds_ratio_corrected = .data$odds_ratio,
      # Pooled sets restate cells already counted in their members, so they are
      # kept out of the multiplicity family rather than inflating it.
      p_adj = {
        adjusted <- rep(NA_real_, length(.data$p_value))
        keep <- !.data$pooled
        adjusted[keep] <- stats::p.adjust(.data$p_value[keep], "BH")
        adjusted
      }
    ) |>
    dplyr::arrange(dplyr::desc(.data$odds_ratio))
}

#' Per-cluster odds ratio of carrying a disease-enriched clone.
#'
#' The stacked composition bars show where each set of cells sits, but reading
#' an enrichment off two stacked bars is hard because every cluster's height
#' depends on all the others. This states the same comparison directly, as one
#' odds ratio per cluster with a confidence interval, so the question "is this
#' cluster over-represented among disease-enriched clones" has a number and an
#' uncertainty attached to it.
#'
#' The CD8TEM clusters are additionally tested pooled. They are one compartment
#' biologically, and pooling raises the cell count behind the estimate.
#'
#' @param composition Output of `reverse_phenotype_cluster_composition()`.
#' @param pool Named list of cluster sets to test as single units, in addition
#'   to every individual cluster.
#' @return Tibble of one row per cluster and per pooled set.
reverse_phenotype_cluster_enrichment <- function(
    composition,
    pool = NULL) {
  stopifnot(
    is.data.frame(composition),
    all(c("clone_group", "cluster", "n_cells") %in% colnames(composition)),
    "disease_enriched" %in% composition$clone_group
  )

  counts <- composition |>
    dplyr::mutate(cluster = as.character(.data$cluster)) |>
    dplyr::select("clone_group", "cluster", "n_cells")

  totals <- counts |>
    dplyr::group_by(.data$clone_group) |>
    dplyr::summarise(total = sum(.data$n_cells), .groups = "drop")
  total_of <- function(group) {
    value <- totals$total[totals$clone_group == group]
    if (length(value) == 0L) 0L else value
  }
  total_enriched <- total_of("disease_enriched")
  total_matched <- total_of("size_matched")
  stopifnot(total_enriched > 0L, total_matched > 0L)

  cells_in <- function(clusters, group) {
    sum(counts$n_cells[counts$cluster %in% clusters &
      counts$clone_group == group])
  }

  sets <- c(
    stats::setNames(as.list(sort(unique(counts$cluster))), sort(unique(counts$cluster))),
    pool
  )

  purrr::imap_dfr(sets, function(clusters, label) {
    in_enriched <- cells_in(clusters, "disease_enriched")
    in_matched <- cells_in(clusters, "size_matched")
    test <- stats::fisher.test(matrix(
      c(
        in_enriched, total_enriched - in_enriched,
        in_matched, total_matched - in_matched
      ),
      nrow = 2L
    ))
    # A cluster absent from one arm gives an infinite or zero odds ratio, which
    # cannot go on a log axis. CD8TEM_3 is exactly that case, so carry a
    # Haldane-Anscombe corrected estimate alongside the exact one for plotting.
    corrected <- ((in_enriched + 0.5) / (total_enriched - in_enriched + 0.5)) /
      ((in_matched + 0.5) / (total_matched - in_matched + 0.5))
    tibble::tibble(
      cluster = label,
      pooled = label %in% names(pool),
      n_enriched = in_enriched,
      n_matched = in_matched,
      total_enriched = total_enriched,
      total_matched = total_matched,
      fraction_enriched = in_enriched / total_enriched,
      fraction_matched = in_matched / total_matched,
      odds_ratio = unname(test$estimate),
      odds_ratio_corrected = corrected,
      ci_low = test$conf.int[[1L]],
      ci_high = test$conf.int[[2L]],
      p_value = test$p.value
    )
  }) |>
    # Pooled sets restate cells already counted in their members, so they are
    # excluded from the multiplicity correction rather than inflating it.
    dplyr::mutate(p_adj = ifelse(
      .data$pooled, NA_real_,
      stats::p.adjust(ifelse(.data$pooled, NA_real_, .data$p_value), "BH")
    )) |>
    dplyr::arrange(dplyr::desc(.data$odds_ratio))
}

#' Cells eligible for the patient-paired reverse-phenotype contrast.
#'
#' Returns the per-patient arm sizes plus a `retained` flag, so the power
#' statement in the reply is backed by the actual table rather than asserted.
reverse_phenotype_cohort <- function(sc_tcr, clone_set,
                                     cluster_name = "CD8TEM_3",
                                     min_cells_per_arm = 10L) {
  stopifnot(
    inherits(sc_tcr, "Seurat"),
    all(c("CTaa", "selected") %in% colnames(clone_set)),
    length(cluster_name) >= 1L
  )
  enriched <- clone_set$CTaa[clone_set$selected]
  stopifnot(length(enriched) > 0L)

  counts <- sc_tcr[[]] |>
    dplyr::filter(
      !is.na(.data$CTaa),
      as.character(.data$cluster) %in% cluster_name
    ) |>
    dplyr::mutate(clonal_status = ifelse(
      .data$CTaa %in% enriched, "disease_enriched", "other"
    )) |>
    dplyr::count(.data$patient, .data$diagnosis, .data$clonal_status, name = "n_cells") |>
    tidyr::pivot_wider(
      names_from = "clonal_status", values_from = "n_cells", values_fill = 0L
    )

  # pivot_wider only creates the columns it saw, so a cohort with no enriched
  # clones in this cluster would otherwise fail on the retained comparison.
  for (column in c("disease_enriched", "other")) {
    if (!column %in% names(counts)) {
      counts[[column]] <- 0L
    }
  }

  counts |>
    dplyr::mutate(retained = .data$disease_enriched >= min_cells_per_arm &
      .data$other >= min_cells_per_arm) |>
    dplyr::arrange(dplyr::desc(.data$retained), dplyr::desc(.data$disease_enriched))
}

#' Patient-paired pseudobulk contrast of clonal versus non-clonal cells.
#'
#' Cell-level FindMarkers on this comparison recovers CD8 effector genes, which
#' is circular, and its p-values follow whichever patient contributed the most
#' cells. Aggregating to patient x clonal-status and blocking on patient removes
#' both problems. Same Libra -> edgeR -> limma-voom path as run_deg_limma().
reverse_phenotype_pseudobulk <- function(sc_tcr, clone_set,
                                         cluster_name = "CD8TEM_3",
                                         min_cells_per_arm = 10L) {
  cohort <- reverse_phenotype_cohort(
    sc_tcr, clone_set, cluster_name, min_cells_per_arm
  )
  retained <- cohort$patient[cohort$retained]

  if (length(retained) < 3L) {
    return(list(
      results = tibble::tibble(),
      cohort = cohort,
      n_patients = length(retained),
      status = "insufficient_paired_patients"
    ))
  }

  enriched <- clone_set$CTaa[clone_set$selected]
  metadata <- sc_tcr[[]]
  keep <- !is.na(metadata$CTaa) &
    as.character(metadata$cluster) %in% cluster_name &
    metadata$patient %in% retained

  subset_object <- subset(sc_tcr, cells = colnames(sc_tcr)[keep])
  subset_object$clonal_status <- ifelse(
    subset_object$CTaa %in% enriched, "disease_enriched", "other"
  )
  # Libra keys the pseudobulk on cell_type_col, so collapse the requested
  # clusters into one bucket.
  unit_label <- paste(cluster_name, collapse = "+")
  subset_object$reverse_unit <- unit_label

  pseudobulk <- Libra::to_pseudobulk(
    subset_object,
    cell_type_col = "reverse_unit",
    label_col = "clonal_status",
    replicate_col = "patient"
  )
  stopifnot(is.list(pseudobulk), unit_label %in% names(pseudobulk))
  counts <- pseudobulk[[unit_label]]

  design_data <- tibble::tibble(column = colnames(counts)) |>
    dplyr::mutate(
      patient = sub(":.+", "", .data$column),
      status = sub(".+:", "", .data$column)
    ) |>
    dplyr::filter(.data$patient %in% retained)
  # Only patients contributing both arms can inform a blocked contrast.
  complete <- design_data |>
    dplyr::count(.data$patient) |>
    dplyr::filter(.data$n == 2L) |>
    dplyr::pull(.data$patient)
  design_data <- dplyr::filter(design_data, .data$patient %in% complete)

  if (length(complete) < 3L) {
    return(list(
      results = tibble::tibble(),
      cohort = cohort,
      n_patients = length(complete),
      status = "insufficient_paired_patients"
    ))
  }

  counts <- counts[, design_data$column, drop = FALSE]
  design_data <- dplyr::mutate(
    design_data,
    patient = factor(.data$patient),
    status = factor(.data$status, levels = c("other", "disease_enriched"))
  )

  # Groups are defined by clonotype, so the clones' own V, D, J and C segments
  # separate the arms by construction. Leaving them in makes TRBV7-2 the top
  # hit, which is circular rather than biological.
  receptor_genes <- grepl("^(TR[ABGD][VDJC]|IG[HKL][VDJC])", rownames(counts))
  counts <- counts[!receptor_genes, , drop = FALSE]

  dge <- edgeR::DGEList(counts = counts)
  keep_genes <- rowSums(edgeR::cpm(dge) > 1) > 2
  stopifnot(any(keep_genes))
  dge <- edgeR::calcNormFactors(dge[keep_genes, ], method = "TMM")

  design <- stats::model.matrix(~ patient + status, data = design_data)
  coefficient <- "statusdisease_enriched"
  stopifnot(coefficient %in% colnames(design))

  fit <- limma::voomWithQualityWeights(dge, design, plot = FALSE) |>
    limma::lmFit(design = design, block = NULL) |>
    limma::eBayes(robust = TRUE)

  results <- limma::topTable(
    fit, coef = coefficient, n = Inf, adjust.method = "BH"
  ) |>
    tibble::rownames_to_column("gene") |>
    tibble::as_tibble() |>
    dplyr::arrange(dplyr::desc(logFC)) |>
    dplyr::rename(avg_log2FC = logFC, p_val_adj = adj.P.Val)

  list(
    results = results,
    cohort = cohort,
    n_patients = length(complete),
    status = "modeled"
  )
}

#' Patient-paired comparison on a pre-specified gene panel.
#'
#' The honest fallback when too few patients carry both arms for a
#' transcriptome-wide fit, and worth reporting either way: a null on a panel
#' fixed in advance says more at this n than a genome-wide scan does.
reverse_phenotype_panel_scores <- function(sc_tcr, clone_set,
                                           genes = reverse_phenotype_gene_panel(),
                                           cluster_name = "CD8TEM_3",
                                           min_cells_per_arm = 10L) {
  cohort <- reverse_phenotype_cohort(
    sc_tcr, clone_set, cluster_name, min_cells_per_arm
  )
  retained <- cohort$patient[cohort$retained]
  enriched <- clone_set$CTaa[clone_set$selected]

  available <- intersect(genes, rownames(sc_tcr))
  stopifnot(length(available) > 0L)

  metadata <- sc_tcr[[]] |>
    tibble::rownames_to_column("cell_id")
  keep <- !is.na(metadata$CTaa) &
    as.character(metadata$cluster) %in% cluster_name &
    metadata$patient %in% retained

  if (!any(keep) || length(retained) < 3L) {
    return(list(
      scores = tibble::tibble(),
      tests = tibble::tibble(),
      cohort = cohort,
      status = "insufficient_paired_patients"
    ))
  }

  cells <- metadata$cell_id[keep]
  expression <- Seurat::GetAssayData(sc_tcr, layer = "data")[available, cells, drop = FALSE]

  cell_annotation <- metadata |>
    dplyr::filter(.data$cell_id %in% cells) |>
    dplyr::mutate(clonal_status = ifelse(
      .data$CTaa %in% enriched, "disease_enriched", "other"
    )) |>
    dplyr::select("cell_id", "patient", "diagnosis", "clonal_status")
  # Align explicitly rather than trusting that filter() preserved column order.
  expression <- expression[, cell_annotation$cell_id, drop = FALSE]

  scores <- cell_annotation |>
    dplyr::bind_cols(as.data.frame(as.matrix(t(expression)))) |>
    tidyr::pivot_longer(
      dplyr::all_of(available), names_to = "gene", values_to = "expression"
    ) |>
    dplyr::group_by(.data$patient, .data$diagnosis, .data$clonal_status, .data$gene) |>
    dplyr::summarise(
      mean_expression = mean(.data$expression),
      n_cells = dplyr::n(),
      .groups = "drop"
    )

  paired <- scores |>
    dplyr::select(-"n_cells") |>
    tidyr::pivot_wider(
      names_from = "clonal_status", values_from = "mean_expression"
    ) |>
    tidyr::drop_na("disease_enriched", "other") |>
    dplyr::mutate(difference = .data$disease_enriched - .data$other)

  tests <- paired |>
    dplyr::group_by(.data$gene) |>
    dplyr::filter(dplyr::n() >= 3L) |>
    dplyr::summarise(
      n_patients = dplyr::n(),
      mean_difference = mean(.data$difference),
      p_value = tryCatch(
        stats::wilcox.test(
          .data$disease_enriched, .data$other, paired = TRUE, exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(p_adj = stats::p.adjust(.data$p_value, "BH")) |>
    dplyr::arrange(.data$p_value)

  list(scores = scores, paired = paired, tests = tests, cohort = cohort,
       status = "modeled")
}

write_reverse_phenotype_workbook <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(data, path)
  path
}

write_disease_enriched_clone_plot <- function(clone_set, colors, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data <- clone_set |>
    dplyr::mutate(
      neg_log10_fdr = -log10(pmax(.data$p_adj, 1e-300)),
      primary_diagnosis = sub(",.*", "", .data$diagnoses)
    )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$n_target, y = .data$neg_log10_fdr)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$primary_diagnosis, shape = .data$selected),
      alpha = 0.8
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_manual(values = colors, name = "Diagnosis") +
    ggplot2::scale_shape_manual(
      values = c("TRUE" = 17, "FALSE" = 1), name = "Selected"
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Cells in GBS/CIDP pool", y = "-log10(FDR)",
      title = "Per-clone selection tracks clone size",
      subtitle = "Significance rises with size; every clone rests on one patient"
    )
  ggplot2::ggsave(path, plot = plot, width = 7, height = 4.5)
  path
}

#' How many patients support each selected specificity group.
#'
#' This is the panel that answers the reviewer's objection directly: the
#' selection no longer rests on clones private to one person.
write_specificity_group_plot <- function(clone_set, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  stopifnot(
    all(c("gliph_clusters", "group_patients", "n_cells") %in% colnames(clone_set))
  )

  groups <- clone_set |>
    dplyr::group_by(.data$gliph_clusters) |>
    dplyr::summarise(
      group_patients = max(.data$group_patients),
      n_clonotypes = dplyr::n(),
      n_cells = sum(.data$n_cells),
      .groups = "drop"
    )

  plot <- ggplot2::ggplot(
    groups,
    ggplot2::aes(x = .data$group_patients, y = .data$n_clonotypes)
  ) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n_cells), alpha = 0.75) +
    ggplot2::scale_size_continuous(name = "Cells") +
    ggplot2::scale_x_continuous(breaks = scales::breaks_width(1)) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Patients contributing to the specificity group",
      y = "Clonotypes in the group",
      title = "Patient support for enriched specificity groups",
      subtitle = "Selection unit is the group, not the individual clonotype"
    )
  ggplot2::ggsave(path, plot = plot, width = 7, height = 4.5)
  path
}

write_clone_cluster_composition_plot <- function(composition, colors, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  plot <- ggplot2::ggplot(
    composition,
    ggplot2::aes(x = .data$clone_group, y = .data$fraction, fill = .data$cluster)
  ) +
    ggplot2::geom_col(width = 0.7, color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = NULL, y = "Fraction of cells",
      title = "Cluster composition of clones",
      subtitle = "Disease-enriched vs size-matched clones"
    )
  ggplot2::ggsave(path, plot = plot, width = 7, height = 4.5)
  path
}

#' Forest plot of the per-cluster enrichment odds ratios.
#'
#' The companion to the stacked bars. Reads left to right as depleted to
#' enriched, with the pooled CD8TEM estimate marked separately.
write_cluster_enrichment_plot <- function(
    enrichment, path,
    highlight = c("CD8TEM_pooled", "CD8TEM_1", "CD8TEM_2", "CD8TEM_3", "CD4TEM"),
    highlight_colors = c(
      CD8TEM_pooled = "#1F78B4", CD8TEM_1 = "#E41A1C", CD8TEM_2 = "#FB9A99",
      CD8TEM_3 = "#6A3D9A", CD4TEM = "#33A02C"
    )) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  stopifnot("odds_ratio_corrected" %in% colnames(enrichment))
  data <- enrichment |>
    # Plot the continuity-corrected estimate: CD8TEM_3 has no size-matched
    # cells, so its exact odds ratio is infinite and would drop the one cluster
    # this panel most needs to show. Open bounds are drawn as arrows.
    dplyr::filter(.data$n_enriched > 0L) |>
    dplyr::arrange(dplyr::desc(.data$odds_ratio_corrected)) |>
    dplyr::mutate(
      cluster = factor(.data$cluster, levels = rev(.data$cluster)),
      unbounded = !is.finite(.data$odds_ratio),
      ci_high_plot = ifelse(
        is.finite(.data$ci_high), .data$ci_high,
        max(.data$odds_ratio_corrected) * 3
      ),
      ci_low_plot = pmax(.data$ci_low, min(.data$odds_ratio_corrected) / 3),
      role = dplyr::if_else(
        as.character(.data$cluster) %in% highlight,
        as.character(.data$cluster), "other"
      ),
      # Pooled sets are held out of the multiplicity family and carry
      # p_adj = NA, so fall back to the interval for those rows only.
      significant = dplyr::if_else(
        is.na(.data$p_adj),
        .data$ci_low > 1 | .data$ci_high < 1,
        .data$p_adj < 0.05
      )
    )

  plot <- ggplot2::ggplot(
    data, ggplot2::aes(x = .data$odds_ratio_corrected, y = .data$cluster)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$ci_low_plot, xmax = .data$ci_high_plot),
      height = 0.25, colour = "grey45"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = .data$odds_ratio_corrected, colour = .data$role,
        shape = .data$significant),
      size = 2.8
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(
      values = c(highlight_colors[highlight], other = "grey35"),
      breaks = c(highlight, "other"),
      name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1), name = "FDR < 0.05"
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Odds ratio, disease-enriched vs size-matched clones (log scale)",
      y = NULL,
      title = "Which clusters carry the disease-enriched clones",
      subtitle = paste(
        "Right of the line means over-represented among disease-enriched",
        "clones.\nEstimated over repeated size-matched draws"
      )
    )
  ggplot2::ggsave(path, plot = plot, width = 7.5, height = 5)
  path
}

write_reverse_phenotype_plot <- function(pseudobulk, path, n_label = 15L,
                                         seed = 42L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  results <- pseudobulk$results

  if (nrow(results) == 0L) {
    plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = 0, y = 0,
        label = paste(
          "Not modeled:", pseudobulk$status,
          "\npatients with both arms:", pseudobulk$n_patients
        )
      ) +
      ggplot2::theme_void()
    ggplot2::ggsave(path, plot = plot, width = 6, height = 4.5)
    return(path)
  }

  data <- dplyr::mutate(
    results,
    significant = .data$p_val_adj < 0.05,
    neg_log10_p = -log10(pmax(.data$P.Value, 1e-300))
  )
  labels <- data |>
    dplyr::filter(.data$significant) |>
    dplyr::slice_max(abs(.data$avg_log2FC), n = n_label, with_ties = FALSE)

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = .data$avg_log2FC, y = .data$neg_log10_p)
  ) +
    ggplot2::geom_point(ggplot2::aes(color = .data$significant), alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "grey70"), guide = "none"
    ) +
    ggrepel::geom_text_repel(
      data = labels, ggplot2::aes(label = .data$gene),
      size = 2.6, max.overlaps = 20, seed = seed
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "log2FC (disease-enriched clones vs other cells)",
      y = "-log10(p)",
      title = "Reverse phenotyping of disease-enriched clones",
      subtitle = paste0(
        "Patient-paired pseudobulk, n = ", pseudobulk$n_patients, " patients"
      )
    )
  ggplot2::ggsave(path, plot = plot, width = 6, height = 4.5)
  path
}

write_reverse_phenotype_panel_plot <- function(panel, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (nrow(panel$paired) == 0L) {
    plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = paste("Not modeled:", panel$status)) +
      ggplot2::theme_void()
    ggplot2::ggsave(path, plot = plot, width = 7, height = 5)
    return(path)
  }

  order <- panel$tests |>
    dplyr::arrange(dplyr::desc(.data$mean_difference)) |>
    dplyr::pull(.data$gene)
  data <- dplyr::mutate(panel$paired, gene = factor(.data$gene, levels = order))

  plot <- ggplot2::ggplot(
    data, ggplot2::aes(x = .data$gene, y = .data$difference)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey92") +
    ggplot2::geom_jitter(width = 0.15, size = 1.3, alpha = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      x = NULL,
      y = "Mean expression difference (enriched - other)",
      title = "Pre-specified panel, patient-paired",
      subtitle = "Effect sizes only; the signed-rank p-value has no resolution at this n"
    )
  ggplot2::ggsave(path, plot = plot, width = 8.5, height = 4.5)
  path
}
