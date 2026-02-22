# ==============================================================================
# TCR GLIPH2 Manuscript Figures
# ==============================================================================
#
# Reads cached analysis results from data/cache/gliph_results.qs and generates
# ~20 individual standalone PDF figures for the manuscript.
#
# Depends on: scripts/tcr_comparison.R (produces the cache)
#
# Figures (all individual standalone PDFs):
#   Disease-specific TCR convergence:
#     - disease_enrichment_volcano, cluster_composition_heatmap,
#       cluster_specificity_private_shared, cdr3_length_enriched_vs_background,
#       vgene_usage_enriched_clusters, public_clone_rate_by_diagnosis,
#       cdr3_motif_enrichment_by_diagnosis
#   GBS CSF tropism + Sukenikova cross-validation:
#     - tissue_compartmentalization_scatter, cross_cohort_cluster_sharing_alluvial,
#       sukenikova_acute_vs_recovery, cross_cohort_top_cluster_composition
#   Myelin reactivity:
#     - myelin_antigen_reactivity_dotplot, myelin_cluster_network,
#       myelin_exact_hamming_match_table
#   Supplementary:
#     - enriched_clusters_overview_bubble, diagnosis_cluster_sharing_heatmap,
#       tissue_bias_fraction_by_diagnosis, csf_vs_pbmc_motif_enrichment,
#       myelin_reactive_fraction_by_diagnosis, sukenikova_shared_clusters_by_diagnosis
#
# ==============================================================================

# ==============================================================================
# 0. Setup: Libraries, Cache Loading, Theme, Helpers
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(tidyverse)
  library(ggrepel)
  library(viridis)
  library(circlize)
  library(gridExtra)
  library(grid)
})

set.seed(42)

# -- Load cache ---------------------------------------------------------------
message("Loading GLIPH analysis cache...")
cache_path <- "data/cache/gliph_results.qs"
stopifnot(file.exists(cache_path))

gliph_cache <- qs::qread(cache_path, nthreads = 6)
list2env(gliph_cache, envir = environment())
message(sprintf("  Loaded %d cached objects (created %s)",
                length(gliph_cache), cache_timestamp))

# -- Output directory ---------------------------------------------------------
fig_dir <- "results/tcr_comparison/manuscript_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# -- Publication theme --------------------------------------------------------
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "grey95")
    )
}

# -- Helper: format p-values for annotation -----------------------------------
format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

# -- Helper: safe ggsave with message -----------------------------------------
safe_save <- function(filename, plot, width, height) {
  out <- file.path(fig_dir, filename)
  ggsave(out, plot = plot, width = width, height = height, device = "pdf")
  message(sprintf("  Saved: %s (%g x %g)", filename, width, height))
}

message("\n============================================================")
message("Generating manuscript figures")
message("============================================================\n")


# ==============================================================================
# STORY A: GLIPH Enrichment, Composition, CDR3 Features
# ==============================================================================

# ------------------------------------------------------------------------------
# A1: Enrichment Volcano (faceted by diagnosis)
# ------------------------------------------------------------------------------
message("--- Figure A1: Enrichment volcano ---")

if (exists("diagnosis_enrichment") && nrow(diagnosis_enrichment) > 0 &&
    exists("tissue_enrichment") && nrow(tissue_enrichment) > 0) {

  volcano_data <- diagnosis_enrichment |>
    dplyr::filter(diagnosis %in% neuropathy_dx) |>
    dplyr::left_join(
      tissue_enrichment |> dplyr::select(cluster, tissue_bias),
      by = "cluster"
    ) |>
    dplyr::mutate(
      tissue_bias = tidyr::replace_na(tissue_bias, "No bias"),
      neg_log10_fdr = -log10(pmax(p.adj, 1e-300)),
      log2_OR_plot = log2(pmax(odds_ratio, 1e-4) + 0.01),
      significant = p.adj < 0.05 & abs(log2_OR_plot) > 1
    )

  # Label top 3 significant per diagnosis
  top_labels_a1 <- volcano_data |>
    dplyr::filter(significant) |>
    dplyr::group_by(diagnosis) |>
    dplyr::slice_max(order_by = neg_log10_fdr, n = 3, with_ties = FALSE) |>
    dplyr::ungroup()

  tissue_bias_colors <- c(
    "CSF-enriched"  = "#E41A1C",
    "PBMC-enriched" = "#377EB8",
    "No bias"       = "grey60"
  )

  p_a1 <- ggplot(volcano_data,
                 aes(x = log2_OR_plot, y = neg_log10_fdr)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_point(aes(color = tissue_bias, size = cluster_size), alpha = 0.65) +
    geom_text_repel(
      data = top_labels_a1,
      aes(label = paste0("C", cluster)),
      size = 2.5, max.overlaps = 15, segment.color = "grey40",
      segment.size = 0.3, fontface = "italic", seed = 42
    ) +
    facet_wrap(~ diagnosis, ncol = 3, scales = "free") +
    scale_color_manual(values = tissue_bias_colors, name = "Tissue bias") +
    scale_size_continuous(range = c(0.8, 5), name = "Cluster size") +
    theme_pub() +
    labs(
      title = "Disease enrichment of GLIPH2 specificity groups",
      x = expression(log[2](Odds~Ratio)),
      y = expression(-log[10](FDR))
    )

  safe_save("fig_gliph_disease_enrichment_volcano.pdf", p_a1, 10, 7)
} else {
  message("  SKIP: insufficient data for A1")
}


# ------------------------------------------------------------------------------
# A2: Cluster Composition Heatmap (ComplexHeatmap)
# ------------------------------------------------------------------------------
message("--- Figure A2: Cluster composition heatmap ---")

if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    exists("cluster_composition") && nrow(cluster_composition) > 0 &&
    exists("enriched_cluster_ids") && length(enriched_cluster_ids) > 0) {

  hm_data <- cluster_composition |>
    dplyr::filter(cluster %in% enriched_cluster_ids) |>
    dplyr::select(cluster, diagnosis, fraction) |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = fraction, values_fill = 0) |>
    tibble::column_to_rownames("cluster") |>
    as.matrix()

  # Ensure all dx_order columns present
  for (dx in dx_order) {
    if (!dx %in% colnames(hm_data)) hm_data <- cbind(hm_data, setNames(rep(0, nrow(hm_data)), dx))
  }
  hm_data <- hm_data[, intersect(dx_order, colnames(hm_data)), drop = FALSE]

  # Row annotations
  row_info <- cluster_entropy |>
    dplyr::filter(cluster %in% as.integer(rownames(hm_data))) |>
    dplyr::arrange(match(cluster, as.integer(rownames(hm_data))))

  if (nrow(row_info) > 0 && nrow(row_info) == nrow(hm_data)) {
    spec_colors <- c("Private" = "#1B9E77", "Semi-private" = "#D95F02", "Promiscuous" = "#7570B3")
    ha_right <- ComplexHeatmap::rowAnnotation(
      Size = ComplexHeatmap::anno_barplot(row_info$cluster_total, width = unit(2, "cm")),
      Specificity = row_info$specificity_class,
      col = list(Specificity = spec_colors),
      annotation_name_side = "top"
    )
  } else {
    ha_right <- NULL
  }

  col_fun <- circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                                  c("white", "#FEE5D9", "#FCAE91", "#FB6A4A", "#CB181D"))

  pdf(file.path(fig_dir, "fig_gliph_cluster_composition_heatmap.pdf"), width = 12, height = 10)
  ht <- ComplexHeatmap::Heatmap(
    hm_data,
    name = "Fraction",
    col = col_fun,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    row_title = "Enriched GLIPH2 clusters",
    column_title = "Cluster composition by diagnosis",
    right_annotation = ha_right,
    column_names_rot = 45,
    border = TRUE,
    heatmap_legend_param = list(title = "Fraction")
  )
  ComplexHeatmap::draw(ht)
  dev.off()
  message(sprintf("  Saved: fig_gliph_cluster_composition_heatmap.pdf (12 x 10)"))
}


# ------------------------------------------------------------------------------
# A3: Private vs Shared Clusters (stacked bar)
# ------------------------------------------------------------------------------
message("--- Figure A3: Private vs shared clusters ---")

if (exists("specificity_summary") && nrow(specificity_summary) > 0) {

  spec_long <- specificity_summary |>
    tidyr::pivot_longer(
      cols = -enriched_for,
      names_to = "specificity_class",
      values_to = "n_clusters"
    ) |>
    dplyr::mutate(
      specificity_class = factor(specificity_class,
                                 levels = c("Private", "Semi-private", "Promiscuous")),
      enriched_for = factor(enriched_for, levels = intersect(c(dx_order, "GBS_Sukenikova"),
                                                              unique(enriched_for)))
    )

  spec_colors <- c("Private" = "#1B9E77", "Semi-private" = "#D95F02", "Promiscuous" = "#7570B3")

  p_a3 <- ggplot(spec_long, aes(x = enriched_for, y = n_clusters, fill = specificity_class)) +
    geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = spec_colors, name = "Specificity class") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Disease-enriched cluster specificity",
      x = "Enriched diagnosis",
      y = "Number of clusters"
    )

  safe_save("fig_gliph_cluster_specificity_private_shared.pdf", p_a3, 8, 5.5)
} else {
  message("  SKIP: no specificity_summary data")
}


# ------------------------------------------------------------------------------
# A4: CDR3 Length Constraint (density, faceted by diagnosis)
# ------------------------------------------------------------------------------
message("--- Figure A4: CDR3 length constraint ---")

if (exists("enriched_cdr3_lengths") && nrow(enriched_cdr3_lengths) > 0 &&
    exists("bg_cdr3_lengths") && nrow(bg_cdr3_lengths) > 0) {

  bg_for_plot <- bg_cdr3_lengths |>
    dplyr::select(cdr3_length) |>
    dplyr::mutate(group = "Background")

  enriched_for_plot <- enriched_cdr3_lengths |>
    dplyr::filter(!is.na(enriched_for)) |>
    dplyr::select(cdr3_length, enriched_for)

  dx_with_data <- unique(enriched_for_plot$enriched_for)

  if (length(dx_with_data) > 0) {
    # Build combined data with background duplicated per facet
    combined_density <- dplyr::bind_rows(
      lapply(dx_with_data, function(dx) {
        dplyr::bind_rows(
          bg_for_plot |> dplyr::mutate(enriched_for = dx),
          enriched_for_plot |>
            dplyr::filter(enriched_for == dx) |>
            dplyr::mutate(group = "Enriched")
        )
      })
    )

    # Add KS test annotations
    ks_annotations <- cdr3_length_tests |>
      dplyr::filter(diagnosis %in% dx_with_data) |>
      dplyr::mutate(
        label = paste0("KS p", ifelse(ks_padj < 0.001, "<0.001",
                                       sprintf("=%.3f", ks_padj))),
        enriched_for = diagnosis
      )

    p_a4 <- ggplot(combined_density, aes(x = cdr3_length, fill = group, color = group)) +
      geom_density(alpha = 0.35, linewidth = 0.6) +
      facet_wrap(~ enriched_for, ncol = 3, scales = "free_y") +
      scale_fill_manual(values = c("Background" = "grey70", "Enriched" = "#E41A1C"),
                        name = "") +
      scale_color_manual(values = c("Background" = "grey50", "Enriched" = "#C10020"),
                         name = "") +
      geom_text(
        data = ks_annotations,
        aes(x = Inf, y = Inf, label = label),
        inherit.aes = FALSE,
        hjust = 1.1, vjust = 1.5, size = 3, fontface = "italic"
      ) +
      theme_pub() +
      labs(
        title = "CDR3 length distributions: enriched clusters vs background",
        x = "CDR3 length (amino acids)",
        y = "Density"
      )

    safe_save("fig_gliph_cdr3_length_enriched_vs_background.pdf", p_a4, 9, 6)
  }
} else {
  message("  SKIP: no CDR3 length data")
}


# ------------------------------------------------------------------------------
# A5: V-Gene Within Clusters Dot Plot
# ------------------------------------------------------------------------------
message("--- Figure A5: V-gene within clusters dot plot ---")

if (exists("vgene_enrichment_tests") && nrow(vgene_enrichment_tests) > 0) {

  vgene_plot <- vgene_enrichment_tests |>
    dplyr::filter(!is.na(fold_enrichment), !is.na(n_in_clusters)) |>
    dplyr::mutate(
      log2_fe = log2(fold_enrichment + 1e-4),
      sig_label = sapply(p.adj, format_pval)
    )

  # Keep V-genes with at least one significant enrichment or top 15 by count
  sig_vgenes <- vgene_plot |>
    dplyr::filter(p.adj < 0.1) |>
    dplyr::pull(TRBV) |> unique()

  top_vgenes <- vgene_plot |>
    dplyr::group_by(TRBV) |>
    dplyr::summarize(total_n = sum(n_in_clusters), .groups = "drop") |>
    dplyr::slice_max(order_by = total_n, n = 15) |>
    dplyr::pull(TRBV)

  keep_vgenes <- unique(c(sig_vgenes, top_vgenes))

  vgene_plot_filt <- vgene_plot |>
    dplyr::filter(TRBV %in% keep_vgenes)

  if (nrow(vgene_plot_filt) > 0) {
    p_a5 <- ggplot(vgene_plot_filt,
                   aes(x = enriched_for, y = TRBV,
                       size = n_in_clusters, color = log2_fe)) +
      geom_point() +
      geom_text(aes(label = sig_label), size = 3, vjust = -0.8, color = "black") +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                            midpoint = 0, name = expression(log[2](Fold~Enrich.))) +
      scale_size_continuous(range = c(1, 8), name = "n in clusters") +
      theme_pub() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
      ) +
      labs(
        title = "V-gene usage within disease-enriched GLIPH2 clusters",
        x = "Enriched diagnosis",
        y = "TRBV gene"
      )

    safe_save("fig_gliph_vgene_usage_enriched_clusters.pdf", p_a5, 10, 8)
  }
} else {
  message("  SKIP: no V-gene enrichment data")
}


# ------------------------------------------------------------------------------
# A6: Public Clone Rates
# ------------------------------------------------------------------------------
message("--- Figure A6: Public clone rates ---")

if (exists("public_clone_rates") && nrow(public_clone_rates) > 0) {

  pcr_plot <- public_clone_rates |>
    dplyr::mutate(
      public_pct = public_rate * 100,
      diagnosis = factor(diagnosis, levels = intersect(c(dx_order, "GBS_Sukenikova"),
                                                        diagnosis))
    )

  p_a6 <- ggplot(pcr_plot, aes(x = diagnosis, y = public_pct, fill = diagnosis)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%", public_pct)),
              vjust = -0.5, size = 3) +
    scale_fill_manual(values = diagnosis_col) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Public clone rate by diagnosis",
      x = "",
      y = "Public clone rate (%)",
      subtitle = "Fraction of unique CDR3b shared by 2+ patients"
    ) +
    expand_limits(y = max(pcr_plot$public_pct, na.rm = TRUE) * 1.15)

  safe_save("fig_gliph_public_clone_rate_by_diagnosis.pdf", p_a6, 8, 5.5)
} else {
  message("  SKIP: no public clone rate data")
}


# ------------------------------------------------------------------------------
# A7: Motif Enrichment (faceted horizontal bar, top 5 per dx)
# ------------------------------------------------------------------------------
message("--- Figure A7: Motif enrichment faceted ---")

if (exists("motif_driving") && !is.null(motif_driving) && nrow(motif_driving) > 0) {

  motif_top <- motif_driving |>
    dplyr::filter(!is.na(log2_fc), !is.na(p.adj)) |>
    dplyr::group_by(diagnosis) |>
    dplyr::slice_max(order_by = abs(log2_fc), n = 5, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      sig_label = sapply(p.adj, format_pval),
      motif_label = paste0(motif, " (", k, "-mer)")
    )

  dx_with_motifs <- unique(motif_top$diagnosis)

  if (length(dx_with_motifs) > 0) {
    # Order motifs within each facet by log2_fc
    motif_top <- motif_top |>
      dplyr::arrange(diagnosis, log2_fc) |>
      dplyr::mutate(motif_label = forcats::fct_inorder(motif_label))

    p_a7 <- ggplot(motif_top,
                   aes(x = log2_fc, y = motif_label, fill = diagnosis)) +
      geom_col(show.legend = FALSE, width = 0.7) +
      geom_text(aes(label = sig_label),
                hjust = ifelse(motif_top$log2_fc > 0, -0.3, 1.3),
                size = 3) +
      facet_wrap(~ diagnosis, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = diagnosis_col) +
      geom_vline(xintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.4) +
      theme_pub() +
      labs(
        title = "CDR3 motif enrichment in disease-associated GLIPH2 clusters",
        x = expression(log[2](fold~change~vs~background)),
        y = ""
      )

    safe_save("fig_gliph_cdr3_motif_enrichment_by_diagnosis.pdf", p_a7, 12, 8)
  }
} else {
  message("  SKIP: no motif data")
}


# ==============================================================================
# STORY B: Tissue Compartmentalization & Cross-Cohort
# ==============================================================================

# ------------------------------------------------------------------------------
# B1: CSF Fraction vs Enrichment Scatter
# ------------------------------------------------------------------------------
message("--- Figure B1: CSF vs enrichment scatter ---")

if (exists("cluster_tissue_dx") && nrow(cluster_tissue_dx) > 0) {

  scatter_data <- cluster_tissue_dx |>
    dplyr::filter(!is.na(csf_fraction), !is.na(log2_OR)) |>
    dplyr::mutate(
      diagnosis = factor(diagnosis, levels = intersect(c(neuropathy_dx, "GBS_Sukenikova"),
                                                        unique(diagnosis)))
    )

  if (nrow(scatter_data) > 0) {
    # Label outliers
    top_scatter <- scatter_data |>
      dplyr::filter(p.adj < 0.01 | abs(csf_fraction - 0.5) > 0.35) |>
      dplyr::group_by(diagnosis) |>
      dplyr::slice_max(order_by = abs(log2_OR), n = 2, with_ties = FALSE) |>
      dplyr::ungroup()

    p_b1 <- ggplot(scatter_data,
                   aes(x = csf_fraction, y = log2_OR, color = diagnosis)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
      geom_point(aes(size = cluster_size), alpha = 0.7) +
      geom_text_repel(
        data = top_scatter,
        aes(label = paste0("C", cluster)),
        size = 2.5, max.overlaps = 12, segment.color = "grey40",
        segment.size = 0.3, seed = 42
      ) +
      scale_color_manual(values = diagnosis_col, name = "Diagnosis") +
      scale_size_continuous(range = c(1, 6), name = "Cluster size") +
      theme_pub() +
      labs(
        title = "Tissue compartmentalization of disease-enriched clusters",
        x = "CSF fraction",
        y = expression(log[2](Odds~Ratio~vs~CTRL))
      )

    safe_save("fig_gliph_tissue_compartmentalization_scatter.pdf", p_b1, 8, 6)
  }
} else {
  message("  SKIP: no cluster_tissue_dx data")
}


# ------------------------------------------------------------------------------
# B2: Cross-Cohort Alluvial (ggalluvial)
# ------------------------------------------------------------------------------
message("--- Figure B2: Cross-cohort alluvial ---")

if (exists("cross_cohort") && nrow(cross_cohort) > 0 &&
    sum(cross_cohort$cross_cohort) > 0) {

  if (requireNamespace("ggalluvial", quietly = TRUE)) {

    cross_clusters <- cross_cohort |>
      dplyr::filter(cross_cohort) |>
      dplyr::pull(cluster)

    # Build alluvial data: cluster -> Heming dx and Sukenikova patient
    alluvial_data <- cluster_meta |>
      dplyr::filter(cluster %in% cross_clusters, !is.na(source), !is.na(diagnosis)) |>
      dplyr::mutate(
        axis_value = dplyr::if_else(source == "Heming", diagnosis, patient)
      ) |>
      dplyr::count(cluster, source, axis_value) |>
      tidyr::pivot_wider(names_from = source, values_from = axis_value,
                         values_fn = function(x) paste(unique(x), collapse = "/")) |>
      dplyr::filter(!is.na(Heming), !is.na(Sukenikova))

    if (nrow(alluvial_data) > 0) {
      p_b2 <- ggplot(alluvial_data,
                     aes(axis1 = Heming, axis2 = Sukenikova, y = n)) +
        ggalluvial::geom_alluvium(aes(fill = Heming), alpha = 0.6, width = 1/6) +
        ggalluvial::geom_stratum(width = 1/6, fill = "grey90", color = "grey40") +
        geom_text(stat = ggalluvial::StatStratum, aes(label = after_stat(stratum)),
                  size = 2.8) +
        scale_x_discrete(limits = c("Heming diagnosis", "Sukenikova patient"),
                         expand = c(0.15, 0.05)) +
        scale_fill_manual(values = diagnosis_col, name = "Heming diagnosis") +
        theme_pub() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()
        ) +
        labs(
          title = "Cross-cohort GLIPH2 cluster sharing",
          y = "Number of shared TCR sequences"
        )

      safe_save("fig_gliph_cross_cohort_cluster_sharing_alluvial.pdf", p_b2, 10, 7)
    }
  } else {
    message("  SKIP: ggalluvial package not available")
  }
} else {
  message("  SKIP: no cross-cohort data or no cross-cohort clusters")
}


# ------------------------------------------------------------------------------
# B3: Sukenikova Timepoint Comparison (grouped bar)
# ------------------------------------------------------------------------------
message("--- Figure B3: Sukenikova timepoint comparison ---")

if (exists("timepoint_summary") && nrow(timepoint_summary) > 0) {

  tp_plot <- timepoint_summary |>
    dplyr::filter(total_AC > 0 | total_REC > 0) |>
    dplyr::select(diagnosis, total_AC, total_REC) |>
    tidyr::pivot_longer(
      cols = c(total_AC, total_REC),
      names_to = "timepoint",
      values_to = "count"
    ) |>
    dplyr::mutate(
      timepoint = gsub("total_", "", timepoint),
      timepoint = factor(timepoint, levels = c("AC", "REC")),
      diagnosis = factor(diagnosis, levels = intersect(c(dx_order, "GBS_Sukenikova"),
                                                        unique(diagnosis)))
    )

  if (nrow(tp_plot) > 0) {
    tp_colors <- c("AC" = "#FC8D62", "REC" = "#66C2A5")

    p_b3 <- ggplot(tp_plot, aes(x = diagnosis, y = count, fill = timepoint)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_text(
        aes(label = count),
        position = position_dodge(width = 0.7),
        vjust = -0.5, size = 3
      ) +
      scale_fill_manual(values = tp_colors, name = "Timepoint",
                        labels = c("AC" = "Acute", "REC" = "Recovery")) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Sukenikova GBS TCRs in disease-enriched clusters: acute vs recovery",
        x = "Heming diagnosis (enriched)",
        y = "Number of Sukenikova TCR sequences"
      ) +
      expand_limits(y = max(tp_plot$count, na.rm = TRUE) * 1.15)

    safe_save("fig_gliph_sukenikova_acute_vs_recovery.pdf", p_b3, 8, 5.5)
  }
} else {
  message("  SKIP: no timepoint data")
}


# ------------------------------------------------------------------------------
# B4: Cross-Cohort Cluster Profiles (top 10 stacked bar)
# ------------------------------------------------------------------------------
message("--- Figure B4: Cross-cohort cluster profiles ---")

if (exists("cross_cohort") && nrow(cross_cohort) > 0 &&
    exists("cluster_composition") && nrow(cluster_composition) > 0) {

  top10_cross <- cross_cohort |>
    dplyr::filter(cross_cohort) |>
    dplyr::slice_max(order_by = n_total, n = 10, with_ties = FALSE) |>
    dplyr::pull(cluster)

  cc_comp <- cluster_composition |>
    dplyr::filter(cluster %in% top10_cross) |>
    dplyr::mutate(
      cluster_label = paste0("C", cluster),
      diagnosis = factor(diagnosis, levels = c(dx_order, "GBS_Sukenikova"))
    )

  if (nrow(cc_comp) > 0) {
    # Order clusters by total size
    cluster_order <- cc_comp |>
      dplyr::distinct(cluster, cluster_label, cluster_total) |>
      dplyr::arrange(desc(cluster_total)) |>
      dplyr::pull(cluster_label)

    cc_comp$cluster_label <- factor(cc_comp$cluster_label, levels = cluster_order)

    p_b4 <- ggplot(cc_comp, aes(x = cluster_label, y = fraction, fill = diagnosis)) +
      geom_col(width = 0.75, color = "white", linewidth = 0.2) +
      scale_fill_manual(values = diagnosis_col, name = "Diagnosis") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Top 10 cross-cohort GLIPH2 cluster composition",
        x = "Cluster",
        y = "Fraction of cluster"
      )

    safe_save("fig_gliph_cross_cohort_top_cluster_composition.pdf", p_b4, 10, 8)
  }
} else {
  message("  SKIP: no cross-cohort cluster data")
}


# ==============================================================================
# STORY C: Myelin Reactivity
# ==============================================================================

# ------------------------------------------------------------------------------
# C1: Myelin Reactivity Dot Plot (diagnosis x antigen)
# ------------------------------------------------------------------------------
message("--- Figure C1: Myelin reactivity dot plot ---")

if (exists("antigen_cluster_enrichment") && nrow(antigen_cluster_enrichment) > 0) {

  # Aggregate: per diagnosis, per antigen
  myelin_dot_data <- antigen_cluster_enrichment |>
    dplyr::filter(!is.na(dx_enriched), !is.na(antigen)) |>
    dplyr::group_by(dx_enriched, antigen) |>
    dplyr::summarize(
      n_reactive = sum(n_antigen_in_cluster, na.rm = TRUE),
      mean_OR = mean(dx_OR, na.rm = TRUE),
      min_padj = min(dx_padj, na.rm = TRUE),
      n_clusters = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      log2_OR = log2(pmax(mean_OR, 1e-4)),
      sig_label = sapply(min_padj, format_pval)
    )

  if (nrow(myelin_dot_data) > 0) {
    p_c1 <- ggplot(myelin_dot_data,
                   aes(x = antigen, y = dx_enriched,
                       size = n_reactive, color = log2_OR)) +
      geom_point() +
      geom_text(aes(label = sig_label), size = 3, vjust = -1, color = "black") +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                            midpoint = 0, name = expression(log[2](OR))) +
      scale_size_continuous(range = c(2, 10), name = "n reactive TCRs") +
      theme_pub() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3)
      ) +
      labs(
        title = "Myelin antigen reactivity in disease-enriched GLIPH2 clusters",
        x = "Myelin antigen",
        y = "Enriched diagnosis"
      )

    safe_save("fig_gliph_myelin_antigen_reactivity_dotplot.pdf", p_c1, 8, 6)
  }
} else {
  message("  SKIP: no antigen cluster enrichment data")
}


# ------------------------------------------------------------------------------
# C2: Myelin Cluster Network (ggraph)
# ------------------------------------------------------------------------------
message("--- Figure C2: Myelin cluster network ---")

if (exists("myelin_clusters") && length(myelin_clusters) > 0 &&
    exists("cluster_meta") && nrow(cluster_meta) > 0) {

  if (requireNamespace("ggraph", quietly = TRUE) &&
      requireNamespace("tidygraph", quietly = TRUE) &&
      requireNamespace("igraph", quietly = TRUE)) {

    # Build edges: CDR3b sequences shared between myelin clusters
    myelin_cl_data <- cluster_meta |>
      dplyr::filter(cluster %in% myelin_clusters, !is.na(CDR3b), !is.na(diagnosis))

    # Nodes: clusters
    node_info <- myelin_cl_data |>
      dplyr::group_by(cluster) |>
      dplyr::summarize(
        size = dplyr::n(),
        dominant_dx = names(sort(table(diagnosis), decreasing = TRUE))[1],
        has_myelin = any(CDR3b %in% myelin_cdr3_set),
        n_diagnoses = dplyr::n_distinct(diagnosis),
        .groups = "drop"
      ) |>
      dplyr::mutate(name = as.character(cluster))

    # Edges: clusters sharing CDR3 sequences
    shared_cdr3 <- myelin_cl_data |>
      dplyr::select(CDR3b, cluster) |>
      dplyr::distinct()

    if (nrow(shared_cdr3) > 0 && dplyr::n_distinct(shared_cdr3$cluster) > 1) {
      edge_list <- shared_cdr3 |>
        dplyr::inner_join(shared_cdr3, by = "CDR3b", suffix = c("_from", "_to")) |>
        dplyr::filter(cluster_from < cluster_to) |>
        dplyr::count(cluster_from, cluster_to, name = "weight") |>
        dplyr::rename(from = cluster_from, to = cluster_to) |>
        dplyr::mutate(from = as.character(from), to = as.character(to))

      if (nrow(edge_list) > 0) {
        graph <- tidygraph::tbl_graph(
          nodes = node_info,
          edges = edge_list,
          directed = FALSE
        )

        p_c2 <- ggraph::ggraph(graph, layout = "fr") +
          ggraph::geom_edge_link(aes(width = weight), alpha = 0.3, color = "grey60") +
          ggraph::geom_node_point(aes(size = size, color = dominant_dx,
                                       shape = has_myelin)) +
          ggraph::geom_node_text(aes(label = name), size = 2.5, repel = TRUE) +
          scale_color_manual(values = diagnosis_col, name = "Dominant diagnosis") +
          scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16),
                             name = "Contains myelin TCR",
                             labels = c("TRUE" = "Yes", "FALSE" = "No")) +
          ggraph::scale_edge_width_continuous(range = c(0.3, 2.5), name = "Shared CDR3s") +
          scale_size_continuous(range = c(2, 10), name = "Cluster size") +
          theme_void() +
          theme(
            plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            legend.position = "right"
          ) +
          labs(title = "Network of myelin-reactive GLIPH2 clusters")

        safe_save("fig_gliph_myelin_cluster_network.pdf", p_c2, 9, 8)
      } else {
        message("  SKIP: no edges between myelin clusters")
      }
    } else {
      message("  SKIP: insufficient shared CDR3 data for network")
    }
  } else {
    message("  SKIP: ggraph/tidygraph/igraph not available")
  }
} else {
  message("  SKIP: no myelin cluster data")
}


# ------------------------------------------------------------------------------
# C3: Exact and Near Match Table
# ------------------------------------------------------------------------------
message("--- Figure C3: Exact/near match table ---")

has_exact <- exists("exact_matches") && nrow(exact_matches) > 0
has_near  <- exists("near_matches_annotated") && nrow(near_matches_annotated) > 0

if (has_exact || has_near) {

  table_rows <- list()

  if (has_exact) {
    exact_tbl <- exact_matches |>
      dplyr::select(TRB_CDR3, patient, tissue, diagnosis, specificity) |>
      dplyr::mutate(match_type = "Exact") |>
      dplyr::rename(CDR3b = TRB_CDR3) |>
      dplyr::distinct()
    table_rows <- c(table_rows, list(exact_tbl))
  }

  if (has_near) {
    near_tbl <- near_matches_annotated |>
      dplyr::filter(!is.na(diagnosis)) |>
      dplyr::select(heming_CDR3b, myelin_CDR3b, patient, tissue, diagnosis, specificity) |>
      dplyr::mutate(
        match_type = "Hamming=1",
        CDR3b = paste0(heming_CDR3b, " -> ", myelin_CDR3b)
      ) |>
      dplyr::select(CDR3b, patient, tissue, diagnosis, specificity, match_type) |>
      dplyr::distinct()
    table_rows <- c(table_rows, list(near_tbl))
  }

  combined_table <- dplyr::bind_rows(table_rows) |>
    dplyr::arrange(match_type, diagnosis, CDR3b)

  col_names <- c("CDR3b", "Patient", "Tissue", "Diagnosis", "Specificity", "Match type")
  names(combined_table) <- col_names

  tg <- gridExtra::tableGrob(
    combined_table,
    rows = NULL,
    theme = gridExtra::ttheme_default(
      core = list(
        fg_params = list(cex = 0.7),
        bg_params = list(fill = c("white", "#F0F0F0"), alpha = 0.8)
      ),
      colhead = list(
        fg_params = list(cex = 0.8, fontface = "bold"),
        bg_params = list(fill = "#D9D9D9")
      )
    )
  )

  # Calculate dynamic height: header + rows
  tbl_height <- max(4, 0.35 * nrow(combined_table) + 1.5)

  pdf(file.path(fig_dir, "fig_gliph_myelin_exact_hamming_match_table.pdf"),
      width = 12, height = tbl_height)
  grid::grid.newpage()
  grid::grid.draw(tg)
  dev.off()
  message(sprintf("  Saved: fig_gliph_myelin_exact_hamming_match_table.pdf (12 x %.1f)", tbl_height))
} else {
  message("  SKIP: no exact or near matches")
}


# ==============================================================================
# SUPPLEMENTARY FIGURES
# ==============================================================================

# ------------------------------------------------------------------------------
# S1: GLIPH Overview Bubble (enriched clusters jitter)
# ------------------------------------------------------------------------------
message("--- Figure S1: GLIPH overview bubble ---")

if (exists("diagnosis_enrichment") && nrow(diagnosis_enrichment) > 0) {

  overview_data <- diagnosis_enrichment |>
    dplyr::filter(p.adj < 0.1, odds_ratio > 1) |>
    dplyr::mutate(
      neg_log10_fdr = -log10(pmax(p.adj, 1e-300)),
      diagnosis = factor(diagnosis, levels = intersect(c(neuropathy_dx, "GBS_Sukenikova"),
                                                        unique(diagnosis)))
    )

  if (nrow(overview_data) > 0) {
    # Add tissue bias if available
    if (exists("tissue_enrichment") && nrow(tissue_enrichment) > 0) {
      overview_data <- overview_data |>
        dplyr::left_join(
          tissue_enrichment |> dplyr::select(cluster, tissue_bias),
          by = "cluster"
        ) |>
        dplyr::mutate(tissue_bias = tidyr::replace_na(tissue_bias, "No bias"))
    } else {
      overview_data$tissue_bias <- "No bias"
    }

    tissue_shapes <- c("CSF-enriched" = 24, "No bias" = 21, "PBMC-enriched" = 25)

    # Top 2 labels per dx by cluster_size
    top_labels_s1 <- overview_data |>
      dplyr::group_by(diagnosis) |>
      dplyr::slice_max(order_by = cluster_size, n = 2, with_ties = FALSE) |>
      dplyr::ungroup()

    p_s1 <- ggplot(overview_data,
                   aes(x = diagnosis, y = cluster_size)) +
      geom_jitter(aes(fill = diagnosis, size = neg_log10_fdr, shape = tissue_bias),
                  width = 0.2, alpha = 0.75, stroke = 0.4) +
      scale_fill_manual(values = diagnosis_col, guide = "none") +
      scale_shape_manual(values = tissue_shapes, name = "Tissue compartment") +
      scale_size_continuous(range = c(2, 8), name = "-log10(FDR)") +
      geom_text_repel(
        data = top_labels_s1,
        aes(label = paste0("C", cluster)),
        size = 2.5, max.overlaps = 15, fontface = "italic",
        segment.color = "grey50", segment.size = 0.3, seed = 42
      ) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Disease-enriched GLIPH2 specificity groups overview",
        x = "",
        y = "Cluster size (# clonotypes)"
      )

    safe_save("fig_gliph_enriched_clusters_overview_bubble.pdf", p_s1, 9, 6.5)
  }
} else {
  message("  SKIP: no enrichment data for S1")
}


# ------------------------------------------------------------------------------
# S2: Diagnosis Sharing Heatmap (ComplexHeatmap)
# ------------------------------------------------------------------------------
message("--- Figure S2: Diagnosis sharing heatmap ---")

if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    exists("cluster_meta") && nrow(cluster_meta) > 0) {

  # Count shared clusters between each pair of diagnoses
  dx_pairs <- cluster_meta |>
    dplyr::filter(!is.na(diagnosis), !is.na(cluster)) |>
    dplyr::distinct(cluster, diagnosis)

  dx_list <- unique(dx_pairs$diagnosis)

  # Build sharing matrix
  share_mat <- matrix(0, nrow = length(dx_list), ncol = length(dx_list),
                      dimnames = list(dx_list, dx_list))

  for (i in seq_along(dx_list)) {
    for (j in seq_along(dx_list)) {
      cl_i <- dx_pairs$cluster[dx_pairs$diagnosis == dx_list[i]]
      cl_j <- dx_pairs$cluster[dx_pairs$diagnosis == dx_list[j]]
      share_mat[i, j] <- length(intersect(cl_i, cl_j))
    }
  }

  # Order by dx_order
  order_idx <- intersect(c(dx_order, "GBS_Sukenikova"), rownames(share_mat))
  share_mat <- share_mat[order_idx, order_idx]

  col_fun_s2 <- circlize::colorRamp2(
    c(0, max(share_mat[upper.tri(share_mat)]) / 2, max(share_mat[upper.tri(share_mat)])),
    c("white", "#FEE08B", "#D73027")
  )

  pdf(file.path(fig_dir, "fig_gliph_diagnosis_cluster_sharing_heatmap.pdf"), width = 8, height = 7)
  ht_s2 <- ComplexHeatmap::Heatmap(
    share_mat,
    name = "Shared\nclusters",
    col = col_fun_s2,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = "Diagnosis",
    column_title = "GLIPH2 cluster sharing between diagnoses",
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid::grid.text(share_mat[i, j], x, y, gp = grid::gpar(fontsize = 8))
    },
    border = TRUE
  )
  ComplexHeatmap::draw(ht_s2)
  dev.off()
  message("  Saved: fig_gliph_diagnosis_cluster_sharing_heatmap.pdf (8 x 7)")
}


# ------------------------------------------------------------------------------
# S3: Tissue Bias Fractions (fraction-based stacked bar)
# ------------------------------------------------------------------------------
message("--- Figure S3: Tissue bias fractions ---")

if (exists("tissue_bias_summary") && nrow(tissue_bias_summary) > 0) {

  tbs_long <- tissue_bias_summary |>
    tidyr::pivot_longer(
      cols = -diagnosis,
      names_to = "tissue_bias",
      values_to = "n_clusters"
    ) |>
    dplyr::group_by(diagnosis) |>
    dplyr::mutate(fraction = n_clusters / sum(n_clusters)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      tissue_bias = factor(tissue_bias, levels = c("CSF-enriched", "No bias", "PBMC-enriched")),
      diagnosis = factor(diagnosis, levels = intersect(c(dx_order, "GBS_Sukenikova"),
                                                        unique(diagnosis)))
    )

  tb_colors <- c("CSF-enriched" = "#E41A1C", "No bias" = "grey70", "PBMC-enriched" = "#377EB8")

  p_s3 <- ggplot(tbs_long, aes(x = diagnosis, y = fraction, fill = tissue_bias)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = tb_colors, name = "Tissue bias") +
    scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.02))) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Tissue compartmentalization of disease-enriched clusters",
      x = "Diagnosis",
      y = "Fraction of enriched clusters"
    )

  safe_save("fig_gliph_tissue_bias_fraction_by_diagnosis.pdf", p_s3, 8, 5.5)
} else {
  message("  SKIP: no tissue bias summary data")
}


# ------------------------------------------------------------------------------
# S4: Tissue Motifs Mirrored Bar (CSF vs PBMC)
# ------------------------------------------------------------------------------
message("--- Figure S4: Tissue motifs mirrored ---")

if (exists("tissue_motif_comparison") && !is.null(tissue_motif_comparison) &&
    nrow(tissue_motif_comparison) > 0) {

  # Top 10 most significant motifs
  top_tissue_motifs <- tissue_motif_comparison |>
    dplyr::filter(!is.na(log2_fc), !is.na(p.adj)) |>
    dplyr::arrange(p.adj) |>
    dplyr::slice_head(n = 20) |>
    dplyr::mutate(
      motif_label = paste0(motif, " (", k, "-mer)"),
      direction = ifelse(log2_fc > 0, "CSF", "PBMC"),
      sig_label = sapply(p.adj, format_pval)
    ) |>
    dplyr::arrange(log2_fc) |>
    dplyr::mutate(motif_label = forcats::fct_inorder(motif_label))

  if (nrow(top_tissue_motifs) > 0) {
    p_s4 <- ggplot(top_tissue_motifs,
                   aes(x = log2_fc, y = motif_label, fill = direction)) +
      geom_col(width = 0.7, show.legend = FALSE) +
      geom_text(aes(label = sig_label),
                hjust = ifelse(top_tissue_motifs$log2_fc > 0, -0.3, 1.3),
                size = 2.8) +
      geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
      scale_fill_manual(values = c("CSF" = "#E41A1C", "PBMC" = "#377EB8")) +
      annotate("text", x = min(top_tissue_motifs$log2_fc) * 0.9, y = Inf,
               label = "PBMC-enriched", color = "#377EB8", fontface = "bold",
               hjust = 0, vjust = 1.5, size = 3.5) +
      annotate("text", x = max(top_tissue_motifs$log2_fc) * 0.9, y = Inf,
               label = "CSF-enriched", color = "#E41A1C", fontface = "bold",
               hjust = 1, vjust = 1.5, size = 3.5) +
      theme_pub() +
      labs(
        title = "CDR3 motifs: CSF-enriched vs PBMC-enriched clusters",
        x = expression(log[2](CSF / PBMC ~ frequency)),
        y = ""
      )

    safe_save("fig_gliph_csf_vs_pbmc_motif_enrichment.pdf", p_s4, 8, 6)
  }
} else {
  message("  SKIP: no tissue motif comparison data")
}


# ------------------------------------------------------------------------------
# S5: Myelin Enrichment Bar (per diagnosis)
# ------------------------------------------------------------------------------
message("--- Figure S5: Myelin enrichment bar ---")

if (exists("myelin_enrichment") && nrow(myelin_enrichment) > 0) {

  me_plot <- myelin_enrichment |>
    dplyr::mutate(
      log2_OR = log2(pmax(odds_ratio, 1e-4)),
      sig_label = sapply(p.adj, format_pval),
      diagnosis = factor(diagnosis, levels = intersect(neuropathy_dx, diagnosis))
    ) |>
    dplyr::filter(!is.na(diagnosis))

  if (nrow(me_plot) > 0) {
    p_s5 <- ggplot(me_plot, aes(x = diagnosis, y = fraction * 100, fill = diagnosis)) +
      geom_col(width = 0.7, show.legend = FALSE) +
      geom_text(aes(label = paste0(sig_label, "\n", sprintf("OR=%.1f", odds_ratio))),
                vjust = -0.3, size = 2.8) +
      scale_fill_manual(values = diagnosis_col) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Myelin-reactive TCR fraction by diagnosis",
        subtitle = "Exact + cluster-level matches to Sukenikova myelin-reactive clones",
        x = "",
        y = "Myelin-reactive TCRs (%)"
      ) +
      expand_limits(y = max(me_plot$fraction * 100, na.rm = TRUE) * 1.4)

    safe_save("fig_gliph_myelin_reactive_fraction_by_diagnosis.pdf", p_s5, 7, 5.5)
  }
} else {
  message("  SKIP: no myelin enrichment data")
}


# ------------------------------------------------------------------------------
# S6: Cross-Cohort Overlap Bar (shared clusters per Heming dx)
# ------------------------------------------------------------------------------
message("--- Figure S6: Cross-cohort overlap ---")

if (exists("cross_dx") && nrow(cross_dx) > 0) {

  cc_bar <- cross_dx |>
    dplyr::filter(heming_diagnoses != "CTRL") |>
    dplyr::mutate(
      heming_diagnoses = factor(heming_diagnoses,
                                levels = intersect(c(neuropathy_dx, "GBS_Sukenikova"),
                                                    heming_diagnoses))
    )

  if (nrow(cc_bar) > 0) {
    p_s6 <- ggplot(cc_bar, aes(x = heming_diagnoses, y = n_shared_clusters,
                                fill = heming_diagnoses)) +
      geom_col(width = 0.7, show.legend = FALSE) +
      geom_text(aes(label = n_shared_clusters), vjust = -0.5, size = 3.5) +
      scale_fill_manual(values = diagnosis_col) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "GLIPH2 clusters shared with Sukenikova GBS cohort",
        x = "Heming diagnosis",
        y = "Number of shared clusters"
      ) +
      expand_limits(y = max(cc_bar$n_shared_clusters, na.rm = TRUE) * 1.15)

    safe_save("fig_gliph_sukenikova_shared_clusters_by_diagnosis.pdf", p_s6, 7, 5)
  }
} else {
  message("  SKIP: no cross_dx data")
}


# ==============================================================================
# Done
# ==============================================================================
message("\n============================================================")
message("All manuscript figures generated!")
message(sprintf("  Output directory: %s", fig_dir))
message(sprintf("  Total PDF files: %d",
                length(list.files(fig_dir, pattern = "\\.pdf$"))))
message("============================================================")
