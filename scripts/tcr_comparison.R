# TCR repertoire comparison across neuropathy diagnoses
# Compares TRA, TRB, TRAB chains with mixed models (patient random effect, tissue fixed effect)

# Deactivate renv to use system/CRAN packages
if (nzchar(Sys.getenv("RENV_PROJECT"))) {
  Sys.setenv(RENV_PROJECT = "")
  message("Note: renv deactivated, using system packages")
}

set.seed(42)

# Libraries
library(qs)
library(Seurat)
library(tidyverse)
library(scRepertoire)
library(writexl)
library(readxl)
library(patchwork)
library(ggpubr)
library(ggsignif)
library(Peptides)
library(viridis)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(immApex)
library(vegan)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(emmeans)
library(ggrepel)
library(cluster)
library(umap)
library(uwot)

devtools::load_all("/Users/nick/Documents/GitHub/immLynx")

# Constants
chains <- c("TRA", "TRB", "TRAB")
aa_codes <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
neuropathy_dx <- c("CIAP", "CIDP", "GBS", "MAG", "MFS", "PNC", "CAN", "PPN")

# Output directories
output_dirs <- c(
  "results/tcr_comparison",
  "results/tcr_comparison/gene_usage",
  "results/tcr_comparison/cdr3_features",
  "results/tcr_comparison/physicochemical",
  "results/tcr_comparison/clonality",
  "results/tcr_comparison/similarity",
  "results/tcr_comparison/embeddings",
  "results/tcr_comparison/motifs",
  "results/tcr_comparison/tcrdist",
  "results/tcr_comparison/gliph",
  "results/tcr_comparison/tables"
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Load data
sc_tcr <- qs::qread("sc_tcr.qs", nthreads = 6)

# Color palettes
diagnosis_col <- sc_tcr@misc$diagnosis_col
tissue_diagnosis_col <- sc_tcr@misc$tissue_diagnosis_col
cluster_col <- sc_tcr@misc$cluster_col
tissue_col <- c("CSF" = "#E41A1C", "PBMC" = "#377EB8")
group_col <- c("CTRL" = "#4DAF4A", "PNP" = "#984EA3")

# Extract TCR metadata
tcr_metadata <- sc_tcr@meta.data |>
  tibble::rownames_to_column("cell_id") |>
  dplyr::filter(!is.na(CTaa)) |>
  dplyr::mutate(
    TRA_CDR3 = gsub("_.*", "", CTaa),
    TRB_CDR3 = gsub(".*_", "", CTaa),
    TRA_length = nchar(TRA_CDR3),
    TRB_length = nchar(TRB_CDR3)
  )

# Unique clonotypes per sample
tcr_clonotypes <- tcr_metadata |>
  dplyr::group_by(CTaa, sample, patient, tissue, group, diagnosis, tissue_group, tissue_diagnosis) |>
  dplyr::summarize(
    n_cells = n(),
    TRA_CDR3 = first(TRA_CDR3),
    TRB_CDR3 = first(TRB_CDR3),
    TRA_length = first(TRA_length),
    TRB_length = first(TRB_length),
    .groups = "drop"
  )

# Sample-level summary for aggregated analyses
sample_summary <- tcr_metadata |>
  dplyr::distinct(sample, patient, tissue, group, diagnosis)

# ============================================================================
# Helper functions
# ============================================================================

calc_aa_composition <- function(seq) {
  if (is.na(seq) || nchar(seq) == 0) return(setNames(rep(NA_real_, 20), aa_codes))
  aa <- strsplit(seq, "")[[1]]
  props <- as.numeric(table(factor(aa, levels = aa_codes))) / length(aa)
  setNames(props, aa_codes)
}

calc_trab_composition <- function(tra, trb) {
  combined <- paste0(ifelse(is.na(tra), "", tra), ifelse(is.na(trb), "", trb))
  if (nchar(combined) == 0) return(setNames(rep(NA_real_, 20), aa_codes))
  calc_aa_composition(combined)
}

# Statistical testing: diagnosis vs CTRL with patient random effect
run_diagnosis_stats <- function(data, value_col, include_tissue = TRUE) {
  data_clean <- data |>
    dplyr::filter(!is.na(.data[[value_col]]), !is.na(diagnosis), !is.na(patient)) |>
    dplyr::mutate(diagnosis = factor(diagnosis, levels = c("CTRL", neuropathy_dx)))

  if (nrow(data_clean) < 20 || length(unique(data_clean$patient)) < 5) {
    return(data.frame(
      diagnosis = neuropathy_dx, estimate = NA, SE = NA, t.ratio = NA,
      p.value = NA, p.adj = NA, method = "insufficient_data"
    ))
  }

  has_tissues <- length(unique(data_clean$tissue)) == 2

  tryCatch({
    if (include_tissue && has_tissues) {
      formula_str <- paste0(value_col, " ~ diagnosis * tissue + (1|patient)")
    } else {
      formula_str <- paste0(value_col, " ~ diagnosis + (1|patient)")
    }

    model <- lmerTest::lmer(as.formula(formula_str), data = data_clean)
    emm <- emmeans::emmeans(model, "diagnosis")
    contrasts <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1)

    result <- as.data.frame(contrasts) |>
      dplyr::mutate(
        diagnosis = gsub(" - CTRL", "", contrast),
        method = "emmeans_dunnett"
      ) |>
      dplyr::select(diagnosis, estimate, SE, t.ratio, p.value, method)

    means <- data_clean |>
      dplyr::group_by(diagnosis) |>
      dplyr::summarize(mean = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

    ctrl_mean <- means$mean[means$diagnosis == "CTRL"]
    result <- result |>
      dplyr::left_join(means, by = "diagnosis") |>
      dplyr::rename(mean_dx = mean) |>
      dplyr::mutate(mean_ctrl = ctrl_mean)

    result

  }, error = function(e) {
    ctrl_data <- data_clean |> dplyr::filter(diagnosis == "CTRL")

    results <- lapply(neuropathy_dx, function(dx) {
      dx_data <- data_clean |> dplyr::filter(diagnosis == dx)
      if (nrow(dx_data) < 3) {
        return(data.frame(diagnosis = dx, estimate = NA, SE = NA, t.ratio = NA,
                          p.value = NA, mean_dx = NA, mean_ctrl = NA, method = "insufficient"))
      }

      test <- tryCatch(
        wilcox.test(ctrl_data[[value_col]], dx_data[[value_col]]),
        error = function(e) list(p.value = NA)
      )

      data.frame(
        diagnosis = dx,
        estimate = mean(dx_data[[value_col]], na.rm = TRUE) - mean(ctrl_data[[value_col]], na.rm = TRUE),
        SE = NA, t.ratio = NA,
        p.value = test$p.value,
        mean_dx = mean(dx_data[[value_col]], na.rm = TRUE),
        mean_ctrl = mean(ctrl_data[[value_col]], na.rm = TRUE),
        method = "wilcoxon"
      )
    }) |> dplyr::bind_rows()

    results
  })
}

# PERMANOVA with patient stratification (macOS safe)
run_permanova <- function(feature_matrix, metadata, formula_str, strata_var = "patient", n_perm = 999) {
  feature_matrix <- as.matrix(feature_matrix)
  storage.mode(feature_matrix) <- "numeric"

  # Filter to complete cases
  complete_rows <- complete.cases(feature_matrix)
  feature_matrix <- feature_matrix[complete_rows, , drop = FALSE]
  metadata <- metadata[complete_rows, , drop = FALSE]

  # Remove zero variance columns
  col_vars <- apply(feature_matrix, 2, var, na.rm = TRUE)
  feature_matrix <- feature_matrix[, col_vars > 1e-10, drop = FALSE]

  # Check dimensions
  if (nrow(feature_matrix) < 10 || ncol(feature_matrix) < 2) {
    return(list(success = FALSE, message = paste("Too few samples:", nrow(feature_matrix), "or features:", ncol(feature_matrix))))
  }

  # Ensure metadata factors are valid
  metadata <- as.data.frame(metadata)
  for (var in all.vars(as.formula(formula_str))) {
    if (var %in% colnames(metadata)) {
      metadata[[var]] <- droplevels(as.factor(metadata[[var]]))
    }
  }

  use_parallel <- if (Sys.info()["sysname"] == "Darwin") NULL else max(1, parallel::detectCores() - 1)

  tryCatch({
    # Create distance matrix first
    dist_mat <- vegan::vegdist(feature_matrix, method = "euclidean")

    result <- vegan::adonis2(
      dist_mat ~ diagnosis * tissue,
      data = metadata,
      permutations = n_perm,
      strata = metadata[[strata_var]],
      parallel = use_parallel
    )
    list(success = TRUE, result = result)
  }, error = function(e) {
    list(success = FALSE, message = e$message)
  })
}

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

boxplot_with_signif <- function(data, y_var, x_var = "diagnosis", stats_df = NULL,
                                 colors = diagnosis_col, title = "", y_lab = NULL,
                                 p_threshold = 0.05) {
  y_lab <- y_lab %||% y_var

  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 0.8) +
    scale_fill_manual(values = colors) +
    theme_pub() +
    labs(title = title, x = "", y = y_lab) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.null(stats_df)) {
    sig_dx <- stats_df |> dplyr::filter(p.adj < p_threshold)
    if (nrow(sig_dx) > 0) {
      comparisons <- lapply(sig_dx$diagnosis, function(dx) c("CTRL", dx))
      annotations <- sapply(sig_dx$p.adj, function(p) {
        if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
      })
      p <- p + ggsignif::geom_signif(
        comparisons = comparisons, annotations = annotations,
        step_increase = 0.08, textsize = 3, vjust = 0.5
      )
    }
  }
  p
}

format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

# ============================================================================
# Section 1: CDR3 amino acid composition
# ============================================================================
message("Section 1: CDR3 amino acid composition")

# Calculate AA composition for each chain
tcr_aa_data <- tcr_clonotypes |>
  dplyr::mutate(
    TRA_aa = purrr::map(TRA_CDR3, calc_aa_composition),
    TRB_aa = purrr::map(TRB_CDR3, calc_aa_composition),
    TRAB_aa = purrr::map2(TRA_CDR3, TRB_CDR3, calc_trab_composition)
  )

tcr_tra_aa <- tcr_aa_data |>
  dplyr::select(-TRB_aa, -TRAB_aa) |>
  tidyr::unnest_wider(TRA_aa, names_sep = "_") |>
  dplyr::rename_with(~ gsub("TRA_aa_", "TRA_", .x), starts_with("TRA_aa_"))

tcr_trb_aa <- tcr_aa_data |>
  dplyr::select(-TRA_aa, -TRAB_aa) |>
  tidyr::unnest_wider(TRB_aa, names_sep = "_") |>
  dplyr::rename_with(~ gsub("TRB_aa_", "TRB_", .x), starts_with("TRB_aa_"))

tcr_trab_aa <- tcr_aa_data |>
  dplyr::select(-TRA_aa, -TRB_aa) |>
  tidyr::unnest_wider(TRAB_aa, names_sep = "_") |>
  dplyr::rename_with(~ gsub("TRAB_aa_", "TRAB_", .x), starts_with("TRAB_aa_"))

# Run stats for each chain
aa_stats <- list()
for (chain in chains) {
  message(paste("  Running stats for", chain))
  data_use <- switch(chain, TRA = tcr_tra_aa, TRB = tcr_trb_aa, TRAB = tcr_trab_aa)

  chain_results <- lapply(aa_codes, function(aa) {
    col_name <- paste0(chain, "_", aa)
    if (!col_name %in% colnames(data_use)) return(NULL)

    result <- run_diagnosis_stats(data_use, col_name)
    result$amino_acid <- aa
    result$chain <- chain
    result
  }) |> dplyr::bind_rows()

  chain_results <- chain_results |>
    dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

  aa_stats[[chain]] <- chain_results
}

aa_stats_combined <- dplyr::bind_rows(aa_stats)

# AA composition heatmap showing mean difference from CTRL
message("  Creating AA composition heatmap")
aa_heatmap_data <- aa_stats_combined |>
  dplyr::filter(!is.na(estimate)) |>
  dplyr::select(chain, amino_acid, diagnosis, estimate, p.adj) |>
  dplyr::mutate(
    sig = ifelse(p.adj < 0.05, "*", ""),
    label = paste0(round(estimate, 3), sig)
  )

for (ch in chains) {
  ch_data <- aa_heatmap_data |> dplyr::filter(chain == ch)
  if (nrow(ch_data) == 0) next

  mat <- ch_data |>
    dplyr::select(amino_acid, diagnosis, estimate) |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = estimate) |>
    tibble::column_to_rownames("amino_acid") |>
    as.matrix()

  sig_mat <- ch_data |>
    dplyr::select(amino_acid, diagnosis, sig) |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = sig) |>
    tibble::column_to_rownames("amino_acid") |>
    as.matrix()

  pdf(paste0("results/tcr_comparison/cdr3_features/aa_composition_heatmap_", ch, ".pdf"),
      width = 10, height = 8)
  pheatmap::pheatmap(
    mat,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-max(abs(mat), na.rm = TRUE), max(abs(mat), na.rm = TRUE), length.out = 101),
    main = paste(ch, "AA composition difference from CTRL"),
    display_numbers = sig_mat,
    cluster_cols = FALSE,
    fontsize_number = 12
  )
  dev.off()
}

# PERMANOVA for multivariate composition
message("  Running PERMANOVA")
permanova_results <- list()
for (chain in chains) {
  data_use <- switch(chain, TRA = tcr_tra_aa, TRB = tcr_trb_aa, TRAB = tcr_trab_aa)
  aa_col_names <- paste0(chain, "_", aa_codes)
  aa_col_names <- aa_col_names[aa_col_names %in% colnames(data_use)]

  # Aggregate to sample level for PERMANOVA
  sample_aa <- data_use |>
    dplyr::group_by(sample, patient, tissue, diagnosis) |>
    dplyr::summarize(
      dplyr::across(all_of(aa_col_names), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::filter(complete.cases(dplyr::across(all_of(aa_col_names))))

  aa_matrix <- sample_aa |> dplyr::select(all_of(aa_col_names)) |> as.matrix()

  message(paste("   ", chain, "- samples:", nrow(aa_matrix), ", features:", ncol(aa_matrix)))

  perm_result <- run_permanova(aa_matrix, sample_aa, "~ diagnosis * tissue", "patient")
  if (perm_result$success) {
    permanova_results[[chain]] <- perm_result$result
    diag_p <- perm_result$result["diagnosis", "Pr(>F)"]
    diag_r2 <- perm_result$result["diagnosis", "R2"]
    tissue_p <- perm_result$result["tissue", "Pr(>F)"]
    tissue_r2 <- perm_result$result["tissue", "R2"]
    message(paste("     diagnosis: p =", round(diag_p, 4), ", R2 =", round(diag_r2, 4)))
    message(paste("     tissue: p =", round(tissue_p, 4), ", R2 =", round(tissue_r2, 4)))
  } else {
    message(paste("     PERMANOVA failed:", perm_result$message))
  }
}

# PCA/UMAP for AA composition
message("  Creating PCA and UMAP visualizations")
for (chain in chains) {
  data_use <- switch(chain, TRA = tcr_tra_aa, TRB = tcr_trb_aa, TRAB = tcr_trab_aa)
  aa_col_names <- paste0(chain, "_", aa_codes)
  aa_col_names <- aa_col_names[aa_col_names %in% colnames(data_use)]

  # Sample-level aggregation
  sample_aa <- data_use |>
    dplyr::group_by(sample, patient, tissue, diagnosis, group) |>
    dplyr::summarize(
      dplyr::across(all_of(aa_col_names), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::filter(complete.cases(dplyr::across(all_of(aa_col_names))))

  aa_mat <- sample_aa |> dplyr::select(all_of(aa_col_names)) |> as.matrix()

  # PCA
  pca <- prcomp(aa_mat, scale. = TRUE, center = TRUE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100

  pca_df <- data.frame(
    PC1 = pca$x[,1], PC2 = pca$x[,2],
    diagnosis = sample_aa$diagnosis, tissue = sample_aa$tissue
  )

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diagnosis, shape = tissue)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = diagnosis), level = 0.95, linetype = "dashed") +
    scale_color_manual(values = diagnosis_col) +
    theme_pub() +
    labs(
      title = paste(chain, "AA composition PCA"),
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    )

  ggsave(paste0("results/tcr_comparison/cdr3_features/aa_pca_", chain, ".pdf"), p_pca, width = 9, height = 7)

  # UMAP
  set.seed(42)
  umap_result <- uwot::umap(aa_mat, n_neighbors = min(15, nrow(aa_mat) - 1), min_dist = 0.1)

  umap_df <- data.frame(
    UMAP1 = umap_result[,1], UMAP2 = umap_result[,2],
    diagnosis = sample_aa$diagnosis, tissue = sample_aa$tissue
  )

  p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = diagnosis, shape = tissue)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = diagnosis_col) +
    theme_pub() +
    labs(title = paste(chain, "AA composition UMAP"))

  ggsave(paste0("results/tcr_comparison/cdr3_features/aa_umap_", chain, ".pdf"), p_umap, width = 9, height = 7)
}

# Top differential AA boxplots with stats
message("  Creating boxplots for top differential AAs")
top_aa <- aa_stats_combined |>
  dplyr::filter(!is.na(p.adj)) |>
  dplyr::group_by(chain, amino_acid) |>
  dplyr::summarize(min_p = min(p.adj, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(min_p) |>
  dplyr::slice(1:9)

aa_boxplots <- lapply(seq_len(nrow(top_aa)), function(i) {
  ch <- top_aa$chain[i]
  aa <- top_aa$amino_acid[i]
  col_name <- paste0(ch, "_", aa)

  data_use <- switch(ch, TRA = tcr_tra_aa, TRB = tcr_trb_aa, TRAB = tcr_trab_aa)
  stats_use <- aa_stats[[ch]] |> dplyr::filter(amino_acid == aa)

  boxplot_with_signif(data_use, y_var = col_name, stats_df = stats_use,
                      colors = diagnosis_col, title = paste(ch, aa))
})

aa_boxplots_combined <- wrap_plots(aa_boxplots, ncol = 3)
ggsave("results/tcr_comparison/cdr3_features/aa_top_differential_boxplots.pdf",
       aa_boxplots_combined, width = 14, height = 12)

writexl::write_xlsx(
  list(TRA = aa_stats$TRA, TRB = aa_stats$TRB, TRAB = aa_stats$TRAB, combined = aa_stats_combined),
  "results/tcr_comparison/cdr3_features/aa_composition_statistics.xlsx"
)

# ============================================================================
# Section 2: Positional motif analysis
# ============================================================================
message("Section 2: Positional motif analysis")

motif_results <- list()

for (chain in c("TRA", "TRB")) {
  seq_col <- paste0(chain, "_CDR3")

  for (dx in c("CTRL", neuropathy_dx)) {
    seqs <- tcr_clonotypes |>
      dplyr::filter(diagnosis == dx) |>
      dplyr::pull(!!sym(seq_col)) |>
      na.omit()

    if (length(seqs) < 100) next

    motifs <- immApex::calculateMotif(
      seqs, motif.lengths = 3:4, min.depth = 5, discontinuous = FALSE
    )

    motifs$diagnosis <- dx
    motifs$chain <- chain
    motifs$prop <- motifs$frequency / length(seqs)

    motif_results[[paste(chain, dx, sep = "_")]] <- motifs
  }
}

motif_combined <- dplyr::bind_rows(motif_results)

# Compare motif enrichment: each diagnosis vs CTRL
message("  Comparing motif enrichment vs CTRL")
motif_enrichment <- list()

for (chain in c("TRA", "TRB")) {
  ctrl_motifs <- motif_combined |>
    dplyr::filter(chain == !!chain, diagnosis == "CTRL") |>
    dplyr::select(motif, prop_ctrl = prop)

  for (dx in neuropathy_dx) {
    dx_motifs <- motif_combined |>
      dplyr::filter(chain == !!chain, diagnosis == dx) |>
      dplyr::select(motif, prop_dx = prop, freq_dx = frequency)

    if (nrow(dx_motifs) == 0) next

    comparison <- dplyr::full_join(ctrl_motifs, dx_motifs, by = "motif") |>
      dplyr::mutate(
        prop_ctrl = tidyr::replace_na(prop_ctrl, 0),
        prop_dx = tidyr::replace_na(prop_dx, 0),
        log2FC = log2((prop_dx + 0.0001) / (prop_ctrl + 0.0001)),
        diagnosis = dx,
        chain = chain
      ) |>
      dplyr::filter(!is.na(freq_dx), freq_dx >= 10)

    motif_enrichment[[paste(chain, dx, sep = "_")]] <- comparison
  }
}

motif_enrichment_combined <- dplyr::bind_rows(motif_enrichment)

# Top enriched/depleted motifs per diagnosis
top_motifs_per_dx <- motif_enrichment_combined |>
  dplyr::group_by(chain, diagnosis) |>
  dplyr::arrange(desc(abs(log2FC))) |>
  dplyr::slice(1:10) |>
  dplyr::ungroup()

# Motif heatmap
for (ch in c("TRA", "TRB")) {
  ch_motifs <- top_motifs_per_dx |> dplyr::filter(chain == ch)
  if (nrow(ch_motifs) < 5) next

  top_motif_list <- unique(ch_motifs$motif)

  heatmap_data <- motif_enrichment_combined |>
    dplyr::filter(chain == ch, motif %in% top_motif_list) |>
    dplyr::select(motif, diagnosis, log2FC) |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = log2FC, values_fill = 0) |>
    tibble::column_to_rownames("motif") |>
    as.matrix()

  if (nrow(heatmap_data) < 3) next

  pdf(paste0("results/tcr_comparison/motifs/motif_enrichment_heatmap_", ch, ".pdf"),
      width = 10, height = 8)
  pheatmap::pheatmap(
    heatmap_data,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = paste(ch, "motif enrichment (log2FC vs CTRL)"),
    cluster_cols = FALSE
  )
  dev.off()
}

writexl::write_xlsx(
  list(
    all_motifs = motif_combined,
    enrichment = motif_enrichment_combined,
    top_motifs = top_motifs_per_dx
  ),
  "results/tcr_comparison/motifs/motif_analysis.xlsx"
)

# ============================================================================
# Section 3: V/J gene usage
# ============================================================================
message("Section 3: V/J gene usage")

gene_usage_stats <- list()
gene_usage_props <- list()

for (chain in c("TRA", "TRB")) {
  for (gene_type in c("V", "J")) {
    message(paste("  Analyzing", chain, gene_type, "genes"))

    gene_data <- tcr_metadata |>
      dplyr::mutate(
        gene = case_when(
          chain == "TRA" & gene_type == "V" ~ gsub("\\..*", "", gsub("_.*", "", CTgene)),
          chain == "TRA" & gene_type == "J" ~ gsub(".*\\.", "", gsub("_.*", "", CTgene)),
          chain == "TRB" & gene_type == "V" ~ gsub("\\..*", "", gsub(".*_", "", CTgene)),
          chain == "TRB" & gene_type == "J" ~ gsub(".*\\.", "", gsub(".*_", "", CTgene))
        )
      ) |>
      dplyr::filter(!is.na(gene), gene != "", nchar(gene) > 2)

    top_genes <- gene_data |>
      dplyr::count(gene, sort = TRUE) |>
      dplyr::slice(1:25) |>
      dplyr::pull(gene)

    sample_gene_props <- gene_data |>
      dplyr::filter(gene %in% top_genes) |>
      dplyr::count(sample, patient, tissue, diagnosis, group, gene) |>
      dplyr::group_by(sample) |>
      dplyr::mutate(prop = n / sum(n)) |>
      dplyr::ungroup()

    gene_usage_props[[paste(chain, gene_type, sep = "_")]] <- sample_gene_props

    gene_results <- lapply(top_genes, function(g) {
      gene_df <- sample_gene_props |>
        dplyr::filter(gene == g) |>
        tidyr::complete(
          tidyr::nesting(sample, patient, tissue, diagnosis, group),
          fill = list(n = 0, prop = 0)
        )

      if (nrow(gene_df) < 10) return(NULL)

      result <- run_diagnosis_stats(gene_df, "prop")
      result$gene <- g
      result$chain <- chain
      result$gene_type <- gene_type
      result
    }) |> dplyr::bind_rows()

    if (nrow(gene_results) > 0) {
      gene_results <- gene_results |>
        dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))
      gene_usage_stats[[paste(chain, gene_type, sep = "_")]] <- gene_results
    }
  }
}

gene_usage_combined <- dplyr::bind_rows(gene_usage_stats, .id = "comparison")

# Gene usage PCA per chain/gene_type
message("  Creating gene usage PCA")
for (key in names(gene_usage_props)) {
  props_data <- gene_usage_props[[key]]

  wide_data <- props_data |>
    tidyr::pivot_wider(
      id_cols = c(sample, patient, tissue, diagnosis, group),
      names_from = gene, values_from = prop, values_fill = 0
    )

  gene_cols <- setdiff(colnames(wide_data), c("sample", "patient", "tissue", "diagnosis", "group"))
  if (length(gene_cols) < 3) next

  gene_mat <- wide_data |> dplyr::select(all_of(gene_cols)) |> as.matrix()

  pca <- prcomp(gene_mat, scale. = TRUE, center = TRUE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100

  pca_df <- data.frame(
    PC1 = pca$x[,1], PC2 = pca$x[,2],
    diagnosis = wide_data$diagnosis, tissue = wide_data$tissue
  )

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diagnosis, shape = tissue)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(group = diagnosis), level = 0.95, linetype = "dashed") +
    scale_color_manual(values = diagnosis_col) +
    theme_pub() +
    labs(
      title = paste(key, "gene usage PCA"),
      x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
      y = paste0("PC2 (", round(var_exp[2], 1), "%)")
    )

  ggsave(paste0("results/tcr_comparison/gene_usage/pca_", key, ".pdf"), p, width = 9, height = 7)
}

# Volcano plots for gene usage
message("  Creating volcano plots")
for (dx in neuropathy_dx) {
  dx_data <- gene_usage_combined |>
    dplyr::filter(diagnosis == dx, !is.na(estimate), !is.na(p.value))

  if (nrow(dx_data) < 5) next

  p <- ggplot(dx_data, aes(x = estimate, y = -log10(p.value))) +
    geom_point(aes(color = p.adj < 0.1), size = 2, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(
      data = dx_data |> dplyr::filter(p.adj < 0.1),
      aes(label = gene), size = 3, max.overlaps = 20
    ) +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red3")) +
    facet_wrap(~ comparison, scales = "free") +
    theme_pub() +
    labs(title = paste(dx, "vs CTRL: V/J gene usage"),
         x = "Difference in proportion", y = "-log10(p-value)") +
    theme(legend.position = "none")

  ggsave(paste0("results/tcr_comparison/gene_usage/volcano_", dx, ".pdf"), p, width = 10, height = 8)
}

# Gene usage heatmap
sig_genes <- gene_usage_combined |>
  dplyr::filter(p.adj < 0.1) |>
  dplyr::distinct(gene) |>
  dplyr::pull(gene)

if (length(sig_genes) > 3) {
  heatmap_data <- gene_usage_combined |>
    dplyr::filter(gene %in% sig_genes) |>
    dplyr::select(comparison, gene, diagnosis, estimate) |>
    dplyr::mutate(gene_comp = paste(comparison, gene, sep = "_")) |>
    dplyr::select(gene_comp, diagnosis, estimate) |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = estimate, values_fill = 0) |>
    tibble::column_to_rownames("gene_comp") |>
    as.matrix()

  pdf("results/tcr_comparison/gene_usage/significant_genes_heatmap.pdf", width = 12, height = 10)
  pheatmap::pheatmap(
    heatmap_data,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "Significant V/J genes (difference from CTRL)",
    cluster_cols = FALSE
  )
  dev.off()
}

writexl::write_xlsx(gene_usage_stats, "results/tcr_comparison/gene_usage/gene_usage_statistics.xlsx")

# ============================================================================
# Section 4: CDR3 length
# ============================================================================
message("Section 4: CDR3 length")

length_stats <- list()

for (chain in c("TRA", "TRB")) {
  length_col <- paste0(chain, "_length")

  result <- run_diagnosis_stats(tcr_clonotypes, length_col)
  result$chain <- chain
  result <- result |> dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

  length_stats[[chain]] <- result
}

length_stats_combined <- dplyr::bind_rows(length_stats)

# Boxplots with stats
length_plots <- lapply(c("TRA", "TRB"), function(ch) {
  boxplot_with_signif(
    tcr_clonotypes, y_var = paste0(ch, "_length"),
    stats_df = length_stats[[ch]],
    colors = diagnosis_col, title = paste(ch, "CDR3 length"),
    y_lab = "Length (amino acids)"
  )
})

length_combined <- wrap_plots(length_plots, ncol = 2)
ggsave("results/tcr_comparison/cdr3_features/cdr3_length_boxplots.pdf", length_combined, width = 12, height = 5)

# Length by tissue and diagnosis
p_length_facet <- tcr_clonotypes |>
  tidyr::pivot_longer(c(TRA_length, TRB_length), names_to = "chain", values_to = "length") |>
  dplyr::mutate(chain = gsub("_length", "", chain)) |>
  ggplot(aes(x = diagnosis, y = length, fill = diagnosis)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_grid(chain ~ tissue) +
  scale_fill_manual(values = diagnosis_col) +
  theme_pub() +
  labs(title = "CDR3 length by diagnosis and tissue", y = "Length (amino acids)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/tcr_comparison/cdr3_features/cdr3_length_by_tissue.pdf", p_length_facet, width = 12, height = 6)

writexl::write_xlsx(
  list(TRA = length_stats$TRA, TRB = length_stats$TRB),
  "results/tcr_comparison/cdr3_features/cdr3_length_stats.xlsx"
)

# ============================================================================
# Section 5: Physicochemical properties
# ============================================================================
message("Section 5: Physicochemical properties")

unique_seqs <- tcr_clonotypes |> dplyr::distinct(CTaa, TRA_CDR3, TRB_CDR3)

# Kidera factors
message("  Computing Kidera factors")
kidera_tra <- immApex::sequenceEncoder(unique_seqs$TRA_CDR3, mode = "property", property.set = "kideraFactors")[[2]]
rownames(kidera_tra) <- unique_seqs$CTaa
colnames(kidera_tra) <- paste0("TRA_KF", seq_len(ncol(kidera_tra)))

kidera_trb <- immApex::sequenceEncoder(unique_seqs$TRB_CDR3, mode = "property", property.set = "kideraFactors")[[2]]
rownames(kidera_trb) <- unique_seqs$CTaa
colnames(kidera_trb) <- paste0("TRB_KF", seq_len(ncol(kidera_trb)))

# Atchley factors
message("  Computing Atchley factors")
atchley_tra <- immApex::sequenceEncoder(unique_seqs$TRA_CDR3, mode = "property", property.set = "atchleyFactors")[[2]]
rownames(atchley_tra) <- unique_seqs$CTaa
colnames(atchley_tra) <- paste0("TRA_AF", seq_len(ncol(atchley_tra)))

atchley_trb <- immApex::sequenceEncoder(unique_seqs$TRB_CDR3, mode = "property", property.set = "atchleyFactors")[[2]]
rownames(atchley_trb) <- unique_seqs$CTaa
colnames(atchley_trb) <- paste0("TRB_AF", seq_len(ncol(atchley_trb)))

# Peptide properties
message("  Computing peptide properties")
compute_peptide_features <- function(seq) {
  if (is.na(seq) || nchar(seq) < 3) return(rep(NA, 8))
  tryCatch({
    c(
      mw = Peptides::mw(seq),
      charge = Peptides::charge(seq, pH = 7.4),
      pI = Peptides::pI(seq),
      aindex = Peptides::aIndex(seq),
      instaindex = Peptides::instaIndex(seq),
      boman = Peptides::boman(seq),
      hydrophobicity = Peptides::hydrophobicity(seq, scale = "KyteDoolittle"),
      hmoment = Peptides::hmoment(seq, angle = 100, window = 5)
    )
  }, error = function(e) rep(NA, 8))
}

peptide_tra <- t(sapply(unique_seqs$TRA_CDR3, compute_peptide_features))
colnames(peptide_tra) <- paste0("TRA_", c("mw", "charge", "pI", "aindex", "instaindex", "boman", "hydrophobicity", "hmoment"))
rownames(peptide_tra) <- unique_seqs$CTaa

peptide_trb <- t(sapply(unique_seqs$TRB_CDR3, compute_peptide_features))
colnames(peptide_trb) <- paste0("TRB_", c("mw", "charge", "pI", "aindex", "instaindex", "boman", "hydrophobicity", "hmoment"))
rownames(peptide_trb) <- unique_seqs$CTaa

physico_features <- cbind(
  kidera_tra, kidera_trb,
  atchley_tra, atchley_trb,
  peptide_tra, peptide_trb
) |>
  as.data.frame() |>
  tibble::rownames_to_column("CTaa")

physico_data <- tcr_clonotypes |>
  dplyr::left_join(physico_features, by = "CTaa")

# Stats for each feature
feature_cols <- colnames(physico_features)[-1]
message("  Running statistics for ", length(feature_cols), " features")

physico_stats <- lapply(feature_cols, function(feat) {
  result <- run_diagnosis_stats(physico_data, feat)
  result$feature <- feat
  result$chain <- gsub("_.*", "", feat)
  result
}) |> dplyr::bind_rows()

physico_stats <- physico_stats |>
  dplyr::group_by(feature) |>
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) |>
  dplyr::ungroup()

# Top differential features boxplots
top_features <- physico_stats |>
  dplyr::filter(!is.na(p.adj)) |>
  dplyr::group_by(feature) |>
  dplyr::summarize(min_p = min(p.adj, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(min_p) |>
  dplyr::slice(1:9) |>
  dplyr::pull(feature)

physico_plots <- lapply(top_features, function(feat) {
  feat_stats <- physico_stats |> dplyr::filter(feature == feat)
  boxplot_with_signif(
    physico_data, y_var = feat, stats_df = feat_stats,
    colors = diagnosis_col, title = gsub("_", " ", feat)
  )
})

physico_combined <- wrap_plots(physico_plots, ncol = 3)
ggsave("results/tcr_comparison/physicochemical/physicochemical_boxplots.pdf",
       physico_combined, width = 14, height = 12)

# PCA by diagnosis (sample-level)
message("  Creating PCA")
sample_physico <- physico_data |>
  dplyr::group_by(sample, patient, tissue, diagnosis, group) |>
  dplyr::summarize(
    dplyr::across(all_of(feature_cols), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::filter(complete.cases(dplyr::across(all_of(feature_cols))))

physico_mat <- sample_physico |> dplyr::select(all_of(feature_cols)) |> as.matrix()
var_cols <- apply(physico_mat, 2, var, na.rm = TRUE) > 0
physico_mat <- physico_mat[, var_cols]

pca_result <- prcomp(physico_mat, scale. = TRUE, center = TRUE)
pca_df <- data.frame(
  PC1 = pca_result$x[,1], PC2 = pca_result$x[,2],
  diagnosis = sample_physico$diagnosis, tissue = sample_physico$tissue
)
var_exp <- summary(pca_result)$importance[2, 1:2] * 100

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diagnosis, shape = tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = diagnosis), level = 0.95, linetype = "dashed") +
  scale_color_manual(values = diagnosis_col) +
  theme_pub() +
  labs(
    title = "Physicochemical properties PCA (sample-level)",
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")
  )

ggsave("results/tcr_comparison/physicochemical/pca_by_diagnosis.pdf", p_pca, width = 9, height = 7)

# Physico heatmap of significant features
sig_physico <- physico_stats |>
  dplyr::filter(p.adj < 0.1) |>
  dplyr::select(feature, diagnosis, estimate) |>
  tidyr::pivot_wider(names_from = diagnosis, values_from = estimate, values_fill = 0) |>
  tibble::column_to_rownames("feature") |>
  as.matrix()

if (nrow(sig_physico) > 2) {
  pdf("results/tcr_comparison/physicochemical/significant_features_heatmap.pdf", width = 10, height = 12)
  pheatmap::pheatmap(
    sig_physico,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "Significant physicochemical features (estimate vs CTRL)",
    cluster_cols = FALSE
  )
  dev.off()
}

writexl::write_xlsx(list(statistics = physico_stats), "results/tcr_comparison/physicochemical/physicochemical_stats.xlsx")

# ============================================================================
# Section 6: Repertoire diversity
# ============================================================================
message("Section 6: Repertoire diversity")

# Compute diversity metrics per sample
message("  Computing diversity metrics")
diversity_data <- tryCatch({
  div_table <- clonalDiversity(sc_tcr, cloneCall = "aa", group.by = "sample", exportTable = TRUE)

  div_table |>
    dplyr::mutate(
      tissue = gsub("_.*", "", sample),
      patient = gsub(".*_", "", sample)
    ) |>
    dplyr::left_join(
      tcr_metadata |> dplyr::distinct(patient, diagnosis, group),
      by = "patient"
    )
}, error = function(e) {
  message("  clonalDiversity failed: ", e$message)
  NULL
})

if (!is.null(diversity_data)) {
  metrics_available <- intersect(c("shannon", "inv.simpson", "chao1"), colnames(diversity_data))

  # Stats by diagnosis (mixed model with tissue)
  message("  Testing diversity by diagnosis")
  diversity_stats_dx <- list()
  for (metric in metrics_available) {
    result <- run_diagnosis_stats(diversity_data, metric, include_tissue = TRUE)
    result$metric <- metric
    result <- result |> dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))
    diversity_stats_dx[[metric]] <- result
  }
  diversity_stats_dx_combined <- dplyr::bind_rows(diversity_stats_dx)

  # Stats by tissue (paired test within patients)
  message("  Testing diversity by tissue")
  diversity_stats_tissue <- lapply(metrics_available, function(metric) {
    paired_data <- diversity_data |>
      dplyr::select(patient, tissue, diagnosis, !!sym(metric)) |>
      tidyr::pivot_wider(names_from = tissue, values_from = !!sym(metric)) |>
      dplyr::filter(!is.na(CSF), !is.na(PBMC))

    if (nrow(paired_data) < 5) {
      return(data.frame(metric = metric, p.value = NA, method = "insufficient_pairs"))
    }

    test <- wilcox.test(paired_data$CSF, paired_data$PBMC, paired = TRUE)
    data.frame(
      metric = metric,
      mean_CSF = mean(paired_data$CSF, na.rm = TRUE),
      mean_PBMC = mean(paired_data$PBMC, na.rm = TRUE),
      diff = mean(paired_data$CSF - paired_data$PBMC, na.rm = TRUE),
      p.value = test$p.value,
      n_pairs = nrow(paired_data),
      method = "paired_wilcoxon"
    )
  }) |> dplyr::bind_rows()

  # Tissue × diagnosis interaction
  message("  Testing tissue × diagnosis interaction")
  diversity_interaction <- lapply(metrics_available, function(metric) {
    tryCatch({
      model <- lmerTest::lmer(
        as.formula(paste0(metric, " ~ diagnosis * tissue + (1|patient)")),
        data = diversity_data
      )
      anova_res <- as.data.frame(anova(model))

      data.frame(
        metric = metric,
        diagnosis_F = anova_res["diagnosis", "F value"],
        diagnosis_p = anova_res["diagnosis", "Pr(>F)"],
        tissue_F = anova_res["tissue", "F value"],
        tissue_p = anova_res["tissue", "Pr(>F)"],
        interaction_F = anova_res["diagnosis:tissue", "F value"],
        interaction_p = anova_res["diagnosis:tissue", "Pr(>F)"]
      )
    }, error = function(e) {
      data.frame(metric = metric, diagnosis_F = NA, diagnosis_p = NA,
                 tissue_F = NA, tissue_p = NA, interaction_F = NA, interaction_p = NA)
    })
  }) |> dplyr::bind_rows()

  # Boxplots by diagnosis with stats
  if (length(metrics_available) > 0) {
    div_dx_plots <- lapply(metrics_available, function(metric) {
      tryCatch({
        boxplot_with_signif(
          diversity_data, y_var = metric,
          stats_df = diversity_stats_dx[[metric]],
          colors = diagnosis_col,
          title = paste(metric, "by diagnosis")
        )
      }, error = function(e) NULL)
    })
    div_dx_plots <- Filter(Negate(is.null), div_dx_plots)

    if (length(div_dx_plots) > 0) {
      div_dx_combined <- wrap_plots(div_dx_plots, ncol = length(div_dx_plots))
      ggsave("results/tcr_comparison/clonality/diversity_by_diagnosis.pdf",
             div_dx_combined, width = 4 * length(div_dx_plots), height = 5)
    }

    # Boxplots by tissue (paired)
    div_tissue_plots <- lapply(metrics_available, function(metric) {
      tryCatch({
        ggplot(diversity_data, aes(x = tissue, y = .data[[metric]], fill = tissue)) +
          geom_boxplot(alpha = 0.7, outlier.shape = NA) +
          geom_line(aes(group = patient), alpha = 0.3, color = "grey50") +
          geom_point(aes(color = diagnosis), size = 2, alpha = 0.7) +
          scale_fill_manual(values = tissue_col) +
          scale_color_manual(values = diagnosis_col) +
          theme_pub() +
          labs(title = paste(metric, "CSF vs PBMC"), y = metric) +
          theme(legend.position = "right")
      }, error = function(e) NULL)
    })
    div_tissue_plots <- Filter(Negate(is.null), div_tissue_plots)

    if (length(div_tissue_plots) > 0) {
      div_tissue_combined <- wrap_plots(div_tissue_plots, ncol = length(div_tissue_plots))
      ggsave("results/tcr_comparison/clonality/diversity_by_tissue_paired.pdf",
             div_tissue_combined, width = 5 * length(div_tissue_plots), height = 5)
    }

    # Faceted plot: diagnosis × tissue
    div_facet_plots <- lapply(metrics_available, function(metric) {
      tryCatch({
        ggplot(diversity_data, aes(x = diagnosis, y = .data[[metric]], fill = diagnosis)) +
          geom_boxplot(alpha = 0.7, outlier.shape = NA) +
          geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
          facet_wrap(~ tissue) +
          scale_fill_manual(values = diagnosis_col) +
          theme_pub() +
          labs(title = paste(metric, "by diagnosis and tissue"), y = metric) +
          theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
      }, error = function(e) NULL)
    })
    div_facet_plots <- Filter(Negate(is.null), div_facet_plots)

    if (length(div_facet_plots) > 0) {
      div_facet_combined <- wrap_plots(div_facet_plots, ncol = 1)
      ggsave("results/tcr_comparison/clonality/diversity_diagnosis_tissue_facet.pdf",
             div_facet_combined, width = 12, height = 4 * length(div_facet_plots))
    }
  }

  # Print summary
  message("\n  Diversity summary:")
  for (metric in metrics_available) {
    int_res <- diversity_interaction |> dplyr::filter(metric == !!metric)
    message(paste("   ", metric, ":"))
    message(paste("      diagnosis p =", round(int_res$diagnosis_p, 4)))
    message(paste("      tissue p =", round(int_res$tissue_p, 4)))
    message(paste("      interaction p =", round(int_res$interaction_p, 4)))
  }

  # Export
  writexl::write_xlsx(
    list(
      by_diagnosis = diversity_stats_dx_combined,
      by_tissue = diversity_stats_tissue,
      interaction = diversity_interaction,
      raw_data = diversity_data
    ),
    "results/tcr_comparison/clonality/diversity_statistics.xlsx"
  )
}

# ============================================================================
# Section 7: Sequence similarity networks
# ============================================================================
message("Section 7: Sequence similarity")

similarity_results <- list()

for (chain in chains) {
  message(paste("  Building", chain, "similarity network"))

  seq_col <- switch(chain, TRA = "TRA_CDR3", TRB = "TRB_CDR3", TRAB = "CTaa")

  chain_data <- tcr_metadata |> dplyr::filter(!is.na(.data[[seq_col]]))

  if (nrow(chain_data) < 100) next

  edge_list <- tryCatch({
    immApex::buildNetwork(chain_data, seq_col = seq_col, threshold = 3)
  }, error = function(e) NULL)

  if (is.null(edge_list) || nrow(edge_list) == 0) next

  edge_patients <- data.frame(
    from = chain_data[[seq_col]][as.numeric(edge_list$from)],
    to = chain_data[[seq_col]][as.numeric(edge_list$to)],
    dist = edge_list$dist
  ) |>
    dplyr::left_join(
      chain_data |> dplyr::distinct(.data[[seq_col]], patient) |> dplyr::rename(patient_from = patient),
      by = setNames(seq_col, "from")
    ) |>
    dplyr::left_join(
      chain_data |> dplyr::distinct(.data[[seq_col]], patient) |> dplyr::rename(patient_to = patient),
      by = setNames(seq_col, "to")
    ) |>
    dplyr::filter(patient_from != patient_to)

  patient_sim <- edge_patients |>
    dplyr::count(patient_from, patient_to, name = "n_similar") |>
    tidyr::pivot_wider(names_from = patient_to, values_from = n_similar, values_fill = 0) |>
    tibble::column_to_rownames("patient_from") |>
    as.matrix()

  all_patients <- sort(unique(c(rownames(patient_sim), colnames(patient_sim))))
  full_matrix <- matrix(0, length(all_patients), length(all_patients),
                        dimnames = list(all_patients, all_patients))
  full_matrix[rownames(patient_sim), colnames(patient_sim)] <- patient_sim
  full_matrix <- full_matrix + t(full_matrix)
  diag(full_matrix) <- 0

  similarity_results[[chain]] <- as.data.frame(full_matrix) |>
    tibble::rownames_to_column("patient")

  patient_anno <- tcr_metadata |>
    dplyr::distinct(patient, diagnosis, group) |>
    dplyr::filter(patient %in% all_patients) |>
    tibble::column_to_rownames("patient")

  pdf(paste0("results/tcr_comparison/similarity/patient_similarity_", chain, ".pdf"),
      width = 10, height = 8)
  pheatmap::pheatmap(
    full_matrix,
    color = colorRampPalette(c("white", "red3"))(50),
    main = paste(chain, "similar sequences shared between patients"),
    display_numbers = FALSE,
    annotation_row = patient_anno,
    annotation_col = patient_anno,
    annotation_colors = list(diagnosis = diagnosis_col, group = group_col)
  )
  dev.off()
}

if (length(similarity_results) > 0) {
  writexl::write_xlsx(similarity_results, "results/tcr_comparison/similarity/patient_similarity.xlsx")
}

# ============================================================================
# Section 8: Protein embeddings
# ============================================================================
tryCatch({
  # --- Step 1: Generate Single Chain Embeddings ---
  message("  Generating embeddings with immLynx")
  
  # Ensure we have single chain embeddings first
  for (chain_type in c("TRB", "TRA")) {
    reduction_name <- paste0("esm_", tolower(chain_type))
    
    # Check if already computed to save time
    # Use function without namespace prefix to ensure devtools::load_all version is used
    if (!reduction_name %in% names(sc_tcr@reductions)) {
      sc_tcr <- runEmbeddings(
        sc_tcr,
        chains = chain_type,
        model_name = "facebook/esm2_t12_35M_UR50D",
        pool = "mean",
        chunk_size = 32,
        reduction_name = reduction_name
      )
      message(paste("    ", chain_type, "embeddings generated"))
    }
  }
  
  # --- Step 2: Create Paired (TRA+TRB) Embedding ---
  message("  Creating Joint TRA+TRB Embeddings")
  if ("esm_tra" %in% names(sc_tcr@reductions) & "esm_trb" %in% names(sc_tcr@reductions)) {
    
    # Extract matrices
    emb_tra <- Embeddings(sc_tcr, "esm_tra")
    emb_trb <- Embeddings(sc_tcr, "esm_trb")
    
    # Ensure cell matching (intersection)
    common_cells <- intersect(rownames(emb_tra), rownames(emb_trb))
    
    # Concatenate columns (Joint Embedding)
    joint_emb <- cbind(emb_tra[common_cells, ], emb_trb[common_cells, ])
    
    # Add to Seurat object
    sc_tcr[["esm_joint"]] <- joint_emb
    message("    Joint embeddings created")
  }
  
  # --- Step 3: Loop Through All Modalities (TRA, TRB, JOINT) ---
  analysis_targets <- c("esm_tra", "esm_trb", "esm_joint")
  esm_permanova_results <- list()

  for (red_name in analysis_targets) {
    if (!red_name %in% names(sc_tcr@reductions)) next
    
    message(paste("  Analyzing:", red_name))
    
    # Extract data
    embeddings <- na.omit(Embeddings(sc_tcr, reduction = red_name))
    meta_sub <- sc_tcr[[]][rownames(embeddings), ]
    
    # A. QUANTITATIVE: Distance to CTRL Centroid
    # We calculate the mean vector (centroid) for CTRL and each diagnosis
    # Then measure Euclidean distance in the high-dim ESM space
    df_emb <- as.data.frame(embeddings)
    df_emb$diagnosis <- meta_sub$diagnosis
    
    centroids <- df_emb %>%
      group_by(diagnosis) %>%
      summarise(across(everything(), mean))
    
    # Isolate CTRL centroid
    if ("CTRL" %in% centroids$diagnosis) {
      ctrl_vec <- as.numeric(centroids[centroids$diagnosis == "CTRL", -1])
      
      # Calculate distance of other diagnoses to CTRL
      dist_res <- centroids %>%
        filter(diagnosis != "CTRL") %>%
        rowwise() %>%
        mutate(
          dist_to_ctrl = sqrt(sum((c_across(-diagnosis) - ctrl_vec)^2)),
          modality = red_name
        ) %>%
        select(diagnosis, dist_to_ctrl, modality)
      
      print(dist_res) # Output distances to console
      write.csv(dist_res, paste0("results/tcr_comparison/embeddings/dist_to_ctrl_", red_name, ".csv"))
    }
    
    # B. STATISTICAL: PERMANOVA (Multivariate ANOVA) - sample-level aggregation
    set.seed(42)

    # Aggregate embeddings to sample level (mean embedding per sample)
    emb_cols <- colnames(embeddings)
    emb_df <- as.data.frame(embeddings)
    emb_df$sample <- meta_sub$sample
    emb_df$patient <- meta_sub$patient
    emb_df$tissue <- meta_sub$tissue
    emb_df$diagnosis <- meta_sub$diagnosis

    sample_emb <- emb_df |>
      dplyr::group_by(sample, patient, tissue, diagnosis) |>
      dplyr::summarize(
        dplyr::across(all_of(emb_cols), ~ mean(.x, na.rm = TRUE)),
        n_cells = n(),
        .groups = "drop"
      ) |>
      dplyr::filter(complete.cases(dplyr::across(all_of(emb_cols))))

    emb_matrix <- sample_emb |> dplyr::select(all_of(emb_cols)) |> as.matrix()

    message(paste("    Sample-level aggregation:", nrow(emb_matrix), "samples,", ncol(emb_matrix), "dimensions"))

    # Run PERMANOVA with patient stratification
    permanova_result <- run_permanova(emb_matrix, sample_emb, "~ diagnosis * tissue", "patient", n_perm = 999)

    if (permanova_result$success) {
      permanova <- permanova_result$result
      diag_r2 <- permanova["diagnosis", "R2"]
      diag_p <- permanova["diagnosis", "Pr(>F)"]
      tissue_r2 <- permanova["tissue", "R2"]
      tissue_p <- permanova["tissue", "Pr(>F)"]
      message(paste("    diagnosis: R2 =", round(diag_r2, 4), ", p =", round(diag_p, 4)))
      message(paste("    tissue: R2 =", round(tissue_r2, 4), ", p =", round(tissue_p, 4)))
      esm_permanova_results[[red_name]] <- permanova
    } else {
      message(paste("    PERMANOVA failed:", permanova_result$message))
      permanova <- NULL
      diag_r2 <- NA
    }

    capture.output(permanova, file = paste0("results/tcr_comparison/embeddings/permanova_", red_name, ".txt"))
    
    # C. VISUALIZATION: UMAP with Density (cell-level for visualization)
    set.seed(42)

    # Subsample for UMAP if too many cells (visualization only)
    max_cells_umap <- 5000
    if (nrow(embeddings) > max_cells_umap) {
      idx_vis <- sample(nrow(embeddings), max_cells_umap)
      emb_vis <- embeddings[idx_vis, ]
      meta_vis <- meta_sub[idx_vis, ]
    } else {
      emb_vis <- embeddings
      meta_vis <- meta_sub
    }

    umap_res <- uwot::umap(emb_vis, n_neighbors = 15, min_dist = 0.1)

    plot_df <- data.frame(
      UMAP1 = umap_res[,1],
      UMAP2 = umap_res[,2],
      diagnosis = meta_vis$diagnosis,
      tissue = meta_vis$tissue
    )

    subtitle_text <- if (!is.na(diag_r2)) paste("PERMANOVA R2 =", round(diag_r2, 3), "(sample-level)") else "PERMANOVA failed"

    p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = diagnosis)) +
      geom_point(alpha = 0.4, size = 0.5) +
      stat_density_2d(aes(alpha = after_stat(level)), bins = 5, color = "black", linewidth = 0.2) +
      scale_color_manual(values = diagnosis_col) +
      scale_alpha_continuous(range = c(0.1, 0.3), guide = "none") +
      theme_pub() +
      labs(
        title = paste("ESM2 UMAP:", toupper(gsub("esm_", "", red_name))),
        subtitle = subtitle_text
      ) +
      facet_wrap(~diagnosis)

    ggsave(paste0("results/tcr_comparison/embeddings/umap_density_", red_name, ".pdf"), p, width = 12, height = 8)

    # Also create sample-level PCA for direct comparison with PERMANOVA
    pca_emb <- prcomp(emb_matrix, scale. = TRUE, center = TRUE)
    var_exp <- summary(pca_emb)$importance[2, 1:2] * 100

    pca_df <- data.frame(
      PC1 = pca_emb$x[,1],
      PC2 = pca_emb$x[,2],
      diagnosis = sample_emb$diagnosis,
      tissue = sample_emb$tissue,
      patient = sample_emb$patient
    )

    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diagnosis, shape = tissue)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(aes(group = diagnosis), level = 0.95, linetype = "dashed") +
      scale_color_manual(values = diagnosis_col) +
      theme_pub() +
      labs(
        title = paste("ESM2 PCA (sample-level):", toupper(gsub("esm_", "", red_name))),
        subtitle = subtitle_text,
        x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PC2 (", round(var_exp[2], 1), "%)")
      )

    ggsave(paste0("results/tcr_comparison/embeddings/pca_sample_", red_name, ".pdf"), p_pca, width = 9, height = 7)
  }
  
  message("  Embeddings analysis complete")
  
}, error = function(e) {
  message("  Embeddings failed: ", e$message)
})

# ============================================================================
# Section 9: TCRdist analysis
# ============================================================================
message("Section 9: TCRdist analysis")

dir.create("results/tcr_comparison/tcrdist", showWarnings = FALSE, recursive = TRUE)

tcrdist_results <- tryCatch({
  message("  Computing TCRdist distances")

  # Run TCRdist for beta chain (most informative for antigen specificity)
  # Subsample to max 5000 sequences to avoid memory issues with large distance matrices
  # Use runTCRdist without namespace prefix to ensure devtools::load_all version is used
  dist_trb <- runTCRdist(
    sc_tcr,
    chains = "beta",
    organism = "human",
    compute_distances = TRUE,
    max_sequences = 5000,
    add_to_object = FALSE
  )
  message("    TRB distances computed: ", nrow(dist_trb$tcr_data), " sequences")

  # Run TCRdist for alpha chain
  dist_tra <- runTCRdist(
    sc_tcr,
    chains = "alpha",
    organism = "human",
    compute_distances = TRUE,
    max_sequences = 5000,
    add_to_object = FALSE
  )
  message("    TRA distances computed: ", nrow(dist_tra$tcr_data), " sequences")

  # Get metadata for the cells in the distance matrix
  trb_cells <- dist_trb$barcodes
  tra_cells <- dist_tra$barcodes

  # Sample-level aggregation of TCRdist: mean distance to other samples
  message("  Aggregating TCRdist to sample level")

  aggregate_tcrdist <- function(dist_mat, barcodes, metadata) {
    # Match barcodes to metadata - try multiple strategies
    # Strategy 1: Direct match to rownames
    matched_idx <- match(barcodes, rownames(metadata))

    # Strategy 2: If no matches, barcodes might be in a cell_id column
    if (all(is.na(matched_idx)) && "cell_id" %in% colnames(metadata)) {
      matched_idx <- match(barcodes, metadata$cell_id)
    }

    # Strategy 3: Try matching partial barcodes (strip sample prefix)
    if (all(is.na(matched_idx))) {
      # Extract core barcode (last part after underscore or dash)
      core_barcodes <- sub(".*[_-]", "", barcodes)
      core_rownames <- sub(".*[_-]", "", rownames(metadata))
      matched_idx <- match(core_barcodes, core_rownames)
    }

    # Strategy 4: Try matching barcodes to rownames allowing partial matches
    if (all(is.na(matched_idx))) {
      # Try substring matching
      matched_idx <- sapply(barcodes, function(bc) {
        idx <- grep(bc, rownames(metadata), fixed = TRUE)
        if (length(idx) == 1) idx else NA
      })
    }

    # Report match rate
    match_rate <- sum(!is.na(matched_idx)) / length(barcodes)
    message(paste("    Barcode match rate:", round(match_rate * 100, 1), "%"))

    if (match_rate < 0.1) {
      warning("Very low barcode match rate. Check barcode formats.")
      # Debug: show examples
      message("    Example barcodes from TCRdist: ", paste(head(barcodes, 3), collapse = ", "))
      message("    Example rownames from metadata: ", paste(head(rownames(metadata), 3), collapse = ", "))
    }

    # Filter to matched cells only
    valid_idx <- which(!is.na(matched_idx))
    if (length(valid_idx) < 10) {
      stop("Too few barcodes matched to metadata (n=", length(valid_idx), ")")
    }

    # Subset distance matrix to valid cells
    dist_mat <- dist_mat[valid_idx, valid_idx]
    matched_idx <- matched_idx[valid_idx]

    # Get sample assignment for each barcode (now indexed 1:length(valid_idx))
    cell_samples <- metadata$sample[matched_idx]

    # Remove NA samples
    valid_samples <- !is.na(cell_samples)
    if (sum(valid_samples) < 10) {
      stop("Too few cells with valid sample assignment")
    }
    dist_mat <- dist_mat[valid_samples, valid_samples]
    cell_samples <- cell_samples[valid_samples]

    # Calculate mean within-sample and between-sample distances
    samples <- unique(cell_samples)
    n_samples <- length(samples)
    message(paste("    Aggregating", sum(valid_samples), "cells into", n_samples, "samples"))

    sample_dist <- matrix(0, n_samples, n_samples, dimnames = list(samples, samples))

    for (i in seq_along(samples)) {
      for (j in seq_along(samples)) {
        # Now idx_i and idx_j are relative to the current dist_mat
        idx_i <- which(cell_samples == samples[i])
        idx_j <- which(cell_samples == samples[j])
        if (length(idx_i) > 0 && length(idx_j) > 0) {
          sample_dist[i, j] <- mean(dist_mat[idx_i, idx_j], na.rm = TRUE)
        }
      }
    }
    sample_dist
  }

  # Aggregate TRB distances
  trb_sample_dist <- NULL
  if (!is.null(dist_trb$distances$pw_beta)) {
    message("  Aggregating TRB distances...")
    dist_mat_trb <- as.matrix(dist_trb$distances$pw_beta)
    # Sync barcodes with matrix dims (tcrdist3 may filter some sequences)
    n_dist <- nrow(dist_mat_trb)
    barcodes_trb <- if (length(dist_trb$barcodes) > n_dist) {
      dist_trb$barcodes[1:n_dist]
    } else {
      dist_trb$barcodes
    }
    message("    Matrix: ", n_dist, " x ", ncol(dist_mat_trb), ", Barcodes: ", length(barcodes_trb))
    trb_sample_dist <- aggregate_tcrdist(dist_mat_trb, barcodes_trb, sc_tcr@meta.data)
  } else {
    message("  Warning: No TRB distance matrix available")
  }

  # Aggregate TRA distances
  tra_sample_dist <- NULL
  if (!is.null(dist_tra$distances$pw_alpha)) {
    message("  Aggregating TRA distances...")
    dist_mat_tra <- as.matrix(dist_tra$distances$pw_alpha)
    # Sync barcodes with matrix dims
    n_dist <- nrow(dist_mat_tra)
    barcodes_tra <- if (length(dist_tra$barcodes) > n_dist) {
      dist_tra$barcodes[1:n_dist]
    } else {
      dist_tra$barcodes
    }
    message("    Matrix: ", n_dist, " x ", ncol(dist_mat_tra), ", Barcodes: ", length(barcodes_tra))
    tra_sample_dist <- aggregate_tcrdist(dist_mat_tra, barcodes_tra, sc_tcr@meta.data
    )
  } else {
    message("  Warning: No TRA distance matrix available")
  }

  # PERMANOVA on sample-level TCRdist
  message("  Running PERMANOVA on TCRdist")

  tcrdist_permanova <- list()

  for (chain_name in c("TRB", "TRA")) {
    dist_mat <- if (chain_name == "TRB") trb_sample_dist else tra_sample_dist

    # Skip if distance matrix not available
    if (is.null(dist_mat)) {
      message(paste("    Skipping", chain_name, "- no distance matrix"))
      next
    }

    # Get metadata for samples
    sample_meta <- sample_summary |>
      dplyr::filter(sample %in% rownames(dist_mat)) |>
      dplyr::arrange(match(sample, rownames(dist_mat)))

    if (nrow(sample_meta) < 3) {
      message(paste("    Skipping", chain_name, "- too few samples (n=", nrow(sample_meta), ")"))
      next
    }

    # Ensure order matches
    dist_mat <- dist_mat[sample_meta$sample, sample_meta$sample]

    # Convert to dist object
    dist_obj <- as.dist(dist_mat)

    perm_result <- tryCatch({
      vegan::adonis2(
        dist_obj ~ diagnosis * tissue,
        data = sample_meta,
        permutations = 999,
        strata = sample_meta$patient
      )
    }, error = function(e) {
      message(paste("    PERMANOVA failed for", chain_name, ":", e$message))
      NULL
    })

    if (!is.null(perm_result)) {
      tcrdist_permanova[[chain_name]] <- perm_result
      message(paste("   ", chain_name, "- diagnosis R2:",
                    round(perm_result["diagnosis", "R2"], 4),
                    "p:", round(perm_result["diagnosis", "Pr(>F)"], 4)))
    }
  }

  # Visualize sample-level TCRdist as heatmap
  message("  Creating TCRdist heatmaps")

  for (chain_name in c("TRB", "TRA")) {
    dist_mat <- if (chain_name == "TRB") trb_sample_dist else tra_sample_dist

    # Skip if distance matrix not available
    if (is.null(dist_mat)) next

    # Get sample annotations
    sample_anno <- sample_summary |>
      dplyr::filter(sample %in% rownames(dist_mat)) |>
      dplyr::select(sample, diagnosis, tissue) |>
      tibble::column_to_rownames("sample")

    if (nrow(sample_anno) < 2) next

    # Reorder to match
    dist_mat <- dist_mat[rownames(sample_anno), rownames(sample_anno)]

    pdf(paste0("results/tcr_comparison/tcrdist/sample_distance_", chain_name, ".pdf"),
        width = 10, height = 8)
    pheatmap::pheatmap(
      dist_mat,
      color = colorRampPalette(c("darkblue", "white", "darkred"))(100),
      main = paste(chain_name, "TCRdist (sample-level mean)"),
      annotation_row = sample_anno,
      annotation_col = sample_anno,
      annotation_colors = list(diagnosis = diagnosis_col, tissue = tissue_col),
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize = 8
    )
    dev.off()
  }

  # PCoA visualization
  message("  Creating PCoA plots")

  pcoa_results <- list()
  for (chain_name in c("TRB", "TRA")) {
    dist_mat <- if (chain_name == "TRB") trb_sample_dist else tra_sample_dist

    # Skip if distance matrix not available
    if (is.null(dist_mat)) next

    # Get sample metadata
    sample_meta <- sample_summary |>
      dplyr::filter(sample %in% rownames(dist_mat))

    if (nrow(sample_meta) < 3) next

    dist_mat <- dist_mat[sample_meta$sample, sample_meta$sample]

    # PCoA
    pcoa_result <- tryCatch({
      cmdscale(as.dist(dist_mat), k = min(5, nrow(dist_mat) - 1), eig = TRUE)
    }, error = function(e) NULL)

    if (is.null(pcoa_result)) next

    pos_eig <- pcoa_result$eig[pcoa_result$eig > 0]
    var_exp <- pcoa_result$eig[1:min(5, length(pos_eig))] / sum(pos_eig) * 100

    pcoa_df <- data.frame(
      sample = sample_meta$sample,
      PCo1 = pcoa_result$points[,1],
      PCo2 = pcoa_result$points[,2],
      diagnosis = sample_meta$diagnosis,
      tissue = sample_meta$tissue,
      patient = sample_meta$patient
    )
    pcoa_results[[chain_name]] <- pcoa_df

    # Basic PCoA plot by diagnosis and tissue
    p <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = diagnosis, shape = tissue)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(aes(group = diagnosis), level = 0.95, linetype = "dashed") +
      scale_color_manual(values = diagnosis_col) +
      theme_pub() +
      labs(
        title = paste(chain_name, "TCRdist PCoA (sample-level)"),
        x = paste0("PCo1 (", round(var_exp[1], 1), "%)"),
        y = paste0("PCo2 (", round(var_exp[2], 1), "%)")
      )

    ggsave(paste0("results/tcr_comparison/tcrdist/pcoa_", chain_name, ".pdf"), p, width = 9, height = 7)
  }

  # =====================================================
  # DIAGNOSIS-FOCUSED VISUALIZATIONS
  # =====================================================
  message("  Creating diagnosis-focused visualizations")

  # 1. Distance to controls analysis
  # For each sample, calculate mean TCRdist to control samples
  message("    Calculating distance to controls")

  dist_to_ctrl_results <- list()
  for (chain_name in c("TRB", "TRA")) {
    dist_mat <- if (chain_name == "TRB") trb_sample_dist else tra_sample_dist
    if (is.null(dist_mat)) next

    sample_meta <- sample_summary |>
      dplyr::filter(sample %in% rownames(dist_mat))

    ctrl_samples <- sample_meta$sample[sample_meta$diagnosis == "CTRL"]
    disease_samples <- sample_meta$sample[sample_meta$diagnosis != "CTRL"]

    if (length(ctrl_samples) < 2 || length(disease_samples) < 2) next

    # Calculate mean distance to controls for each sample
    dist_to_ctrl <- sapply(rownames(dist_mat), function(s) {
      if (s %in% ctrl_samples) {
        # For controls, calculate distance to other controls
        other_ctrl <- setdiff(ctrl_samples, s)
        mean(dist_mat[s, other_ctrl], na.rm = TRUE)
      } else {
        # For disease, calculate distance to all controls
        mean(dist_mat[s, ctrl_samples], na.rm = TRUE)
      }
    })

    dist_ctrl_df <- data.frame(
      sample = names(dist_to_ctrl),
      dist_to_ctrl = dist_to_ctrl
    ) |>
      dplyr::left_join(sample_meta, by = "sample") |>
      dplyr::mutate(
        is_ctrl = diagnosis == "CTRL",
        diagnosis_group = ifelse(diagnosis == "CTRL", "CTRL", "Disease")
      )

    dist_to_ctrl_results[[chain_name]] <- dist_ctrl_df

    # Statistical test: Wilcoxon test for each diagnosis vs CTRL
    ctrl_dists <- dist_ctrl_df$dist_to_ctrl[dist_ctrl_df$diagnosis == "CTRL"]
    stat_results <- lapply(setdiff(unique(dist_ctrl_df$diagnosis), "CTRL"), function(dx) {
      dx_dists <- dist_ctrl_df$dist_to_ctrl[dist_ctrl_df$diagnosis == dx]
      if (length(dx_dists) < 2) return(NULL)
      test <- wilcox.test(dx_dists, ctrl_dists, alternative = "greater")
      data.frame(
        chain = chain_name,
        diagnosis = dx,
        n_samples = length(dx_dists),
        median_dist = median(dx_dists),
        ctrl_median = median(ctrl_dists),
        fold_change = median(dx_dists) / median(ctrl_dists),
        p.value = test$p.value
      )
    }) |> dplyr::bind_rows()

    if (nrow(stat_results) > 0) {
      stat_results$p.adj <- p.adjust(stat_results$p.value, method = "BH")
      dist_to_ctrl_results[[paste0(chain_name, "_stats")]] <- stat_results

      # Report significant findings
      sig <- stat_results |> dplyr::filter(p.adj < 0.1)
      if (nrow(sig) > 0) {
        message("      Significantly divergent from CTRL (FDR < 0.1):")
        for (i in seq_len(nrow(sig))) {
          message(sprintf("        %s: FC=%.2f, p.adj=%.4f",
                          sig$diagnosis[i], sig$fold_change[i], sig$p.adj[i]))
        }
      }
    }

    # Use consistent diagnosis ordering (CTRL first)
    dx_levels <- c("CTRL", intersect(neuropathy_dx, unique(dist_ctrl_df$diagnosis)))
    dist_ctrl_df$diagnosis <- factor(dist_ctrl_df$diagnosis, levels = dx_levels)

    # Plot: Distance to controls by diagnosis with significance annotations
    p_dist <- ggplot(dist_ctrl_df, aes(x = diagnosis, y = dist_to_ctrl, fill = diagnosis)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(aes(shape = tissue), width = 0.2, size = 2, alpha = 0.8) +
      scale_fill_manual(values = diagnosis_col) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Add significance brackets from Wilcoxon tests
    if (nrow(stat_results) > 0) {
      sig_stats <- stat_results |> dplyr::filter(p.adj < 0.05)
      if (nrow(sig_stats) > 0) {
        comparisons <- lapply(sig_stats$diagnosis, function(dx) c("CTRL", dx))
        annotations <- sapply(sig_stats$p.adj, format_pval)
        p_dist <- p_dist + ggsignif::geom_signif(
          comparisons = comparisons, annotations = annotations,
          step_increase = 0.08, textsize = 3, vjust = 0.5
        )
      }
    }

    p_dist <- p_dist +
      labs(
        title = paste(chain_name, "- TCR distance to healthy controls"),
        x = "Diagnosis",
        y = "Mean TCRdist to CTRL samples"
      )
    ggsave(paste0("results/tcr_comparison/tcrdist/dist_to_ctrl_", chain_name, ".pdf"),
           p_dist, width = 10, height = 7)

    # Plot: Faceted by tissue with per-tissue statistics
    tissue_stats <- lapply(c("CSF", "PBMC"), function(tis) {
      tis_df <- dist_ctrl_df |> dplyr::filter(tissue == tis)
      ctrl_d <- tis_df$dist_to_ctrl[tis_df$diagnosis == "CTRL"]
      if (length(ctrl_d) < 2) return(NULL)
      lapply(setdiff(levels(tis_df$diagnosis), "CTRL"), function(dx) {
        dx_d <- tis_df$dist_to_ctrl[tis_df$diagnosis == dx]
        if (length(dx_d) < 2) return(NULL)
        test <- tryCatch(wilcox.test(dx_d, ctrl_d, alternative = "greater"), error = function(e) list(p.value = NA))
        data.frame(tissue = tis, diagnosis = dx, p.value = test$p.value)
      }) |> dplyr::bind_rows()
    }) |> dplyr::bind_rows()

    if (nrow(tissue_stats) > 0) {
      tissue_stats$p.adj <- p.adjust(tissue_stats$p.value, method = "BH")
    }

    p_dist_tissue <- ggplot(dist_ctrl_df, aes(x = diagnosis, y = dist_to_ctrl, fill = diagnosis)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      facet_wrap(~tissue, scales = "free_y") +
      scale_fill_manual(values = diagnosis_col) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(
        title = paste(chain_name, "- TCR distance to controls by tissue"),
        x = "Diagnosis",
        y = "Mean TCRdist to CTRL"
      )

    # Add per-tissue significance brackets
    if (nrow(tissue_stats) > 0) {
      for (tis in unique(tissue_stats$tissue)) {
        tis_sig <- tissue_stats |> dplyr::filter(tissue == tis, p.adj < 0.05)
        if (nrow(tis_sig) > 0) {
          comparisons <- lapply(tis_sig$diagnosis, function(dx) c("CTRL", dx))
          annotations <- sapply(tis_sig$p.adj, format_pval)
          p_dist_tissue <- p_dist_tissue + ggsignif::geom_signif(
            comparisons = comparisons, annotations = annotations,
            step_increase = 0.08, textsize = 3, vjust = 0.5,
            data = dist_ctrl_df |> dplyr::filter(tissue == tis)
          )
        }
      }
    }

    ggsave(paste0("results/tcr_comparison/tcrdist/dist_to_ctrl_by_tissue_", chain_name, ".pdf"),
           p_dist_tissue, width = 12, height = 6)
  }

  # 2. Within-diagnosis vs between-diagnosis distances
  message("    Calculating within vs between diagnosis distances")

  within_between_results <- list()
  for (chain_name in c("TRB", "TRA")) {
    dist_mat <- if (chain_name == "TRB") trb_sample_dist else tra_sample_dist
    if (is.null(dist_mat)) next

    sample_meta <- sample_summary |>
      dplyr::filter(sample %in% rownames(dist_mat))

    # Calculate pairwise comparisons
    comparisons <- expand.grid(s1 = rownames(dist_mat), s2 = rownames(dist_mat),
                               stringsAsFactors = FALSE) |>
      dplyr::filter(s1 < s2) |>
      dplyr::mutate(
        distance = mapply(function(a, b) dist_mat[a, b], s1, s2),
        dx1 = sample_meta$diagnosis[match(s1, sample_meta$sample)],
        dx2 = sample_meta$diagnosis[match(s2, sample_meta$sample)],
        tissue1 = sample_meta$tissue[match(s1, sample_meta$sample)],
        tissue2 = sample_meta$tissue[match(s2, sample_meta$sample)],
        patient1 = sample_meta$patient[match(s1, sample_meta$sample)],
        patient2 = sample_meta$patient[match(s2, sample_meta$sample)],
        same_diagnosis = dx1 == dx2,
        same_tissue = tissue1 == tissue2,
        same_patient = patient1 == patient2,
        comparison_type = case_when(
          same_diagnosis & same_patient ~ "Within patient",
          same_diagnosis ~ "Within diagnosis",
          TRUE ~ "Between diagnoses"
        )
      )

    within_between_results[[chain_name]] <- comparisons

    # Wilcoxon test: within vs between diagnosis (excluding same-patient)
    comp_no_self <- comparisons |> dplyr::filter(!same_patient)

    wb_stats <- lapply(c(TRUE, FALSE), function(st) {
      sub <- comp_no_self |> dplyr::filter(same_tissue == st)
      within_d <- sub$distance[sub$comparison_type == "Within diagnosis"]
      between_d <- sub$distance[sub$comparison_type == "Between diagnoses"]
      if (length(within_d) < 3 || length(between_d) < 3) return(NULL)
      test <- wilcox.test(within_d, between_d)
      data.frame(same_tissue = st, p.value = test$p.value,
                 median_within = median(within_d), median_between = median(between_d))
    }) |> dplyr::bind_rows()

    within_between_results[[paste0(chain_name, "_stats")]] <- wb_stats

    # Use consistent colors matching script theme
    wb_fill <- c("Within diagnosis" = "#66C2A5", "Between diagnoses" = "#FC8D62")

    # Plot: Within vs between diagnosis distances with stats
    p_wb <- ggplot(comp_no_self,
                   aes(x = comparison_type, y = distance, fill = comparison_type)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.2, alpha = 0.9, outlier.shape = NA) +
      facet_wrap(~same_tissue, labeller = labeller(same_tissue = c("FALSE" = "Different tissues",
                                                                    "TRUE" = "Same tissue"))) +
      scale_fill_manual(values = wb_fill) +
      theme_pub() +
      theme(legend.position = "none")

    # Add Wilcoxon p-values per facet
    if (nrow(wb_stats) > 0) {
      for (i in seq_len(nrow(wb_stats))) {
        row <- wb_stats[i, ]
        p_wb <- p_wb + ggsignif::geom_signif(
          comparisons = list(c("Within diagnosis", "Between diagnoses")),
          annotations = format_pval(row$p.value),
          textsize = 3.5, vjust = 0.5,
          data = comp_no_self |> dplyr::filter(same_tissue == row$same_tissue)
        )
      }
    }

    p_wb <- p_wb +
      labs(
        title = paste(chain_name, "- Within vs between diagnosis TCR distances"),
        x = "",
        y = "TCRdist"
      )
    ggsave(paste0("results/tcr_comparison/tcrdist/within_between_dx_", chain_name, ".pdf"),
           p_wb, width = 10, height = 6)
  }

  # 3. Diagnosis-specific clustering in PCoA space
  message("    Creating diagnosis-specific PCoA panels")

  for (chain_name in c("TRB", "TRA")) {
    if (!chain_name %in% names(pcoa_results)) next
    pcoa_df <- pcoa_results[[chain_name]]

    # Highlight each diagnosis vs controls
    diagnoses <- setdiff(unique(pcoa_df$diagnosis), "CTRL")

    plots <- lapply(diagnoses, function(dx) {
      df <- pcoa_df |>
        dplyr::mutate(
          highlight = case_when(
            diagnosis == dx ~ dx,
            diagnosis == "CTRL" ~ "CTRL",
            TRUE ~ "Other"
          ),
          highlight = factor(highlight, levels = c(dx, "CTRL", "Other"))
        )

      ggplot(df, aes(x = PCo1, y = PCo2)) +
        geom_point(data = df |> dplyr::filter(highlight == "Other"),
                   color = "gray80", size = 2, alpha = 0.5) +
        geom_point(data = df |> dplyr::filter(highlight == "CTRL"),
                   aes(shape = tissue), color = "gray40", size = 3, alpha = 0.8) +
        geom_point(data = df |> dplyr::filter(highlight == dx),
                   aes(shape = tissue), color = diagnosis_col[dx], size = 3, alpha = 0.9) +
        stat_ellipse(data = df |> dplyr::filter(diagnosis == "CTRL"),
                     color = "gray40", linetype = "dashed") +
        stat_ellipse(data = df |> dplyr::filter(diagnosis == dx),
                     color = diagnosis_col[dx], linetype = "solid") +
        theme_pub() +
        labs(title = paste(dx, "vs CTRL"), x = "PCo1", y = "PCo2") +
        theme(legend.position = "bottom")
    })

    # Combine into multi-panel figure
    p_combined <- patchwork::wrap_plots(plots, ncol = 3) +
      patchwork::plot_annotation(
        title = paste(chain_name, "- Diagnosis-specific TCR profiles vs controls"),
      )
    ggsave(paste0("results/tcr_comparison/tcrdist/pcoa_dx_vs_ctrl_", chain_name, ".pdf"),
           p_combined, width = 14, height = ceiling(length(diagnoses)/3) * 5)
  }

  # 4. Centroid distances between diagnoses
  message("    Calculating diagnosis centroid distances")

  centroid_results <- list()
  for (chain_name in c("TRB", "TRA")) {
    dist_mat <- if (chain_name == "TRB") trb_sample_dist else tra_sample_dist
    if (is.null(dist_mat)) next

    sample_meta <- sample_summary |>
      dplyr::filter(sample %in% rownames(dist_mat))

    diagnoses <- unique(sample_meta$diagnosis)
    n_dx <- length(diagnoses)

    # Calculate mean distance between diagnosis groups
    dx_dist_mat <- matrix(0, n_dx, n_dx, dimnames = list(diagnoses, diagnoses))
    for (i in seq_along(diagnoses)) {
      for (j in seq_along(diagnoses)) {
        samples_i <- sample_meta$sample[sample_meta$diagnosis == diagnoses[i]]
        samples_j <- sample_meta$sample[sample_meta$diagnosis == diagnoses[j]]
        dx_dist_mat[i, j] <- mean(dist_mat[samples_i, samples_j], na.rm = TRUE)
      }
    }

    centroid_results[[chain_name]] <- dx_dist_mat

    # Heatmap of diagnosis-level distances
    pdf(paste0("results/tcr_comparison/tcrdist/diagnosis_centroid_dist_", chain_name, ".pdf"),
        width = 8, height = 7)
    pheatmap::pheatmap(
      dx_dist_mat,
      color = colorRampPalette(c("navy", "white", "firebrick"))(100),
      main = paste(chain_name, "- Mean TCRdist between diagnosis groups"),
      display_numbers = TRUE,
      number_format = "%.1f",
      fontsize_number = 8,
      cluster_rows = TRUE,
      cluster_cols = TRUE
    )
    dev.off()
  }

  # Export results
  message("  Exporting TCRdist results")

  export_list <- list()
  if (!is.null(trb_sample_dist)) {
    export_list$TRB_sample_dist <- as.data.frame(trb_sample_dist) |> tibble::rownames_to_column("sample")
  }
  if (!is.null(tra_sample_dist)) {
    export_list$TRA_sample_dist <- as.data.frame(tra_sample_dist) |> tibble::rownames_to_column("sample")
  }
  if ("TRB" %in% names(dist_to_ctrl_results)) {
    export_list$TRB_dist_to_ctrl <- dist_to_ctrl_results$TRB
  }
  if ("TRA" %in% names(dist_to_ctrl_results)) {
    export_list$TRA_dist_to_ctrl <- dist_to_ctrl_results$TRA
  }
  if ("TRB_stats" %in% names(dist_to_ctrl_results)) {
    export_list$TRB_vs_ctrl_stats <- dist_to_ctrl_results$TRB_stats
  }
  if ("TRA_stats" %in% names(dist_to_ctrl_results)) {
    export_list$TRA_vs_ctrl_stats <- dist_to_ctrl_results$TRA_stats
  }
  if ("TRB" %in% names(centroid_results)) {
    export_list$TRB_dx_centroids <- as.data.frame(centroid_results$TRB) |> tibble::rownames_to_column("diagnosis")
  }
  if ("TRA" %in% names(centroid_results)) {
    export_list$TRA_dx_centroids <- as.data.frame(centroid_results$TRA) |> tibble::rownames_to_column("diagnosis")
  }

  if (length(export_list) > 0) {
    writexl::write_xlsx(export_list, "results/tcr_comparison/tcrdist/tcrdist_analysis_results.xlsx")
  }

  list(
    trb_dist = trb_sample_dist,
    tra_dist = tra_sample_dist,
    permanova = tcrdist_permanova,
    dist_to_ctrl = dist_to_ctrl_results,
    centroid_dist = centroid_results
  )

}, error = function(e) {

  message("  TCRdist analysis failed: ", e$message)
  NULL
})

# ============================================================================
# Section 10: GLIPH2 specificity group analysis
# ============================================================================
message("Section 10: GLIPH2 specificity group analysis")

dir.create("results/tcr_comparison/gliph", showWarnings = FALSE, recursive = TRUE)

gliph_results <- tryCatch({
  message("  Running GLIPH2 clustering")

  # Run GLIPH2 on TRB chains
  gliph_out <- runGLIPH(
    sc_tcr,
    chains = "TRB",
    local_similarities = TRUE,
    global_similarities = TRUE,
    local_method = "fisher",
    motif_length = 3,
    vgene_match = FALSE,
    n_cores = 1,
    return_seurat = FALSE
  )

  # turboGliph returns cluster info in $clusters with members column (space-separated CDR3s)
  # Parse the members column to create proper membership mapping
  n_clusters <- if (!is.null(gliph_out$clusters)) nrow(gliph_out$clusters) else 0
  message("    GLIPH2 identified ", n_clusters, " specificity groups")

  if (n_clusters == 0) {
    stop("No GLIPH clusters found")
  }

  # Extract enriched motifs
  gliph_motifs <- tryCatch({
    extractGLIPHmotifs(gliph_out, fdr_threshold = 0.1)
  }, error = function(e) NULL)
  n_motifs <- if (!is.null(gliph_motifs)) nrow(gliph_motifs) else 0
  message("    Found ", n_motifs, " enriched motifs (FDR < 0.1)")

  # Parse cluster membership from gliph_out$clusters$members column
  # Each row has space-separated CDR3 sequences in the 'members' column
  message("  Parsing cluster membership from turboGliph output")

  cluster_membership <- lapply(seq_len(nrow(gliph_out$clusters)), function(i) {
    row <- gliph_out$clusters[i, ]
    members_str <- row$members
    if (is.null(members_str) || is.na(members_str) || members_str == "") {
      return(NULL)
    }
    cdr3_seqs <- strsplit(as.character(members_str), "\\s+")[[1]]
    cdr3_seqs <- cdr3_seqs[cdr3_seqs != ""]
    if (length(cdr3_seqs) == 0) return(NULL)

    data.frame(
      cluster = i,
      CDR3b = cdr3_seqs,
      cluster_tag = row$tag,
      cluster_size = row$cluster_size,
      stringsAsFactors = FALSE
    )
  }) |> dplyr::bind_rows()

  if (is.null(cluster_membership) || nrow(cluster_membership) == 0) {
    stop("No cluster membership data could be parsed")
  }
  message("    Parsed ", nrow(cluster_membership), " CDR3-to-cluster mappings")

  # Extract TRB CDR3 sequences from Seurat metadata for matching
  tcr_data <- immApex::getIR(sc_tcr, chains = "TRB")
  tcr_data <- tcr_data |>
    dplyr::left_join(
      sc_tcr@meta.data |>
        tibble::rownames_to_column("barcode") |>
        dplyr::select(barcode, sample, patient, tissue, diagnosis),
      by = "barcode"
    )

  # Match GLIPH clusters to cells via CDR3 sequence
  message("  Mapping GLIPH clusters to cell metadata")
  cluster_meta <- cluster_membership |>
    dplyr::left_join(
      tcr_data |> dplyr::select(barcode, cdr3_aa, sample, patient, tissue, diagnosis) |> dplyr::distinct(),
      by = c("CDR3b" = "cdr3_aa")
    )

  n_matched <- sum(!is.na(cluster_meta$diagnosis))
  message("    Matched ", n_matched, " / ", nrow(cluster_meta), " sequences to metadata")

  # =====================================================
  # DIAGNOSIS-FOCUSED GLIPH ANALYSIS
  # =====================================================
  message("  Analyzing diagnosis enrichment in GLIPH clusters")

  # Calculate diagnosis distribution per cluster
  cluster_diagnosis <- cluster_meta |>
    dplyr::filter(!is.na(cluster), !is.na(diagnosis)) |>
    dplyr::group_by(cluster, diagnosis) |>
    dplyr::summarize(n_cells = n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = n_cells, values_fill = 0)

  # Total cells per diagnosis (for expected proportions)
  total_by_diagnosis <- cluster_meta |>
    dplyr::filter(!is.na(diagnosis)) |>
    dplyr::count(diagnosis, name = "total")

  # Test for diagnosis enrichment using Fisher's exact test
  # Focus on Disease vs CTRL enrichment
  diagnosis_enrichment <- lapply(unique(cluster_meta$cluster[!is.na(cluster_meta$cluster)]), function(cl) {
    cl_data <- cluster_meta |> dplyr::filter(cluster == cl, !is.na(diagnosis))
    if (nrow(cl_data) < 3) return(NULL)

    # Test each diagnosis vs CTRL specifically
    results <- lapply(setdiff(unique(cl_data$diagnosis), "CTRL"), function(dx) {
      # Cells in cluster: diagnosis of interest vs CTRL
      in_cluster_dx <- sum(cl_data$diagnosis == dx)
      in_cluster_ctrl <- sum(cl_data$diagnosis == "CTRL")

      # Total cells in dataset
      total_dx <- total_by_diagnosis$total[total_by_diagnosis$diagnosis == dx]
      total_ctrl <- total_by_diagnosis$total[total_by_diagnosis$diagnosis == "CTRL"]

      if (is.na(total_dx) || is.na(total_ctrl) || total_dx == 0 || total_ctrl == 0) return(NULL)

      out_cluster_dx <- total_dx - in_cluster_dx
      out_cluster_ctrl <- total_ctrl - in_cluster_ctrl

      if (any(c(in_cluster_dx, in_cluster_ctrl, out_cluster_dx, out_cluster_ctrl) < 0)) return(NULL)

      mat <- matrix(c(in_cluster_dx, in_cluster_ctrl, out_cluster_dx, out_cluster_ctrl), nrow = 2)

      tryCatch({
        test <- fisher.test(mat)
        data.frame(
          cluster = cl,
          diagnosis = dx,
          n_dx_in_cluster = in_cluster_dx,
          n_ctrl_in_cluster = in_cluster_ctrl,
          cluster_size = nrow(cl_data),
          odds_ratio = test$estimate,
          p.value = test$p.value
        )
      }, error = function(e) NULL)
    }) |> dplyr::bind_rows()

    results
  }) |> dplyr::bind_rows()

  if (nrow(diagnosis_enrichment) > 0) {
    diagnosis_enrichment <- diagnosis_enrichment |>
      dplyr::mutate(
        p.adj = p.adjust(p.value, method = "BH"),
        enrichment_direction = ifelse(odds_ratio > 1, "Disease-enriched", "CTRL-enriched"),
        log2_OR = log2(odds_ratio + 0.01)
      ) |>
      dplyr::arrange(p.adj)

    # Report significant enrichments
    sig_enrichments <- diagnosis_enrichment |> dplyr::filter(p.adj < 0.1, odds_ratio > 1)
    if (nrow(sig_enrichments) > 0) {
      message("    Disease-enriched clusters (vs CTRL, FDR < 0.1):")
      for (i in seq_len(min(10, nrow(sig_enrichments)))) {
        row <- sig_enrichments[i, ]
        message(sprintf("      Cluster %s - %s: OR=%.2f, FDR=%.4f (n=%d vs %d CTRL)",
                        row$cluster, row$diagnosis, row$odds_ratio, row$p.adj,
                        row$n_dx_in_cluster, row$n_ctrl_in_cluster))
      }
    }
  }

  # =====================================================
  # GLIPH VISUALIZATIONS
  # =====================================================
  message("  Creating GLIPH visualizations")

  # 1. Summary of disease-enriched vs control-enriched clusters
  if (nrow(diagnosis_enrichment) > 0) {
    # Volcano plot: log2(OR) vs -log10(p.adj) for each diagnosis
    p_volcano <- ggplot(diagnosis_enrichment, aes(x = log2_OR, y = -log10(p.adj + 1e-10))) +
      geom_point(aes(color = diagnosis, size = cluster_size), alpha = 0.6) +
      geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "gray50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = diagnosis_col) +
      scale_size_continuous(range = c(1, 6), name = "Cluster size") +
      facet_wrap(~diagnosis, scales = "free") +
      theme_pub() +
      labs(
        title = "GLIPH cluster enrichment: Disease vs CTRL",
        x = "log2(Odds Ratio vs CTRL)",
        y = "-log10(FDR)"
      )
    ggsave("results/tcr_comparison/gliph/cluster_enrichment_volcano.pdf", p_volcano, width = 14, height = 10)

    # Count significant clusters per diagnosis
    sig_summary <- diagnosis_enrichment |>
      dplyr::filter(p.adj < 0.1) |>
      dplyr::group_by(diagnosis, enrichment_direction) |>
      dplyr::summarize(n_clusters = n(), .groups = "drop")

    if (nrow(sig_summary) > 0) {
      # Consistent diagnosis ordering
      dx_levels_gliph <- intersect(neuropathy_dx, unique(sig_summary$diagnosis))
      sig_summary$diagnosis <- factor(sig_summary$diagnosis, levels = dx_levels_gliph)

      # Chi-squared test: is the distribution of enriched clusters non-random across diagnoses?
      enrichment_table <- sig_summary |>
        dplyr::filter(enrichment_direction == "Disease-enriched")
      chisq_pval <- NA
      if (nrow(enrichment_table) > 1 && sum(enrichment_table$n_clusters) > 0) {
        chisq_test <- tryCatch(chisq.test(enrichment_table$n_clusters), error = function(e) NULL)
        if (!is.null(chisq_test)) chisq_pval <- chisq_test$p.value
      }

      subtitle_text <- "FDR < 0.1"
      if (!is.na(chisq_pval)) {
        subtitle_text <- paste0(subtitle_text, " | Chi-squared test across diagnoses: p = ",
                                 signif(chisq_pval, 3))
      }

      p_sig_count <- ggplot(sig_summary, aes(x = diagnosis, y = n_clusters, fill = enrichment_direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c("Disease-enriched" = "firebrick", "CTRL-enriched" = "steelblue")) +
        theme_pub() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
          title = "Number of significantly enriched GLIPH clusters per diagnosis",
          x = "Diagnosis",
          y = "Number of clusters",
          fill = ""
        )
      ggsave("results/tcr_comparison/gliph/significant_cluster_counts.pdf", p_sig_count, width = 10, height = 6)
    }
  }

  # 2. Top clusters by size - composition by diagnosis
  top_clusters <- cluster_meta |>
    dplyr::filter(!is.na(cluster)) |>
    dplyr::count(cluster, sort = TRUE) |>
    dplyr::slice(1:30) |>
    dplyr::pull(cluster)

  if (length(top_clusters) > 0) {
    cluster_comp <- cluster_meta |>
      dplyr::filter(cluster %in% top_clusters, !is.na(diagnosis)) |>
      dplyr::count(cluster, diagnosis) |>
      dplyr::group_by(cluster) |>
      dplyr::mutate(prop = n / sum(n), total = sum(n)) |>
      dplyr::ungroup() |>
      dplyr::arrange(desc(total))

    # Reorder clusters by size
    cluster_order <- cluster_comp |> dplyr::distinct(cluster, total) |>
      dplyr::arrange(desc(total)) |> dplyr::pull(cluster)
    cluster_comp$cluster <- factor(cluster_comp$cluster, levels = cluster_order)

    p_comp <- ggplot(cluster_comp, aes(x = cluster, y = prop, fill = diagnosis)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = diagnosis_col) +
      theme_pub() +
      labs(
        title = "GLIPH cluster composition by diagnosis",
        x = "Cluster (ordered by size)", y = "Proportion"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
    ggsave("results/tcr_comparison/gliph/cluster_composition.pdf", p_comp, width = 14, height = 6)

    # Same but faceted by tissue
    cluster_comp_tissue <- cluster_meta |>
      dplyr::filter(cluster %in% top_clusters[1:15], !is.na(diagnosis), !is.na(tissue)) |>
      dplyr::count(cluster, diagnosis, tissue) |>
      dplyr::group_by(cluster, tissue) |>
      dplyr::mutate(prop = n / sum(n)) |>
      dplyr::ungroup()

    if (nrow(cluster_comp_tissue) > 0) {
      p_comp_tissue <- ggplot(cluster_comp_tissue,
                              aes(x = factor(cluster), y = prop, fill = diagnosis)) +
        geom_bar(stat = "identity") +
        facet_wrap(~tissue) +
        scale_fill_manual(values = diagnosis_col) +
        theme_pub() +
        labs(
          title = "GLIPH cluster composition by diagnosis and tissue",
          x = "Cluster", y = "Proportion"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
      ggsave("results/tcr_comparison/gliph/cluster_composition_by_tissue.pdf",
             p_comp_tissue, width = 14, height = 8)
    }
  }

  # 3. Disease-specific clusters: show CDR3 motifs
  if (nrow(diagnosis_enrichment) > 0) {
    # Get top disease-enriched clusters for each diagnosis
    top_enriched <- diagnosis_enrichment |>
      dplyr::filter(p.adj < 0.2, odds_ratio > 2) |>
      dplyr::group_by(diagnosis) |>
      dplyr::slice_max(order_by = odds_ratio, n = 5) |>
      dplyr::ungroup()

    if (nrow(top_enriched) > 0) {
      # Get CDR3 sequences from these clusters
      enriched_cdr3 <- cluster_meta |>
        dplyr::filter(cluster %in% top_enriched$cluster) |>
        dplyr::left_join(top_enriched |> dplyr::select(cluster, diagnosis, odds_ratio, p.adj),
                         by = "cluster", suffix = c("", "_enriched")) |>
        dplyr::filter(!is.na(diagnosis_enriched)) |>
        dplyr::group_by(cluster, diagnosis_enriched) |>
        dplyr::summarize(
          n_sequences = n(),
          example_CDR3s = paste(head(unique(CDR3b), 5), collapse = ", "),
          odds_ratio = first(odds_ratio),
          p.adj = first(p.adj),
          .groups = "drop"
        )

      # Export disease-specific CDR3 motifs
      if (nrow(enriched_cdr3) > 0) {
        writexl::write_xlsx(
          enriched_cdr3,
          "results/tcr_comparison/gliph/disease_enriched_cdr3_motifs.xlsx"
        )
        message("    Exported disease-enriched CDR3 motifs")
      }
    }
  }

  # 3b. Motif-level analysis: what CDR3 motifs drive diagnosis differences?
  message("    Extracting driving motifs from GLIPH clusters")

  # Extract motifs (k-mers of length 3-4) from CDR3 sequences in enriched clusters
  extract_cdr3_motifs <- function(cdr3_seqs, k = 3) {
    motifs <- unlist(lapply(cdr3_seqs, function(seq) {
      if (is.na(seq) || nchar(seq) < k) return(character(0))
      sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
    }))
    table(motifs)
  }

  if (nrow(diagnosis_enrichment) > 0) {
    # For each diagnosis, compare motifs in disease-enriched vs CTRL-enriched clusters
    dx_with_enrichment <- unique(diagnosis_enrichment$diagnosis)

    motif_driving <- lapply(dx_with_enrichment, function(dx) {
      # Disease-enriched clusters for this diagnosis
      dx_enriched_cl <- diagnosis_enrichment |>
        dplyr::filter(diagnosis == dx, odds_ratio > 1, p.adj < 0.2) |>
        dplyr::pull(cluster)
      # CTRL-enriched clusters
      ctrl_enriched_cl <- diagnosis_enrichment |>
        dplyr::filter(diagnosis == dx, odds_ratio < 1, p.adj < 0.2) |>
        dplyr::pull(cluster)

      # Get CDR3 sequences from each set
      dx_cdr3 <- cluster_meta |>
        dplyr::filter(cluster %in% dx_enriched_cl, !is.na(CDR3b)) |>
        dplyr::pull(CDR3b) |> unique()
      ctrl_cdr3 <- cluster_meta |>
        dplyr::filter(cluster %in% ctrl_enriched_cl, !is.na(CDR3b)) |>
        dplyr::pull(CDR3b) |> unique()

      if (length(dx_cdr3) < 3) return(NULL)

      # Extract 3-mer and 4-mer motifs
      lapply(c(3, 4), function(k) {
        dx_motifs <- extract_cdr3_motifs(dx_cdr3, k)
        # Background: all CDR3 sequences
        all_cdr3 <- cluster_meta |> dplyr::filter(!is.na(CDR3b)) |> dplyr::pull(CDR3b) |> unique()
        bg_motifs <- extract_cdr3_motifs(all_cdr3, k)

        # Compute enrichment for each motif (Fisher's exact)
        shared_motifs <- intersect(names(dx_motifs), names(bg_motifs))
        if (length(shared_motifs) == 0) return(NULL)

        motif_stats <- lapply(shared_motifs, function(m) {
          a <- as.integer(dx_motifs[m])        # motif in disease-enriched
          b <- sum(dx_motifs) - a               # other motifs in disease-enriched
          c_val <- as.integer(bg_motifs[m])     # motif in background
          d <- sum(bg_motifs) - c_val           # other in background
          mat <- matrix(c(a, c_val, b, d), nrow = 2)
          ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
          if (is.null(ft)) return(NULL)
          data.frame(
            diagnosis = dx, motif = m, k = k,
            count_enriched = a, count_background = c_val,
            freq_enriched = a / sum(dx_motifs),
            freq_background = c_val / sum(bg_motifs),
            odds_ratio = ft$estimate, p.value = ft$p.value
          )
        }) |> dplyr::bind_rows()

        if (nrow(motif_stats) > 0) {
          motif_stats$p.adj <- p.adjust(motif_stats$p.value, method = "BH")
          motif_stats$log2_fc <- log2((motif_stats$freq_enriched + 1e-6) /
                                       (motif_stats$freq_background + 1e-6))
        }
        motif_stats
      }) |> dplyr::bind_rows()
    }) |> dplyr::bind_rows()

    if (!is.null(motif_driving) && nrow(motif_driving) > 0) {
      # Export full motif table
      writexl::write_xlsx(
        motif_driving |> dplyr::arrange(p.adj),
        "results/tcr_comparison/gliph/motif_enrichment_by_diagnosis.xlsx"
      )

      # Visualize top driving motifs per diagnosis
      top_motifs <- motif_driving |>
        dplyr::filter(p.adj < 0.1, odds_ratio > 1, k == 3) |>
        dplyr::group_by(diagnosis) |>
        dplyr::slice_max(order_by = log2_fc, n = 10) |>
        dplyr::ungroup()

      if (nrow(top_motifs) > 0) {
        # Consistent diagnosis ordering
        top_motifs$diagnosis <- factor(top_motifs$diagnosis,
                                       levels = intersect(neuropathy_dx, unique(top_motifs$diagnosis)))

        p_motifs <- ggplot(top_motifs, aes(x = reorder(motif, log2_fc), y = log2_fc, fill = diagnosis)) +
          geom_bar(stat = "identity") +
          geom_text(aes(label = sapply(p.adj, function(p) {
            if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
          })), hjust = -0.2, size = 3) +
          coord_flip() +
          facet_wrap(~diagnosis, scales = "free_y") +
          scale_fill_manual(values = diagnosis_col) +
          theme_pub() +
          theme(legend.position = "none") +
          labs(
            title = "Top CDR3 3-mer motifs enriched in disease-associated GLIPH clusters",
            x = "CDR3 motif",
            y = "log2(fold change vs background)"
          )
        ggsave("results/tcr_comparison/gliph/driving_motifs_by_diagnosis.pdf",
               p_motifs, width = 14, height = max(6, length(unique(top_motifs$diagnosis)) * 3))
      }

      message("    Exported driving motif analysis")
    }
  }

  # 3c. Tissue-specific motif and cluster analysis
  message("    Analyzing tissue-specificity of GLIPH clusters")

  # Test each cluster for tissue enrichment (CSF vs PBMC) using Fisher's exact
  tissue_enrichment <- lapply(unique(cluster_meta$cluster[!is.na(cluster_meta$cluster)]), function(cl) {
    cl_data <- cluster_meta |> dplyr::filter(cluster == cl, !is.na(tissue))
    if (nrow(cl_data) < 3) return(NULL)

    n_csf <- sum(cl_data$tissue == "CSF")
    n_pbmc <- sum(cl_data$tissue == "PBMC")

    # Background tissue proportions
    total_csf <- sum(cluster_meta$tissue == "CSF", na.rm = TRUE)
    total_pbmc <- sum(cluster_meta$tissue == "PBMC", na.rm = TRUE)

    mat <- matrix(c(n_csf, total_csf - n_csf, n_pbmc, total_pbmc - n_pbmc), nrow = 2)
    ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
    if (is.null(ft)) return(NULL)

    data.frame(
      cluster = cl, n_csf = n_csf, n_pbmc = n_pbmc,
      csf_fraction = n_csf / (n_csf + n_pbmc),
      tissue_OR = ft$estimate, tissue_p = ft$p.value
    )
  }) |> dplyr::bind_rows()

  if (nrow(tissue_enrichment) > 0) {
    tissue_enrichment$tissue_padj <- p.adjust(tissue_enrichment$tissue_p, method = "BH")
    tissue_enrichment <- tissue_enrichment |>
      dplyr::mutate(tissue_bias = case_when(
        tissue_padj < 0.1 & tissue_OR > 1 ~ "CSF-enriched",
        tissue_padj < 0.1 & tissue_OR < 1 ~ "PBMC-enriched",
        TRUE ~ "No bias"
      ))

    # Combine tissue and diagnosis enrichment
    if (nrow(diagnosis_enrichment) > 0) {
      cluster_tissue_dx <- diagnosis_enrichment |>
        dplyr::filter(p.adj < 0.2, odds_ratio > 1) |>
        dplyr::left_join(tissue_enrichment |> dplyr::select(cluster, csf_fraction, tissue_OR,
                                                              tissue_padj, tissue_bias),
                         by = "cluster")

      if (nrow(cluster_tissue_dx) > 0) {
        # Consistent diagnosis ordering
        cluster_tissue_dx$diagnosis <- factor(cluster_tissue_dx$diagnosis,
                                              levels = intersect(neuropathy_dx, unique(cluster_tissue_dx$diagnosis)))

        # Visualize: scatter of disease enrichment OR vs CSF fraction, colored by tissue bias
        p_tissue_dx <- ggplot(cluster_tissue_dx,
                               aes(x = csf_fraction, y = log2(odds_ratio + 0.01),
                                   color = diagnosis, shape = tissue_bias)) +
          geom_point(aes(size = cluster_size), alpha = 0.7) +
          geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
          scale_color_manual(values = diagnosis_col) +
          scale_shape_manual(values = c("CSF-enriched" = 17, "PBMC-enriched" = 15, "No bias" = 16)) +
          scale_size_continuous(range = c(2, 7), name = "Cluster size") +
          theme_pub() +
          labs(
            title = "Disease-enriched GLIPH clusters: tissue specificity",
            x = "CSF fraction of cluster",
            y = "log2(OR vs CTRL)",
            shape = "Tissue bias"
          )
        ggsave("results/tcr_comparison/gliph/cluster_tissue_vs_diagnosis.pdf",
               p_tissue_dx, width = 11, height = 8)

        # Summary: how many disease-enriched clusters are tissue-biased?
        tissue_bias_summary <- cluster_tissue_dx |>
          dplyr::count(diagnosis, tissue_bias) |>
          tidyr::pivot_wider(names_from = tissue_bias, values_from = n, values_fill = 0)

        # Bar chart of tissue bias per diagnosis
        tissue_bias_long <- cluster_tissue_dx |>
          dplyr::count(diagnosis, tissue_bias)

        tissue_bias_long$tissue_bias <- factor(tissue_bias_long$tissue_bias,
                                                levels = c("CSF-enriched", "No bias", "PBMC-enriched"))

        p_tissue_bias <- ggplot(tissue_bias_long,
                                 aes(x = diagnosis, y = n, fill = tissue_bias)) +
          geom_bar(stat = "identity", position = "stack") +
          scale_fill_manual(values = c("CSF-enriched" = tissue_col["CSF"],
                                        "PBMC-enriched" = tissue_col["PBMC"],
                                        "No bias" = "gray70")) +
          theme_pub() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(
            title = "Tissue bias of disease-enriched GLIPH clusters",
            x = "Diagnosis",
            y = "Number of clusters",
            fill = "Tissue bias"
          )
        ggsave("results/tcr_comparison/gliph/tissue_bias_per_diagnosis.pdf",
               p_tissue_bias, width = 10, height = 6)

        # Fisher's exact: is tissue bias associated with diagnosis?
        if (nrow(tissue_bias_summary) > 1) {
          tissue_dx_table <- cluster_tissue_dx |>
            dplyr::mutate(is_csf_biased = tissue_bias == "CSF-enriched") |>
            dplyr::count(diagnosis, is_csf_biased) |>
            tidyr::pivot_wider(names_from = is_csf_biased, values_from = n, values_fill = 0)
          if (ncol(tissue_dx_table) == 3 && nrow(tissue_dx_table) > 1) {
            fisher_tissue <- tryCatch(
              fisher.test(as.matrix(tissue_dx_table[, -1])),
              error = function(e) NULL
            )
            if (!is.null(fisher_tissue)) {
              message(sprintf("      Fisher's test for tissue bias ~ diagnosis: p = %.4f",
                              fisher_tissue$p.value))
            }
          }
        }

        # Export tissue-specificity results
        writexl::write_xlsx(
          list(
            cluster_tissue_dx = cluster_tissue_dx,
            tissue_enrichment = tissue_enrichment,
            tissue_bias_summary = tissue_bias_summary
          ),
          "results/tcr_comparison/gliph/tissue_specificity_analysis.xlsx"
        )
        message("    Exported tissue-specificity analysis")
      }
    }
  }

  # 3d. Tissue-specific motifs: which motifs are enriched in CSF vs PBMC within disease clusters?
  message("    Analyzing tissue-specific motifs")

  if (nrow(tissue_enrichment) > 0) {
    # Get CDR3 sequences from CSF-enriched vs PBMC-enriched disease clusters
    csf_clusters <- tissue_enrichment$cluster[tissue_enrichment$tissue_bias == "CSF-enriched"]
    pbmc_clusters <- tissue_enrichment$cluster[tissue_enrichment$tissue_bias == "PBMC-enriched"]

    csf_cdr3 <- cluster_meta |>
      dplyr::filter(cluster %in% csf_clusters, !is.na(CDR3b)) |>
      dplyr::pull(CDR3b) |> unique()
    pbmc_cdr3 <- cluster_meta |>
      dplyr::filter(cluster %in% pbmc_clusters, !is.na(CDR3b)) |>
      dplyr::pull(CDR3b) |> unique()

    if (length(csf_cdr3) >= 5 && length(pbmc_cdr3) >= 5) {
      tissue_motif_comparison <- lapply(c(3, 4), function(k) {
        csf_m <- extract_cdr3_motifs(csf_cdr3, k)
        pbmc_m <- extract_cdr3_motifs(pbmc_cdr3, k)
        shared <- intersect(names(csf_m), names(pbmc_m))
        if (length(shared) == 0) return(NULL)

        lapply(shared, function(m) {
          a <- as.integer(csf_m[m]); b <- sum(csf_m) - a
          c_val <- as.integer(pbmc_m[m]); d <- sum(pbmc_m) - c_val
          ft <- tryCatch(fisher.test(matrix(c(a, c_val, b, d), nrow = 2)), error = function(e) NULL)
          if (is.null(ft)) return(NULL)
          data.frame(
            motif = m, k = k,
            count_csf = a, count_pbmc = c_val,
            freq_csf = a / sum(csf_m), freq_pbmc = c_val / sum(pbmc_m),
            odds_ratio = ft$estimate, p.value = ft$p.value
          )
        }) |> dplyr::bind_rows()
      }) |> dplyr::bind_rows()

      if (nrow(tissue_motif_comparison) > 0) {
        tissue_motif_comparison$p.adj <- p.adjust(tissue_motif_comparison$p.value, method = "BH")
        tissue_motif_comparison$log2_fc <- log2((tissue_motif_comparison$freq_csf + 1e-6) /
                                                  (tissue_motif_comparison$freq_pbmc + 1e-6))

        # Top tissue-biased motifs
        top_tissue_motifs <- tissue_motif_comparison |>
          dplyr::filter(p.adj < 0.1, k == 3) |>
          dplyr::arrange(desc(abs(log2_fc))) |>
          dplyr::slice(1:min(30, n()))

        if (nrow(top_tissue_motifs) > 0) {
          top_tissue_motifs$tissue_direction <- ifelse(top_tissue_motifs$log2_fc > 0,
                                                        "CSF-enriched", "PBMC-enriched")

          p_tissue_motifs <- ggplot(top_tissue_motifs,
                                     aes(x = reorder(motif, log2_fc), y = log2_fc,
                                         fill = tissue_direction)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = sapply(p.adj, function(p) {
              if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
            })), hjust = ifelse(top_tissue_motifs$log2_fc > 0, -0.2, 1.2), size = 3) +
            coord_flip() +
            scale_fill_manual(values = c("CSF-enriched" = tissue_col["CSF"],
                                          "PBMC-enriched" = tissue_col["PBMC"])) +
            theme_pub() +
            labs(
              title = "CDR3 motifs differentially enriched in CSF vs PBMC GLIPH clusters",
              subtitle = "3-mer motifs (FDR < 0.1)",
              x = "CDR3 motif",
              y = "log2(CSF freq / PBMC freq)",
              fill = ""
            )
          ggsave("results/tcr_comparison/gliph/tissue_specific_motifs.pdf",
                 p_tissue_motifs, width = 10, height = max(6, nrow(top_tissue_motifs) * 0.3))
        }

        writexl::write_xlsx(
          tissue_motif_comparison |> dplyr::arrange(p.adj),
          "results/tcr_comparison/gliph/tissue_specific_motif_comparison.xlsx"
        )
        message("    Exported tissue-specific motif comparison")
      }
    } else {
      message("    Insufficient tissue-biased clusters for motif comparison (CSF: ",
              length(csf_cdr3), ", PBMC: ", length(pbmc_cdr3), " unique CDR3s)")
    }
  }

 
  # 4. Per-diagnosis private cluster analysis
  disease_specific <- cluster_sharing |>
    dplyr::filter(cluster_type == "Disease-specific") |>
    dplyr::separate_rows(diagnoses, sep = ",") |>
    dplyr::count(diagnoses, name = "n_private_clusters")

  if (nrow(disease_specific) > 0) {
    # Consistent diagnosis ordering
    ds_levels <- intersect(c("CTRL", neuropathy_dx), unique(disease_specific$diagnoses))
    disease_specific$diagnoses <- factor(disease_specific$diagnoses, levels = ds_levels)

    # Binomial test: does any diagnosis have more private clusters than expected by chance?
    # Expected proportion based on number of cells per diagnosis
    total_private <- sum(disease_specific$n_private_clusters)
    if (total_private > 0 && nrow(total_by_diagnosis) > 0) {
      disease_specific <- disease_specific |>
        dplyr::left_join(total_by_diagnosis, by = c("diagnoses" = "diagnosis")) |>
        dplyr::mutate(
          expected_prop = total / sum(total_by_diagnosis$total),
          binom_p = mapply(function(k, p) {
            tryCatch(binom.test(k, total_private, p, alternative = "greater")$p.value,
                     error = function(e) NA)
          }, n_private_clusters, expected_prop)
        )
      disease_specific$binom_padj <- p.adjust(disease_specific$binom_p, method = "BH")
    }

    p_private <- ggplot(disease_specific, aes(x = diagnoses, y = n_private_clusters, fill = diagnoses)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = diagnosis_col) +
      theme_pub() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

    # Annotate bars that are significantly enriched
    if ("binom_padj" %in% colnames(disease_specific)) {
      sig_private <- disease_specific |> dplyr::filter(binom_padj < 0.05)
      if (nrow(sig_private) > 0) {
        p_private <- p_private +
          geom_text(data = sig_private,
                    aes(label = sapply(binom_padj, format_pval)),
                    vjust = -0.5, size = 4)
      }
    }

    p_private <- p_private +
      labs(
        title = "Disease-specific GLIPH clusters",
        subtitle = "Clusters found only in one diagnosis group (binomial test vs expected proportion)",
        x = "Diagnosis",
        y = "Number of private clusters"
      )
    ggsave("results/tcr_comparison/gliph/disease_specific_clusters.pdf", p_private, width = 8, height = 6)
  }

  # 6. Heatmap of cluster presence across patient-tissue combinations
  # This helps visualize patient effects
  patient_cluster_mat <- cluster_meta |>
    dplyr::filter(!is.na(cluster), !is.na(patient), !is.na(tissue)) |>
    dplyr::mutate(patient_tissue = paste(patient, tissue, sep = "_")) |>
    dplyr::count(cluster, patient_tissue) |>
    tidyr::pivot_wider(names_from = patient_tissue, values_from = n, values_fill = 0) |>
    tibble::column_to_rownames("cluster") |>
    as.matrix()

  if (nrow(patient_cluster_mat) > 5 && ncol(patient_cluster_mat) > 2) {
    # Select top clusters by variance
    cluster_var <- apply(patient_cluster_mat, 1, var)
    top_var_clusters <- names(sort(cluster_var, decreasing = TRUE))[1:min(50, nrow(patient_cluster_mat))]
    mat_subset <- patient_cluster_mat[top_var_clusters, ]

    # Create annotation for columns
    col_anno <- data.frame(
      patient_tissue = colnames(mat_subset)
    ) |>
      tidyr::separate(patient_tissue, into = c("patient", "tissue"), sep = "_", remove = FALSE) |>
      dplyr::left_join(sample_summary |> dplyr::distinct(patient, diagnosis), by = "patient") |>
      tibble::column_to_rownames("patient_tissue")

    pdf("results/tcr_comparison/gliph/cluster_patient_tissue_heatmap.pdf", width = 14, height = 12)
    pheatmap::pheatmap(
      log10(mat_subset + 1),
      color = colorRampPalette(c("white", "navy"))(50),
      annotation_col = col_anno |> dplyr::select(diagnosis, tissue),
      annotation_colors = list(diagnosis = diagnosis_col, tissue = tissue_col),
      main = "GLIPH cluster presence across patients/tissues",
      fontsize_row = 6,
      fontsize_col = 8,
      cluster_cols = TRUE,
      cluster_rows = TRUE
    )
    dev.off()
  }

  # Export comprehensive results
  message("  Exporting GLIPH results")

  export_list <- list()
  if (!is.null(gliph_out$clusters)) export_list$cluster_properties <- gliph_out$clusters
  if (!is.null(cluster_meta)) export_list$cluster_membership <- cluster_meta
  if (nrow(diagnosis_enrichment) > 0) export_list$diagnosis_enrichment <- diagnosis_enrichment
  if (!is.null(gliph_motifs) && nrow(gliph_motifs) > 0) export_list$motifs <- gliph_motifs
  export_list$cluster_sharing <- cluster_sharing
  if (exists("tissue_enrichment") && nrow(tissue_enrichment) > 0) {
    export_list$tissue_enrichment <- tissue_enrichment
  }
  if (exists("motif_driving") && !is.null(motif_driving) && nrow(motif_driving) > 0) {
    export_list$motif_driving <- motif_driving
  }

  writexl::write_xlsx(export_list, "results/tcr_comparison/gliph/gliph_analysis_results.xlsx")

  message("  GLIPH analysis COMPLETE!")

  list(
    gliph = gliph_out,
    cluster_meta = cluster_meta,
    diagnosis_enrichment = diagnosis_enrichment,
    cluster_sharing = cluster_sharing,
    motifs = gliph_motifs,
    tissue_enrichment = if (exists("tissue_enrichment")) tissue_enrichment else NULL,
    motif_driving = if (exists("motif_driving")) motif_driving else NULL
  )

}, error = function(e) {
  message("  GLIPH2 analysis failed: ", e$message)
  NULL
})


# ============================================================================
# Section 11: Summary export
# ============================================================================
message("Section 11: Exporting summary")

all_stats <- list(
  aa_composition = aa_stats_combined,
  gene_usage = gene_usage_combined,
  cdr3_length = length_stats_combined,
  physicochemical = physico_stats,
  motif_enrichment = motif_enrichment_combined
)

if (exists("diversity_stats_dx_combined")) all_stats$diversity_by_diagnosis <- diversity_stats_dx_combined
if (exists("diversity_stats_tissue")) all_stats$diversity_by_tissue <- diversity_stats_tissue
if (exists("diversity_interaction")) all_stats$diversity_interaction <- diversity_interaction

# Add AA composition PERMANOVA results
permanova_summary <- lapply(names(permanova_results), function(ch) {
  res <- permanova_results[[ch]]
  data.frame(
    analysis = "aa_composition",
    chain = ch,
    term = rownames(res),
    Df = res$Df,
    SumOfSqs = res$SumOfSqs,
    R2 = res$R2,
    F = res$F,
    p.value = res[, "Pr(>F)"]
  )
}) |> dplyr::bind_rows()

# Add ESM embedding PERMANOVA results
if (exists("esm_permanova_results") && length(esm_permanova_results) > 0) {
  esm_summary <- lapply(names(esm_permanova_results), function(red) {
    res <- esm_permanova_results[[red]]
    data.frame(
      analysis = "esm_embeddings",
      chain = gsub("esm_", "", red),
      term = rownames(res),
      Df = res$Df,
      SumOfSqs = res$SumOfSqs,
      R2 = res$R2,
      F = res$F,
      p.value = res[, "Pr(>F)"]
    )
  }) |> dplyr::bind_rows()
  permanova_summary <- dplyr::bind_rows(permanova_summary, esm_summary)
}

# Add TCRdist PERMANOVA results
if (!is.null(tcrdist_results) && length(tcrdist_results$permanova) > 0) {
  tcrdist_summary <- lapply(names(tcrdist_results$permanova), function(ch) {
    res <- tcrdist_results$permanova[[ch]]
    data.frame(
      analysis = "tcrdist",
      chain = ch,
      term = rownames(res),
      Df = res$Df,
      SumOfSqs = res$SumOfSqs,
      R2 = res$R2,
      F = res$F,
      p.value = res[, "Pr(>F)"]
    )
  }) |> dplyr::bind_rows()
  permanova_summary <- dplyr::bind_rows(permanova_summary, tcrdist_summary)
}

if (nrow(permanova_summary) > 0) {
  all_stats$permanova <- permanova_summary
}

# Add GLIPH diagnosis enrichment results
if (!is.null(gliph_results) && !is.null(gliph_results$enrichment) && nrow(gliph_results$enrichment) > 0) {
  all_stats$gliph_diagnosis_enrichment <- gliph_results$enrichment
}
if (!is.null(gliph_results) && !is.null(gliph_results$tissue_enrichment) && nrow(gliph_results$tissue_enrichment) > 0) {
  all_stats$gliph_tissue_enrichment <- gliph_results$tissue_enrichment
}
if (!is.null(gliph_results) && !is.null(gliph_results$motif_driving) && nrow(gliph_results$motif_driving) > 0) {
  all_stats$gliph_driving_motifs <- gliph_results$motif_driving
}

writexl::write_xlsx(all_stats, "results/tcr_comparison/tables/tcr_statistics_comprehensive.xlsx")

message("Analysis complete. Results in results/tcr_comparison/")
