# TCR repertoire comparison across neuropathy diagnoses
# Compares TRA, TRB, TRAB chains with mixed models (patient random effect, tissue fixed effect)

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
    if (!reduction_name %in% names(sc_tcr@reductions)) {
      sc_tcr <- immLynx::runEmbeddings(
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
    
    # B. STATISTICAL: PERMANOVA (Multivariate ANOVA)
    set.seed(42)
    #idx_sample <- sample(nrow(embeddings), min(2000, nrow(embeddings))) 
    
    dist_mat <- dist(embeddings)
    permanova <- adonis2(dist_mat ~ diagnosis, 
                         data = meta_sub, 
                         permutations = 100)
    
    message(paste("    PERMANOVA R2:", round(permanova$R2[1], 4), " P-val:", permanova$`Pr(>F)`[1]))
    capture.output(permanova, file = paste0("results/tcr_comparison/embeddings/stats_", red_name, ".txt"))
    
    # C. VISUALIZATION: UMAP with Density
    set.seed(42)
    umap_res <- uwot::umap(embeddings, n_neighbors = 15, min_dist = 0.1)
    
    plot_df <- data.frame(
      UMAP1 = umap_res[,1],
      UMAP2 = umap_res[,2],
      diagnosis = meta_sub$diagnosis
    )
    
    p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = diagnosis)) +
      geom_point(alpha = 0.4, size = 0.5) +
      stat_density_2d(aes(alpha = ..level..), bins = 5, color = "black", size = 0.2) + # Adds contours
      scale_color_manual(values = diagnosis_col) +
      theme_pub() +
      labs(
        title = paste("ESM2 UMAP:", toupper(gsub("esm_", "", red_name))),
        subtitle = paste("PERMANOVA R2 =", round(permanova$R2[1], 3))
      ) +
      facet_wrap(~diagnosis) # Split view to see distinct clusters
    
    ggsave(paste0("results/tcr_comparison/embeddings/umap_density_", red_name, ".pdf"), p, width = 10, height = 8)
  }
  
  message("  Embeddings analysis complete")
  
}, error = function(e) {
  message("  Embeddings failed: ", e$message)
})

# ============================================================================
# Section 9: Summary export
# ============================================================================
message("Section 9: Exporting summary")

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

# Add PERMANOVA results
permanova_summary <- lapply(names(permanova_results), function(ch) {
  res <- permanova_results[[ch]]
  data.frame(
    chain = ch,
    term = rownames(res),
    Df = res$Df,
    SumOfSqs = res$SumOfSqs,
    R2 = res$R2,
    F = res$F,
    p.value = res[, "Pr(>F)"]
  )
}) |> dplyr::bind_rows()

if (nrow(permanova_summary) > 0) {
  all_stats$permanova <- permanova_summary
}

writexl::write_xlsx(all_stats, "results/tcr_comparison/tables/tcr_statistics_comprehensive.xlsx")

message("Analysis complete. Results in results/tcr_comparison/")
