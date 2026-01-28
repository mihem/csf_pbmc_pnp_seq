##################################################
# TCR Featurization and Comparative Analysis
# Comprehensive comparison across metadata variables
##################################################

# ============== SECTION 1: SETUP ==============

# Load libraries ----
library(qs)
library(Seurat)
library(tidyverse)
library(scRepertoire)
library(writexl)
library(readxl)
library(pheatmap)
library(patchwork)
library(ggpubr)
library(Peptides)
library(viridis)
library(RColorBrewer)
library(uwot)
library(cluster)
library(immApex)
library(ggplot2)
library(vegan)
library(ComplexHeatmap)
library(umap)

devtools::load_all("/Users/nick/Documents/GitHub/immLynx")

# Create output directories ----
output_dirs <- c(
    "results/tcr_comparison",
    "results/tcr_comparison/gene_usage",
    "results/tcr_comparison/cdr3_features",
    "results/tcr_comparison/physicochemical",
    "results/tcr_comparison/clonality",
    "results/tcr_comparison/similarity",
    "results/tcr_comparison/embeddings",
    "results/tcr_comparison/tables"
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Load data ----
sc_tcr <- qs::qread(file.path("sc_tcr.qs"), nthreads = 6)

# Extract color palettes from sc_tcr@misc ----
diagnosis_col <- sc_tcr@misc$diagnosis_col
tissue_diagnosis_col <- sc_tcr@misc$tissue_diagnosis_col
cluster_col <- sc_tcr@misc$cluster_col
tissue_col <- c("CSF" = "#E41A1C", "PBMC" = "#377EB8")
group_col <- c("CTRL" = "#4DAF4A", "PNP" = "#984EA3")

# Define comparison variable groups ----
comparison_vars <- list(
    tissue = list(
        var = "tissue",
        levels = c("CSF", "PBMC"),
        colors = tissue_col,
        paired = TRUE,
        paired_by = "patient"
    ),
    group = list(
        var = "group",
        levels = c("CTRL", "PNP"),
        colors = group_col,
        paired = FALSE
    ),
    diagnosis = list(
        var = "diagnosis",
        colors = diagnosis_col,
        paired = FALSE
    ),
    tissue_group = list(
        var = "tissue_group",
        colors = c(
            "CSF_CTRL" = "#66C2A5", "CSF_PNP" = "#FC8D62",
            "PBMC_CTRL" = "#8DA0CB", "PBMC_PNP" = "#E78AC3"
        ),
        paired = FALSE
    ),
    tissue_diagnosis = list(
        var = "tissue_diagnosis",
        colors = tissue_diagnosis_col,
        paired = FALSE
    ),
    cluster = list(
        var = "cluster",
        colors = cluster_col,
        paired = FALSE
    )
)

# Extract TCR metadata ----
tcr_metadata <- sc_tcr@meta.data |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::filter(!is.na(CTaa)) |>
    dplyr::mutate(
        TRA_CDR3 = gsub("_.*", "", CTaa),
        TRB_CDR3 = gsub(".*_", "", CTaa),
        TRA_length = nchar(TRA_CDR3),
        TRB_length = nchar(TRB_CDR3)
    )

# Get unique clonotypes per sample for sample-level analysis ----
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

# Helper functions ----

#' Compare groups with appropriate statistical test
#' @param data Data frame
#' @param value_col Column name for values
#' @param group_col Column name for grouping
#' @param paired Whether to use paired test
#' @param paired_by Column for pairing (patient)
compare_groups <- function(data, value_col, group_col, paired = FALSE, paired_by = NULL) {
  # Remove NAs
  data <- data |> dplyr::filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]]))
  
  groups <- unique(data[[group_col]])
  
  if (length(groups) != 2) {
    return(list(statistic = NA, p.value = NA, method = "Insufficient groups"))
  }
  
  if (paired && !is.null(paired_by)) {
    # For paired test, keep only subjects with data in BOTH groups
    paired_data <- data |>
      dplyr::group_by(.data[[paired_by]]) |>
      dplyr::filter(dplyr::n_distinct(.data[[group_col]]) == 2) |>
      dplyr::ungroup() |>
      dplyr::arrange(.data[[paired_by]], .data[[group_col]])
    
    if (nrow(paired_data) == 0) {
      return(list(statistic = NA, p.value = NA, method = "No complete pairs"))
    }
    
    # Ensure equal length by taking one observation per subject per group
    paired_data <- paired_data |>
      dplyr::group_by(.data[[paired_by]], .data[[group_col]]) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::arrange(.data[[paired_by]], .data[[group_col]])
    
    x <- paired_data[[value_col]][paired_data[[group_col]] == groups[1]]
    y <- paired_data[[value_col]][paired_data[[group_col]] == groups[2]]
    
    if (length(x) < 3 || length(y) < 3) {
      return(list(statistic = NA, p.value = NA, method = "Too few pairs"))
    }
    
    result <- wilcox.test(x, y, paired = TRUE)
  } else {
    x <- data[[value_col]][data[[group_col]] == groups[1]]
    y <- data[[value_col]][data[[group_col]] == groups[2]]
    
    if (length(x) < 3 || length(y) < 3) {
      return(list(statistic = NA, p.value = NA, method = "Too few observations"))
    }
    
    result <- wilcox.test(x, y, paired = FALSE)
  }
  
  return(result)
}


#' Create boxplot with statistical annotation
create_comparison_boxplot <- function(data, x_var, y_var, fill_var = NULL,
                                      colors = NULL, title = "", paired = FALSE) {
    fill_var <- fill_var %||% x_var

    p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[fill_var]])) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
        stat_compare_means(method = ifelse(paired, "wilcox.test", "kruskal.test"),
                          label = "p.format") +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
        ) +
        labs(title = title, x = "", y = y_var)

    if (!is.null(colors)) {
        p <- p + scale_fill_manual(values = colors)
    }

    return(p)
}


# ============== SECTION 2: V/J GENE USAGE ANALYSIS ==============

message("Starting V/J Gene Usage Analysis...")

# 2.1 Gene usage quantification per group ----

# Calculate gene usage for different groupings
gene_usage_results <- list()

for (gene_type in c("Vgene", "Jgene")) {
    for (chain in c("TRA", "TRB")) {
        for (comp_name in names(comparison_vars)) {
            comp <- comparison_vars[[comp_name]]

            tryCatch({
                gene_df <- percentGenes(
                    sc_tcr,
                    chain = chain,
                    gene = gene_type,
                    exportTable = TRUE,
                    group.by = comp$var
                )

                gene_usage_results[[paste(chain, gene_type, comp_name, sep = "_")]] <- gene_df
            }, error = function(e) {
                message(paste("Skipping", chain, gene_type, comp_name, ":", e$message))
            })
        }
    }
}

# Export gene usage tables
if (length(gene_usage_results) > 0) {
  # Keep only data frame elements
  df_list <- gene_usage_results[sapply(gene_usage_results, is.data.frame)]
  
  if (length(df_list) > 0) {
    writexl::write_xlsx(
      df_list,
      file.path("results", "tcr_comparison", "gene_usage", "gene_usage_tables.xlsx")
    )
  }
}

# 2.2 PCA of gene usage ----

# Create PCA plots for TRBV gene usage across different groupings
pca_plots <- list()

for (comp_name in c("tissue_diagnosis", "tissue_group", "diagnosis", "cluster")) {
    comp <- comparison_vars[[comp_name]]
    
    #TRBV Gene Usage
    tryCatch({
        df_genes <- percentGenes(
            sc_tcr,
            chain = "TRB",
            gene = "Vgene",
            exportTable = TRUE,
            group.by = comp$var
        )

        # Perform PCA
        pca_result <- prcomp(t(df_genes), scale. = TRUE, center = TRUE)

        # Create plot data frame
        df_plot <- as.data.frame(pca_result$x[, 1:2])
        df_plot$Sample <- colnames(df_genes)

        # Calculate variance explained
        var_explained <- summary(pca_result)$importance[2, 1:2] * 100

        p <- ggplot(df_plot, aes(x = PC1, y = PC2, fill = Sample)) +
            geom_point(shape = 21, size = 4, stroke = 0.5) +
            scale_fill_manual(values = comp$colors) +
            theme_bw() +
            labs(
                title = paste("TRBV Gene Usage PCA by", comp_name),
                x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
                y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
                fill = comp_name
            ) +
            theme(legend.position = "right")

        pca_plots[[comp_name]] <- p

        ggsave(
            file.path("results", "tcr_comparison", "gene_usage",
                     paste0("TRBVgene_usage_pca_", comp_name, ".pdf")),
            p, width = 8, height = 6
        )
    }, error = function(e) {
        message(paste("PCA failed for", comp_name, ":", e$message))
    })
    
    #TRBV/J Gene Usage
    tryCatch({
      df_genes <- percentVJ(
        sc_tcr,
        chain = "TRB",
        exportTable = TRUE,
        group.by = comp$var
      )
      
      df_genes  <- df_genes [apply(df_genes , 1, var) != 0,]
      
      # Perform PCA
      pca_result <- prcomp(t(df_genes), scale. = TRUE, center = TRUE)
      
      # Create plot data frame
      df_plot <- as.data.frame(pca_result$x[, 1:2])
      df_plot$Sample <- colnames(df_genes)
      
      # Calculate variance explained
      var_explained <- summary(pca_result)$importance[2, 1:2] * 100
      
      p <- ggplot(df_plot, aes(x = PC1, y = PC2, fill = Sample)) +
        geom_point(shape = 21, size = 4, stroke = 0.5) +
        scale_fill_manual(values = comp$colors) +
        theme_bw() +
        labs(
          title = paste("TRBV/J Gene Usage PCA by", comp_name),
          x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
          y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
          fill = comp_name
        ) +
        theme(legend.position = "right")
      
      pca_plots[[comp_name]] <- p
      
      ggsave(
        file.path("results", "tcr_comparison", "gene_usage",
                  paste0("TRBVJgene_usage_pca_", comp_name, ".pdf")),
        p, width = 8, height = 6
      )
    }, error = function(e) {
      message(paste("PCA failed for", comp_name, ":", e$message))
    })
}



# 2.4 Differential gene usage (Fisher's test) ----

# Compare CTRL vs PNP
differential_genes <- list()

tcr_metadata <- tcr_metadata |>
  dplyr::mutate(
    # Split into TRA and TRB chains
    TRA_genes = stringr::str_extract(CTgene, "^[^_]+"),
    TRB_genes = stringr::str_extract(CTgene, "[^_]+$"),
    # Extract V and J genes from each chain
    TRAV = stringr::str_extract(TRA_genes, "TRAV[^.]+"),
    TRAJ = stringr::str_extract(TRA_genes, "TRAJ[^.]+"),
    TRBV = stringr::str_extract(TRB_genes, "TRBV[^.]+"),
    TRBJ = stringr::str_extract(TRB_genes, "TRBJ[^.]+")
  )

for (chain in c("TRA", "TRB")) {
    for (gene_type in c("V", "J")) {
        gene_col <- paste0(chain, gene_type)

        # Get counts per group
        gene_counts <- tcr_metadata |>
            dplyr::filter(!is.na(.data[[gene_col]])) |>
            dplyr::count(.data[[gene_col]], group) |>
            tidyr::pivot_wider(names_from = group, values_from = n, values_fill = 0)

        # Fisher's test for each gene
        fisher_results <- gene_counts |>
            dplyr::rowwise() |>
            dplyr::mutate(
                total_ctrl = sum(tcr_metadata$group == "CTRL", na.rm = TRUE),
                total_pnp = sum(tcr_metadata$group == "PNP", na.rm = TRUE),
                fisher_p = {
                    mat <- matrix(c(CTRL, PNP, total_ctrl - CTRL, total_pnp - PNP),
                                 nrow = 2)
                    fisher.test(mat)$p.value
                },
                odds_ratio = {
                    mat <- matrix(c(CTRL, PNP, total_ctrl - CTRL, total_pnp - PNP),
                                 nrow = 2)
                    fisher.test(mat)$estimate
                }
            ) |>
            dplyr::ungroup() |>
            dplyr::mutate(p_adj = p.adjust(fisher_p, method = "BH")) |>
            dplyr::arrange(p_adj)

        differential_genes[[paste(chain, gene_type, sep = "_")]] <- fisher_results
    }
}

writexl::write_xlsx(
    differential_genes,
    file.path("results", "tcr_comparison", "gene_usage", "differential_gene_usage.xlsx")
)



# ============== SECTION 3: CDR3 SEQUENCE FEATURES ==============

message("Starting CDR3 Sequence Features Analysis...")

# 3.1 CDR3 length distributions ----

length_plots <- list()
length_stats <- list()

for (chain in c("TRA", "TRB")) {
    length_col <- paste0(chain, "_length")

    for (comp_name in c("tissue", "group", "diagnosis", "tissue_diagnosis")) {
        comp <- comparison_vars[[comp_name]]

        # Create plot
        p <- create_comparison_boxplot(
            tcr_clonotypes,
            x_var = comp$var,
            y_var = length_col,
            colors = comp$colors,
            title = paste(chain, "CDR3 Length by", comp_name),
            paired = comp$paired %||% FALSE
        )

        length_plots[[paste(chain, comp_name, sep = "_")]] <- p

        # Statistical test
        test_result <- compare_groups(
            tcr_clonotypes, length_col, comp$var,
            paired = comp$paired %||% FALSE,
            paired_by = comp$paired_by
        )

        length_stats[[paste(chain, comp_name, sep = "_")]] <- data.frame(
            chain = chain,
            comparison = comp_name,
            statistic = test_result$statistic,
            p_value = test_result$p.value,
            method = test_result$method
        )
    }
}

# Combine and save length plots
length_combined <- wrap_plots(length_plots, ncol = 4)
ggsave(
    file.path("results", "tcr_comparison", "cdr3_features", "cdr3_length_comparison.pdf"),
    length_combined, width = 16, height = 12
)

# Save statistics
length_stats_df <- dplyr::bind_rows(length_stats)
writexl::write_xlsx(
    list(length_stats = length_stats_df),
    file.path("results", "tcr_comparison", "cdr3_features", "cdr3_length_stats.xlsx")
)


# Define comparison variables for AA analysis
aa_comparisons <- list(
  group = list(
    var = "group",
    colors = group_col
  ),
  tissue = list(
    var = "tissue",
    colors = tissue_col
  ),
  diagnosis = list(
    var = "diagnosis",
    colors = diagnosis_col
  ),
  tissue_diagnosis = list(
    var = "tissue_diagnosis",
    colors = tissue_diagnosis_col
  ),
  cluster = list(
    var = "seurat_clusters",
    colors = cluster_col 
  )
)

aa_analysis_results <- list()

for (chain in c("TRA", "TRB")) {
  cdr3_col <- paste0(chain, "_CDR3")
  
  # Calculate AA composition
  aa_data <- tcr_metadata |>
    dplyr::filter(!is.na(.data[[cdr3_col]])) |>
    dplyr::mutate(
      aa_comp = purrr::map(.data[[cdr3_col]], calculate_aa_composition)
    ) |>
    tidyr::unnest_wider(aa_comp)
  
  aa_matrix <- aa_data |>
    dplyr::select(all_of(aa_cols)) |>
    as.matrix()
  
  complete_rows <- complete.cases(aa_matrix)
  aa_matrix <- aa_matrix[complete_rows, ]
  aa_data <- aa_data[complete_rows, ]
  
  # --- 1. PCA ---
  pca_result <- prcomp(aa_matrix, scale. = TRUE, center = TRUE)
  var_explained <- summary(pca_result)$importance[2, ] * 100
  
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    tissue = aa_data$tissue,
    group = aa_data$group,
    diagnosis = aa_data$diagnosis,
    tissue_diagnosis = aa_data$tissue_diagnosis,
    seurat_clusters = aa_data$seurat_clusters
  )
  
  # PCA plots for each comparison
  for (comp_name in names(aa_comparisons)) {
    comp <- aa_comparisons[[comp_name]]
    
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[comp$var]])) +
      geom_point(alpha = 0.5, size = 1) +
      labs(
        title = paste(chain, "CDR3 AA Composition - PCA by", comp_name),
        x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
        y = paste0("PC2 (", round(var_explained[2], 1), "%)")
      ) +
      theme_minimal()
    
    # Add ellipse only for comparisons with few groups
    if (comp_name %in% c("group", "tissue", "tissue_diagnosis")) {
      p <- p + stat_ellipse(level = 0.95)
    }
    
    # Add colors if specified
    if (!is.null(comp$colors)) {
      p <- p + scale_color_manual(values = comp$colors)
    }
    
    ggsave(file.path("results", "tcr_comparison", "cdr3_features", 
                     paste0("aa_pca_", chain, "_", comp_name, ".pdf")), p, width = 8, height = 6)
  }
  
  # PCA loadings
  loadings_df <- data.frame(
    AA = aa_cols,
    PC1 = pca_result$rotation[, 1],
    PC2 = pca_result$rotation[, 2]
  )
  
  p_loadings <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = AA)) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                 arrow = arrow(length = unit(0.2, "cm")), color = "gray50") +
    geom_text(size = 4, fontface = "bold") +
    labs(title = paste(chain, "AA Loadings"), x = "PC1", y = "PC2") +
    theme_minimal() +
    coord_fixed()
  
  ggsave(file.path("results", "tcr_comparison", "cdr3_features", 
                   paste0("aa_pca_loadings_", chain, ".pdf")), p_loadings, width = 6, height = 6)
  
  # --- 2. PERMANOVA for each comparison ---
  for (comp_name in c("group", "tissue", "diagnosis", "tissue_diagnosis")) {
    comp <- aa_comparisons[[comp_name]]
    
    tryCatch({
      formula <- as.formula(paste("aa_matrix ~", comp$var))
      perm_result <- adonis2(formula, data = aa_data, 
                             method = "euclidean", permutations = 999)
      aa_analysis_results[[paste0(chain, "_permanova_", comp_name)]] <- perm_result
      message(paste("\n", chain, "PERMANOVA Results for", comp_name, ":"))
      print(perm_result)
    }, error = function(e) {
      message(paste("PERMANOVA failed for", chain, comp_name, ":", e$message))
    })
  }
  
  # Interaction PERMANOVA
  tryCatch({
    perm_interaction <- adonis2(aa_matrix ~ tissue * group, data = aa_data, 
                                method = "euclidean", permutations = 999)
    aa_analysis_results[[paste0(chain, "_permanova_interaction")]] <- perm_interaction
    message(paste("\n", chain, "PERMANOVA Results (tissue * group interaction):"))
    print(perm_interaction)
  }, error = function(e) {
    message(paste("Interaction PERMANOVA failed for", chain, ":", e$message))
  })
  
  # --- 3. Heatmaps for each comparison ---
  for (comp_name in c("group", "tissue", "diagnosis", "tissue_diagnosis")) {
    comp <- aa_comparisons[[comp_name]]
    
    tryCatch({
      aa_summary <- aa_data |>
        dplyr::group_by(.data[[comp$var]]) |>
        dplyr::summarize(
          dplyr::across(all_of(aa_cols), mean, na.rm = TRUE),
          .groups = "drop"
        )
      
      mat <- aa_summary |>
        dplyr::select(all_of(aa_cols)) |>
        as.matrix()
      rownames(mat) <- aa_summary[[comp$var]]
      mat_scaled <- scale(mat)
      
      pdf(file.path("results", "tcr_comparison", "cdr3_features", 
                    paste0("aa_heatmap_", chain, "_", comp_name, ".pdf")), 
          width = 10, height = max(4, nrow(mat) * 0.5))
      print(Heatmap(mat_scaled, name = "Z-score", 
                    column_title = paste(chain, "AA Composition by", comp_name),
                    cluster_rows = TRUE, cluster_columns = TRUE))
      dev.off()
    }, error = function(e) {
      message(paste("Heatmap failed for", chain, comp_name, ":", e$message))
    })
  }
  
  # --- 4. Individual AA boxplots with statistics for each comparison ---
  aa_long <- aa_data |>
    tidyr::pivot_longer(cols = all_of(aa_cols), names_to = "amino_acid", values_to = "proportion")
  
  for (comp_name in c("group", "tissue", "diagnosis")) {
    comp <- aa_comparisons[[comp_name]]
    
    tryCatch({
      # Statistical test for each AA
      aa_stats <- aa_long |>
        dplyr::group_by(amino_acid) |>
        dplyr::summarize(
          p_value = tryCatch(
            kruskal.test(proportion ~ .data[[comp$var]])$p.value,
            error = function(e) NA
          ),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          p_adj = p.adjust(p_value, method = "BH"),
          sig = dplyr::case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01 ~ "**",
            p_adj < 0.05 ~ "*",
            TRUE ~ ""
          )
        )
      
      aa_analysis_results[[paste0(chain, "_aa_stats_", comp_name)]] <- aa_stats
      
      p_box <- aa_long |>
        dplyr::left_join(aa_stats, by = "amino_acid") |>
        ggplot(aes(x = .data[[comp$var]], y = proportion, fill = .data[[comp$var]])) +
        geom_boxplot(outlier.size = 0.5) +
        facet_wrap(~ paste0(amino_acid, sig), scales = "free_y", ncol = 5) +
        labs(title = paste(chain, "AA Composition by", comp_name), y = "Proportion") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      if (!is.null(comp$colors)) {
        p_box <- p_box + scale_fill_manual(values = comp$colors)
      }
      
      ggsave(file.path("results", "tcr_comparison", "cdr3_features", 
                       paste0("aa_boxplots_", chain, "_", comp_name, ".pdf")), 
             p_box, width = 12, height = 10)
    }, error = function(e) {
      message(paste("Boxplots failed for", chain, comp_name, ":", e$message))
    })
  }
  
  # --- 5. UMAP ---
  if (nrow(aa_matrix) > 5000) {
    set.seed(42)
    idx <- sample(nrow(aa_matrix), 5000)
    aa_matrix_sub <- aa_matrix[idx, ]
    aa_data_sub <- aa_data[idx, ]
  } else {
    aa_matrix_sub <- aa_matrix
    aa_data_sub <- aa_data
  }
  
  umap_result <- umap(aa_matrix_sub)
  
  umap_df <- data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    tissue = aa_data_sub$tissue,
    group = aa_data_sub$group,
    diagnosis = aa_data_sub$diagnosis,
    tissue_diagnosis = aa_data_sub$tissue_diagnosis,
    seurat_clusters = aa_data_sub$seurat_clusters
  )
  
  # UMAP plots for each comparison
  for (comp_name in names(aa_comparisons)) {
    comp <- aa_comparisons[[comp_name]]
    
    p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = .data[[comp$var]])) +
      geom_point(alpha = 0.3, size = 0.5) +
      labs(title = paste(chain, "CDR3 AA Composition - UMAP by", comp_name)) +
      theme_minimal()
    
    if (!is.null(comp$colors)) {
      p_umap <- p_umap + scale_color_manual(values = comp$colors)
    }
    
    # Facet by tissue for non-tissue comparisons
    if (comp_name %in% c("group", "diagnosis", "cluster")) {
      p_umap <- p_umap + facet_wrap(~ tissue)
    }
    
    ggsave(file.path("results", "tcr_comparison", "cdr3_features", 
                     paste0("aa_umap_", chain, "_", comp_name, ".pdf")), 
           p_umap, width = 10, height = 5)
  }
}

# Export statistics
stats_to_export <- aa_analysis_results[grep("_aa_stats_", names(aa_analysis_results))]
if (length(stats_to_export) > 0) {
  writexl::write_xlsx(
    stats_to_export,
    file.path("results", "tcr_comparison", "cdr3_features", "aa_statistics.xlsx")
  )
}



# 3.3 Positional entropy ----

for (chain in c("TRA", "TRB")) {
    for (comp_name in c("tissue_diagnosis", "group")) {
        comp <- comparison_vars[[comp_name]]

        tryCatch({
            entropy_plot <- positionalEntropy(
                sc_tcr,
                chain = chain,
                group.by = comp$var
            )

            ggsave(
                file.path("results", "tcr_comparison", "cdr3_features",
                         paste0("positional_entropy_", chain, "_", comp_name, ".pdf")),
                entropy_plot, width = 12, height = 6
            )
        }, error = function(e) {
            message(paste("Entropy failed for", chain, comp_name, ":", e$message))
        })
    }
}
