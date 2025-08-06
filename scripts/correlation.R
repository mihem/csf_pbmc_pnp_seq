##################################################
# analyze correlations
##################################################
library(tidyverse)
library(qs)
library(Seurat)
library(rcna)
library(glue)
library(scMisc)
library(viridis)

# load data ----
lookup <- qs::qread(file.path("objects", "lookup.qs"))
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))


# analyze subgroup with high NKCD56bright_1 cells ----
nkbright_high_group <- c("P16", "P17", "P18", "P21", "P23")

lookup_subgroup <- lookup |>
    dplyr::mutate(
        subgroup = dplyr::case_when(
            patient %in% nkbright_high_group ~ "NKhi",
            TRUE ~ as.character(diagnosis)
        )
    )

plot_severity_subgroup <-
    lookup_subgroup |>
    ggplot(
        aes(x = subgroup, y = incat_at_lumbar_puncture, fill = subgroup)
    ) +
    geom_boxplot() +
    geom_jitter(
        width = 0.2,
        height = 0,
        alpha = 0.5
    ) +
    theme_classic() +
    theme(legend.position = "none")

ggsave(
    plot = plot_severity_subgroup,
    filename = file.path("results", "severity", "severity_subgroup.pdf"),
    width = 6,
    height = 4
)

# subset those that have the severity variable
sc_merge$incat_at_lumbar_puncture
sc_merge[["incat_at_lumbar_puncture"]]

rcna_plot <- function(seu_obj, param) {
    # check if file exists
    if (
        file.exists(file.path(
            "results",
            "correlation",
            paste0("rcna_", param, ".pdf")
        ))
    ) {
        return(paste0("File ", param, " already exists, skipping."))
    }
    seu_obj_filter <- seu_obj[, !is.na(seu_obj[[param]])]
    seu_obj_filter$sex_numeric <- as.numeric(seu_obj_filter$sex == "male")

    seu_obj_filter <- association.Seurat(
        seurat_object = seu_obj_filter,
        test_var = param,
        samplem_key = "sample",
        graph_use = "RNA_nn",
        verbose = TRUE,
        batches = NULL,
        # no batch variables to include, only works with matched design https://github.com/immunogenomics/cna/issues/11
        covs = c("sex_numeric", "age")
    )

    title <- gsub(pattern = "_", replacement = " ", x = param)
    rcna_fplot <-
        FeaturePlot(
            seu_obj_filter,
            features = c("cna_ncorrs"),
            reduction = "umap.stacas.ss.all",
            pt.size = 0.1,
            order = FALSE,
            coord.fixed = TRUE,
            raster = FALSE,
            alpha = 0.2
        ) +
        viridis::scale_color_viridis(option = "inferno") +
        theme_rect() +
        labs(
            title = title,
            color = "Correlation",
            x = "UMAP1",
            y = "UMAP2"
        )

    ggsave(
        plot = rcna_fplot,
        filename = file.path(
            "results",
            "correlation",
            paste0("rcna_", param, ".pdf")
        ),
        width = 8,
        height = 6
    )
}

rcna_params <- c(
    "incat_at_lumbar_puncture",
    "incat_follow_up",
    "incat_progress",
    "onls_at_lumbar_puncture",
    "onls_follow_up",
    "onls_progress",
    "mrc_sum_score_60_at_lumbar_puncture",
    "mrc_sum_score_60_follow_up",
    "mrc_sum_score_60_progress",
    "csf_protein",
    "dml_ulnar_motoric",
    "cmap_ulnar_motoric",
    "ncv_ulnar_motoric",
    "min_f_latency_ulnar_motoric",
    "cmap_tibial_motoric",
    "ncv_tibial_motoric"
)

lapply(
    rcna_params,
    function(x) {
        rcna_plot(
            seu_obj = sc_merge,
            param = x
        )
    }
)

# standard correlation plots ----
sc_merge_csf <- subset(sc_merge, tissue == "CSF")
seu_filtered <- sc_merge_csf[,
    !is.na(sc_merge_csf@meta.data[["onls_at_lumbar_puncture"]])
]

table(seu_filtered$patient)

props <- speckle::getTransformedProps(
    clusters = seu_filtered@meta.data[["cluster"]],
    sample = seu_filtered@meta.data[["patient"]],
    transform = "logit"
)

colnames(cd8tem_3_props$TransformedProps)

meta_lookup <-
    tibble::tibble(
        patient = colnames(props$TransformedProps)
    ) |>
    dplyr::left_join(lookup, by = "patient") |>
    dplyr::distinct(patient, .keep_all = TRUE)

formula <- "~0 + onls_at_lumbar_puncture + sex + age"
my_design <- model.matrix(as.formula(formula), data = meta_lookup)
fit <- limma::lmFit(props$TransformedProps, my_design)
fit <- limma::eBayes(fit)

 limma::topTable(
    fit,
        coef = "onls_at_lumbar_puncture",
        number = Inf,
        sort.by = "p"
    ) 
    

cd8tem3_plot <-
    cd8tem_3_props_meta |>
    ggplot(aes(
        x = onls_at_lumbar_puncture,
        y = CD8naive
    )) +
    geom_point() +
    theme_classic() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
        title = "CD8TEM_3",
        x = "INCAT at lumbar puncture",
        y = "Proportion"
    )

ggsave(
    plot = cd8tem3_plot,
    filename = "cd8tem3.pdf",
    width = 6,
    height = 4
)


# Get topTable results
top_results <- limma::topTable(
    fit,
    coef = "onls_at_lumbar_puncture",
    number = Inf,
    sort.by = "p"
)

# Function to create regression plots with statistical results
create_regression_plot <- function(cell_type, props_data, meta_data, stats_results) {
    # Get stats for this cell type
    cell_stats <- stats_results[cell_type, ]
    coef_val <- round(cell_stats$logFC, 3)  # This is the regression coefficient (slope), not logFC
    p_val <- cell_stats$P.Value
    adj_p_val <- cell_stats$adj.P.Val
    
    # Calculate correlation coefficients for comparison
    plot_data <- data.frame(
        onls = meta_data$onls_at_lumbar_puncture,
        proportion = props_data[cell_type, ],
        age = meta_data$age,
        sex_numeric = as.numeric(meta_data$sex == "male"),
        stringsAsFactors = FALSE
    )
    
    # Remove rows with missing data
    complete_data <- plot_data[complete.cases(plot_data), ]
    
    # Calculate simple (unadjusted) correlations
    simple_pearson <- round(cor(complete_data$onls, complete_data$proportion, method = "pearson"), 3)
    
    # Method 1: Partial correlations (residuals approach)
    onls_residuals <- residuals(lm(onls ~ age + sex_numeric, data = complete_data))
    prop_residuals <- residuals(lm(proportion ~ age + sex_numeric, data = complete_data))
    partial_pearson <- round(cor(onls_residuals, prop_residuals, method = "pearson"), 3)
    
    # Method 2: Adjusted proportions approach (more intuitive)
    # Adjust proportions by removing age/sex effects
    adjusted_proportions <- residuals(lm(proportion ~ age + sex_numeric, data = complete_data)) + mean(complete_data$proportion)
    adjusted_onls <- residuals(lm(onls ~ age + sex_numeric, data = complete_data)) + mean(complete_data$onls)
    adjusted_pearson <- round(cor(adjusted_onls, adjusted_proportions, method = "pearson"), 3)
    
    # Format p-value
    p_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
    adj_p_text <- ifelse(adj_p_val < 0.001, "adj.p < 0.001", paste("adj.p =", round(adj_p_val, 3)))
    
    # Create plot with adjusted values
    adjusted_data <- data.frame(
        onls_adj = adjusted_onls,
        prop_adj = adjusted_proportions
    )
    
    # Create plot with adjusted data
    p <- ggplot(adjusted_data, aes(x = onls_adj, y = prop_adj)) +
        geom_point() +
        geom_smooth(method = "lm", se = TRUE) +
        theme_classic() +
        labs(
            title = gsub("_", " ", cell_type),
            x = "ONLS (adjusted for age & sex)",
            y = "Proportion (adjusted for age & sex)",
            subtitle = paste("Slope:", coef_val, "| r =", simple_pearson, "| r.adj =", adjusted_pearson, "|", p_text)
        )
    
    return(p)
}

# Function to calculate correlation statistics for all cell types
calculate_correlation_stats <- function(props_data, meta_data) {
    results <- data.frame(
        cell_type = rownames(props_data),
        # Simple correlations (unadjusted)
        simple_pearson_r = NA,
        simple_pearson_p = NA,
        simple_spearman_rho = NA,
        simple_spearman_p = NA,
        # Partial correlations (residuals approach)
        partial_pearson_r = NA,
        partial_pearson_p = NA,
        partial_spearman_rho = NA,
        partial_spearman_p = NA,
        # Adjusted values correlations (adjusted proportions approach)
        adjusted_pearson_r = NA,
        adjusted_pearson_p = NA,
        adjusted_spearman_rho = NA,
        adjusted_spearman_p = NA,
        stringsAsFactors = FALSE
    )
    
    for (i in 1:nrow(props_data)) {
        cell_type <- rownames(props_data)[i]
        
        # Create complete dataset
        complete_data <- data.frame(
            onls = meta_data$onls_at_lumbar_puncture,
            proportion = props_data[cell_type, ],
            age = meta_data$age,
            sex_numeric = as.numeric(meta_data$sex == "male")
        )
        complete_data <- complete_data[complete.cases(complete_data), ]
        
        if (nrow(complete_data) > 4) {  # Need at least 5 points for partial correlation
            # Simple correlations
            simple_pearson_test <- cor.test(complete_data$onls, complete_data$proportion, method = "pearson")
            results$simple_pearson_r[i] <- simple_pearson_test$estimate
            results$simple_pearson_p[i] <- simple_pearson_test$p.value
            
            simple_spearman_test <- cor.test(complete_data$onls, complete_data$proportion, method = "spearman")
            results$simple_spearman_rho[i] <- simple_spearman_test$estimate
            results$simple_spearman_p[i] <- simple_spearman_test$p.value
            
            # Partial correlations (residuals approach)
            onls_residuals <- residuals(lm(onls ~ age + sex_numeric, data = complete_data))
            prop_residuals <- residuals(lm(proportion ~ age + sex_numeric, data = complete_data))
            
            partial_pearson_test <- cor.test(onls_residuals, prop_residuals, method = "pearson")
            results$partial_pearson_r[i] <- partial_pearson_test$estimate
            results$partial_pearson_p[i] <- partial_pearson_test$p.value
            
            partial_spearman_test <- cor.test(onls_residuals, prop_residuals, method = "spearman")
            results$partial_spearman_rho[i] <- partial_spearman_test$estimate
            results$partial_spearman_p[i] <- partial_spearman_test$p.value
            
            # Adjusted values approach (add back means for interpretability)
            adjusted_onls <- onls_residuals + mean(complete_data$onls)
            adjusted_prop <- prop_residuals + mean(complete_data$proportion)
            
            adjusted_pearson_test <- cor.test(adjusted_onls, adjusted_prop, method = "pearson")
            results$adjusted_pearson_r[i] <- adjusted_pearson_test$estimate
            results$adjusted_pearson_p[i] <- adjusted_pearson_test$p.value
            
            adjusted_spearman_test <- cor.test(adjusted_onls, adjusted_prop, method = "spearman")
            results$adjusted_spearman_rho[i] <- adjusted_spearman_test$estimate
            results$adjusted_spearman_p[i] <- adjusted_spearman_test$p.value
        }
    }
    
    # Add adjusted p-values
    results$simple_pearson_p_adj <- p.adjust(results$simple_pearson_p, method = "BH")
    results$simple_spearman_p_adj <- p.adjust(results$simple_spearman_p, method = "BH")
    results$partial_pearson_p_adj <- p.adjust(results$partial_pearson_p, method = "BH")
    results$partial_spearman_p_adj <- p.adjust(results$partial_spearman_p, method = "BH")
    results$adjusted_pearson_p_adj <- p.adjust(results$adjusted_pearson_p, method = "BH")
    results$adjusted_spearman_p_adj <- p.adjust(results$adjusted_spearman_p, method = "BH")
    
    return(results)
}

# Get topTable results
top_results <- limma::topTable(
    fit,
    coef = "onls_at_lumbar_puncture",
    number = Inf,
    sort.by = "p"
)

# Calculate correlation statistics
correlation_stats <- calculate_correlation_stats(props$TransformedProps, meta_lookup)

# Save correlation results
write.csv(correlation_stats, 
          file.path("results", "correlation", "onls_correlation_statistics.csv"), 
          row.names = FALSE)

# Print significant correlations
cat("Significant simple Pearson correlations (p < 0.05):\n")
sig_simple_pearson <- correlation_stats[correlation_stats$simple_pearson_p < 0.05 & !is.na(correlation_stats$simple_pearson_p), ]
if(nrow(sig_simple_pearson) > 0) {
    print(sig_simple_pearson[order(sig_simple_pearson$simple_pearson_p), c("cell_type", "simple_pearson_r", "simple_pearson_p", "simple_pearson_p_adj")])
}

cat("\nSignificant adjusted Pearson correlations (adjusted for age & sex, p < 0.05):\n")
sig_adjusted_pearson <- correlation_stats[correlation_stats$adjusted_pearson_p < 0.05 & !is.na(correlation_stats$adjusted_pearson_p), ]
if(nrow(sig_adjusted_pearson) > 0) {
    print(sig_adjusted_pearson[order(sig_adjusted_pearson$adjusted_pearson_p), c("cell_type", "adjusted_pearson_r", "adjusted_pearson_p", "adjusted_pearson_p_adj")])
}

cat("\nSignificant partial Pearson correlations (residuals method, p < 0.05):\n")
sig_partial_pearson <- correlation_stats[correlation_stats$partial_pearson_p < 0.05 & !is.na(correlation_stats$partial_pearson_p), ]
if(nrow(sig_partial_pearson) > 0) {
    print(sig_partial_pearson[order(sig_partial_pearson$partial_pearson_p), c("cell_type", "partial_pearson_r", "partial_pearson_p", "partial_pearson_p_adj")])
}

# Create plots for top significant results (e.g., adj.P.Val < 0.05)
# significant_cells <- rownames(top_results)[top_results$adj.P.Val < 0.9]
significant_cells <- rownames(top_results)[top_results$P.Val < 0.05]

# Generate plots for significant cell types
regression_plots <- lapply(significant_cells, function(cell_type) {
    create_regression_plot(cell_type, props$TransformedProps, meta_lookup, top_results)
})

# Save plots
for (i in seq_along(regression_plots)) {
    cell_name <- significant_cells[i]
    ggsave(
        plot = regression_plots[[i]],
        filename = file.path("results", "correlation", paste0("onls_", cell_name, "_regression.pdf")),
        width = 8,
        height = 6
    )
}

# Also create a summary plot showing all significant regressions
if (length(significant_cells) > 0) {
    library(cowplot)
    combined_plot <- plot_grid(plotlist = regression_plots, ncol = 2)
    ggsave(
        plot = combined_plot,
        filename = file.path("results", "correlation", "onls_all_significant_regressions.pdf"),
        width = 16,
        height = 6 * ceiling(length(regression_plots) / 2)
    )
}
