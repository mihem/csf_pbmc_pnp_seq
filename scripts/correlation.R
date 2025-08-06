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

# Function to create regression plots using limma fit object
create_regression_plot_limma <- function(cell_type, fit_object, props_data, meta_data) {
    # Extract statistics from limma fit object
    coef_val <- round(fit_object$coefficients[cell_type, "onls_at_lumbar_puncture"], 3)
    p_val <- fit_object$p.value[cell_type, "onls_at_lumbar_puncture"]
    t_val <- fit_object$t[cell_type, "onls_at_lumbar_puncture"]
    
    # Calculate residuals from the limma model - these are already adjusted for age and sex
    design_matrix <- fit_object$design
    fitted_values <- design_matrix %*% t(fit_object$coefficients)
    residuals_matrix <- props_data - t(fitted_values)
    
    # Get residuals for this specific cell type (already adjusted for age and sex)
    cell_residuals <- residuals_matrix[cell_type, ]
    
    # Calculate simple correlation for comparison (unadjusted)
    simple_pearson <- round(cor(meta_data$onls_at_lumbar_puncture, props_data[cell_type, ], 
                               method = "pearson", use = "complete.obs"), 3)
    
    # For plotting, we need to show the relationship between ONLS and the cell proportions
    # We can use the partial residuals approach: residuals + ONLS effect
    onls_effect <- fit_object$coefficients[cell_type, "onls_at_lumbar_puncture"] * meta_data$onls_at_lumbar_puncture
    partial_residuals <- cell_residuals + onls_effect
    
    # Calculate correlation from the model relationship
    # The correlation can be derived from the t-statistic and degrees of freedom
    df_residual <- fit_object$df.residual[1]  # All cell types have the same df.residual
    adjusted_pearson <- round(t_val / sqrt(t_val^2 + df_residual), 3)
    
    # Format p-value
    p_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
    
    # Create plot data using partial residuals
    plot_data <- data.frame(
        onls = meta_data$onls_at_lumbar_puncture,
        partial_residuals = partial_residuals,
        stringsAsFactors = FALSE
    )
    
    # Remove incomplete cases
    complete_cases <- complete.cases(plot_data)
    plot_data <- plot_data[complete_cases, ]
    
    # Create plot
    p <- ggplot(plot_data, aes(x = onls, y = partial_residuals)) +
        geom_point() +
        geom_smooth(method = "lm", se = TRUE) +
        theme_classic() +
        labs(
            title = gsub("_", " ", cell_type),
            x = "ONLS at lumbar puncture",
            y = "Proportion (adjusted for age & sex)",
            subtitle = paste("Slope:", coef_val, "| r =", simple_pearson, "| r.adj =", adjusted_pearson, "|", p_text, "| t =", round(t_val, 2))
        )
    
    return(p)
}

# Function to extract correlation statistics from limma fit object
extract_limma_correlation_stats <- function(fit_object, props_data, meta_data) {
    cell_types <- rownames(fit_object$coefficients)
    
    results <- data.frame(
        cell_type = cell_types,
        # Limma statistics (these are already adjusted for age and sex)
        limma_coef = fit_object$coefficients[, "onls_at_lumbar_puncture"],
        limma_t_stat = fit_object$t[, "onls_at_lumbar_puncture"],
        limma_p_value = fit_object$p.value[, "onls_at_lumbar_puncture"],
        limma_sigma = fit_object$sigma,
        # Simple correlations (unadjusted)
        simple_pearson_r = NA,
        simple_pearson_p = NA,
        # Adjusted correlation coefficient derived from limma t-statistic
        adjusted_pearson_r = NA,
        stringsAsFactors = FALSE
    )
    
    for (i in 1:length(cell_types)) {
        cell_type <- cell_types[i]
        
        # Simple correlation (unadjusted)
        simple_test <- cor.test(meta_data$onls_at_lumbar_puncture, 
                               props_data[cell_type, ], 
                               method = "pearson", 
                               use = "complete.obs")
        results$simple_pearson_r[i] <- simple_test$estimate
        results$simple_pearson_p[i] <- simple_test$p.value
        
        # Calculate adjusted correlation coefficient from t-statistic and degrees of freedom
        # This represents the partial correlation controlling for age and sex
        t_val <- fit_object$t[cell_type, "onls_at_lumbar_puncture"]
        df_residual <- fit_object$df.residual[1]  # All cell types have the same df.residual
        results$adjusted_pearson_r[i] <- t_val / sqrt(t_val^2 + df_residual)
    }
    
    # Add adjusted p-values
    results$limma_p_adj <- p.adjust(results$limma_p_value, method = "BH")
    results$simple_pearson_p_adj <- p.adjust(results$simple_pearson_p, method = "BH")
    
    return(results)
}

# Get topTable results
top_results <- limma::topTable(
    fit,
    coef = "onls_at_lumbar_puncture",
    number = Inf,
    sort.by = "p"
)

# Extract correlation statistics using limma fit object
correlation_stats <- extract_limma_correlation_stats(fit, props$TransformedProps, meta_lookup)

# Save correlation results
write.csv(correlation_stats, 
          file.path("results", "correlation", "onls_correlation_statistics_limma.csv"), 
          row.names = FALSE)

# Create plots for significant results based on limma p-values
significant_cells <- correlation_stats$cell_type[correlation_stats$limma_p_value < 0.05 & !is.na(correlation_stats$limma_p_value)]

# Generate plots for significant cell types using limma fit object
regression_plots <- lapply(significant_cells, function(cell_type) {
    create_regression_plot_limma(cell_type, fit, props$TransformedProps, meta_lookup)
})

# Save plots
for (i in seq_along(regression_plots)) {
    cell_name <- significant_cells[i]
    ggsave(
        plot = regression_plots[[i]],
        filename = file.path("results", "correlation", paste0("onls_", cell_name, "_regression_limma.pdf")),
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
        filename = file.path("results", "correlation", "onls_all_significant_regressions_limma.pdf"),
        width = 16,
        height = 6 * ceiling(length(regression_plots) / 2)
    )
}
