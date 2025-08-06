########################################################
# correlation analysis helper function
########################################################
library(dplyr)
library(ggplot2)
library(stringr)
library(speckle)
library(limma)
library(patchwork)

# Function to perform correlation analysis for a given severity score
analyze_severity_correlation <- function(
    seurat_obj,
    severity_score,
    lookup_data,
    tissue_type
) {
    # Filter seurat object by tissue type
    seu_tissue <- subset(seurat_obj, tissue == tissue_type)

    # Filter out samples with missing severity scores
    seu_filtered <- seu_tissue[,
        !is.na(seu_tissue@meta.data[[severity_score]])
    ]

    # Check if we have enough samples
    if (ncol(seu_filtered) < 5) {
        warning(paste("Not enough samples with", severity_score, "data"))
        return(NULL)
    }

    # Get transformed proportions
    props <- speckle::getTransformedProps(
        clusters = seu_filtered@meta.data[["cluster"]],
        sample = seu_filtered@meta.data[["patient"]],
        transform = "logit"
    )

    # Create metadata lookup
    meta_lookup <- tibble::tibble(
        patient = colnames(props$TransformedProps)
    ) |>
        dplyr::left_join(lookup_data, by = "patient") |>
        dplyr::distinct(patient, .keep_all = TRUE)

    # Create design matrix
    formula <- paste0("~0 + ", severity_score, " + sex + age")
    my_design <- model.matrix(as.formula(formula), data = meta_lookup)

    # Fit limma model
    fit <- limma::lmFit(props$TransformedProps, my_design)
    fit <- limma::eBayes(fit)

    return(list(
        fit = fit,
        props = props,
        meta_lookup = meta_lookup,
        severity_score = severity_score
    ))
}

# Function to create combined correlation plot for all cell types
create_combined_correlation_plot <- function(
    fit_object,
    props_data,
    meta_data,
    severity_score,
    correlation_threshold = 0.3
) {
    # Get correlation statistics
    correlation_stats <- extract_limma_correlation_stats(
        fit_object,
        props_data,
        meta_data,
        severity_score
    )

    # Filter for high correlation cell types
    high_corr_cells <- correlation_stats$cell_type[
        abs(correlation_stats$adjusted_pearson_r) > correlation_threshold &
            !is.na(correlation_stats$adjusted_pearson_r)
    ]

    if (length(high_corr_cells) == 0) {
        return(NULL)
    }

    # Order cell types by absolute correlation (strongest first)
    cell_order <- correlation_stats[
        correlation_stats$cell_type %in% high_corr_cells,
    ] |>
        dplyr::arrange(desc(abs(adjusted_pearson_r))) |>
        dplyr::pull(cell_type)

    # Create nice labels
    score_label <- gsub("_", " ", severity_score)
    score_label <- stringr::str_to_title(score_label)

    # Create individual plots for each cell type
    individual_plots <- list()

    for (cell_type in cell_order) {
        # Get statistics for this cell type
        cell_stats <- correlation_stats[
            correlation_stats$cell_type == cell_type,
        ]

        # Calculate residuals from the limma model
        design_matrix <- fit_object$design
        fitted_values <- design_matrix %*% t(fit_object$coefficients)
        residuals_matrix <- props_data - t(fitted_values)

        # Get residuals for this specific cell type
        cell_residuals <- residuals_matrix[cell_type, ]

        # Calculate partial residuals
        score_effect <- fit_object$coefficients[cell_type, severity_score] *
            meta_data[[severity_score]]
        partial_residuals <- cell_residuals + score_effect

        # Create data frame for this cell type
        plot_data <- data.frame(
            severity_score = meta_data[[severity_score]],
            partial_residuals = partial_residuals,
            stringsAsFactors = FALSE
        )

        # Remove incomplete cases
        plot_data <- plot_data[complete.cases(plot_data), ]

        # Format p-value
        p_text <- ifelse(
            cell_stats$limma_p_value < 0.001,
            "p < 0.001",
            paste("p =", round(cell_stats$limma_p_value, 3))
        )

        # Create individual plot
        p <- ggplot(plot_data, aes(x = severity_score, y = partial_residuals)) +
            geom_point(alpha = 0.7, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, color = "blue", size = 1) +
            theme_classic() +
            theme(
                plot.title = element_text(size = 11, face = "bold"),
                plot.subtitle = element_text(size = 9),
                axis.text = element_text(size = 9),
                axis.title = element_text(size = 10)
            ) +
            labs(
                title = gsub("_", " ", cell_type),
                subtitle = paste(
                    "r =",
                    round(cell_stats$adjusted_pearson_r, 3)
                ),
                x = score_label,
                y = "Proportion\n(adj. for age & sex)"
            )

        individual_plots[[cell_type]] <- p
    }

    # Combine plots using patchwork
    combined_plot <- patchwork::wrap_plots(individual_plots, ncol = 3)

    return(combined_plot)
}

# Function to extract correlation statistics from limma fit object
extract_limma_correlation_stats <- function(
    fit_object,
    props_data,
    meta_data,
    severity_score
) {
    cell_types <- rownames(fit_object$coefficients)

    results <- data.frame(
        cell_type = cell_types,
        # Limma statistics (these are already adjusted for age and sex)
        limma_coef = fit_object$coefficients[, severity_score],
        limma_t_stat = fit_object$t[, severity_score],
        limma_p_value = fit_object$p.value[, severity_score],
        limma_sigma = fit_object$sigma,
        # Adjusted correlation coefficient derived from limma t-statistic
        adjusted_pearson_r = NA,
        stringsAsFactors = FALSE
    )

    for (i in 1:length(cell_types)) {
        cell_type <- cell_types[i]

        # Calculate adjusted correlation coefficient from t-statistic and degrees of freedom
        # This represents the partial correlation controlling for age and sex
        t_val <- fit_object$t[cell_type, severity_score]
        df_residual <- fit_object$df.residual[1] # All cell types have the same df.residual
        results$adjusted_pearson_r[i] <- t_val / sqrt(t_val^2 + df_residual)
    }

    # Add adjusted p-values
    results$limma_p_adj <- p.adjust(results$limma_p_value, method = "BH")

    return(results)
}

# Function to process correlation results and generate plots
process_correlation_results <- function(
    analysis_result,
    output_dir,
    correlation_threshold,
    tissue_type
) {
    if (is.null(analysis_result)) {
        return(NULL)
    }

    fit <- analysis_result$fit
    props <- analysis_result$props
    meta_lookup <- analysis_result$meta_lookup
    severity_score <- analysis_result$severity_score

    # Get topTable results
    top_results <- limma::topTable(
        fit,
        coef = severity_score,
        number = Inf,
        sort.by = "p"
    )

    # Extract correlation statistics using limma fit object
    correlation_stats <- extract_limma_correlation_stats(
        fit,
        props$TransformedProps,
        meta_lookup,
        severity_score
    )

    # Save correlation results
    output_file <- file.path(
        output_dir,
        paste0(
            "correlation_",
            severity_score,
            "_",
            tissue_type,
            ".csv"
        )
    )
    write.csv(
        correlation_stats,
        output_file,
        row.names = FALSE
    )

    # Create plots for results with high correlation coefficients
    high_corr_cells <- correlation_stats$cell_type[
        abs(correlation_stats$adjusted_pearson_r) > correlation_threshold &
            !is.na(correlation_stats$adjusted_pearson_r)
    ]

    if (length(high_corr_cells) > 0) {
        # Generate combined plot for all high correlation cell types
        combined_plot <- create_combined_correlation_plot(
            fit,
            props$TransformedProps,
            meta_lookup,
            severity_score,
            correlation_threshold
        )

        if (!is.null(combined_plot)) {
            # Calculate plot dimensions based on number of cell types
            # 3 plots per row, so calculate number of rows needed
            n_plots <- length(high_corr_cells)
            n_rows <- ceiling(n_plots / 3)

            # Base height per row (adjust as needed)
            height_per_row <- 4
            plot_height <- n_rows * height_per_row

            # Keep width consistent but ensure minimum dimensions
            plot_width <- 12

            # Save the combined plot
            plot_file <- file.path(
                output_dir,
                paste0(
                    "correlation_",
                    severity_score,
                    "_",
                    tissue_type,
                    ".pdf"
                )
            )
            ggsave(
                plot = combined_plot,
                filename = plot_file,
                width = plot_width,
                height = plot_height
            )

            cat(
                "Generated combined plot for",
                tissue_type,
                severity_score,
                "with",
                length(high_corr_cells),
                "cell types (|r| >",
                correlation_threshold,
                ") - dimensions:",
                plot_width,
                "x",
                plot_height,
                "\n"
            )
        }
    } else {
        cat(
            "No correlations with |r| >",
            correlation_threshold,
            "found for",
            tissue_type,
            severity_score,
            "\n"
        )
    }

    return(list(
        correlation_stats = correlation_stats,
        top_results = top_results,
        high_correlation_cells = high_corr_cells
    ))
}

# Wrapper function to run complete correlation analysis
run_severity_correlation_analysis <- function(
    seurat_obj,
    severity_scores,
    lookup_data,
    tissue_type,
    correlation_threshold,
    output_dir
) {
    results <- list()

    for (score in severity_scores) {
        cat("Analyzing correlations for:", score, "\n")

        # Check if plot file already exists - skip entire analysis if it does
        plot_file <- file.path(
            "results",
            "correlation",
            paste0("correlation_", score, "_", tissue_type, ".pdf")
        )

        if (file.exists(plot_file)) {
            cat(
                "Plot file already exists for",
                tissue_type,
                score,
                "- skipping entire analysis\n"
            )
            results[[score]] <- list(
                correlation_stats = NULL,
                top_results = NULL,
                high_correlation_cells = NULL,
                skipped = TRUE
            )
            next
        }

        # Run analysis
        analysis_result <- analyze_severity_correlation(
            seurat_obj,
            score,
            lookup_data,
            tissue_type
        )

        # Process results and generate plots
        processed_result <- process_correlation_results(
            analysis_result,
            correlation_threshold = correlation_threshold,
            output_dir = output_dir,
            tissue_type = tissue_type
        )

        results[[score]] <- processed_result
    }

    return(results)
}
