##################################################
# Helper functions for deg_hypothesis.R
##################################################

# Function to create and save violin plots
create_save_violin_plot <- function(
    seurat_obj,
    features,
    group_by = "diagnosis",
    filename_suffix = ""
) {
    filename <- file.path(
        "results",
        "de",
        paste0("vln_", filename_suffix, ".png")
    )

    # Skip if file already exists
    if (file.exists(filename)) {
        message(sprintf("Skipping existing file: %s", filename))
        return(invisible(NULL))
    }

    # remove NA values from features
    features <- features[!is.na(features)]

    plot <- VlnPlot(
        seurat_obj,
        features = features,
        group.by = group_by,
        pt.size = 0
    ) +
        NoLegend() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
        )

    height <- ceiling(length(features) / 4 * 3)

    ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = 10,
        height = height,
        limitsize = FALSE
    )
}

# Simple function to create violin plots for a Seurat object with multiple feature sets
create_violin_plots_for_object <- function(
    seurat_obj,
    base_suffix,
    feature_sets
) {
    for (set in feature_sets) {
        suffix <- paste0(set$name, "_", base_suffix)
        create_save_violin_plot(
            seurat_obj = seurat_obj,
            features = set$features,
            filename_suffix = suffix
        )
    }
}

# Function to create expression heatmap
create_expression_heatmap <- function(seurat_obj, features, group_by) {
    # Calculate average expression
    avg_expr <- AverageExpression(
        seurat_obj,
        features = features,
        group.by = group_by
    )$RNA
    return(avg_expr)
}

# Function to create and save heatmap
create_save_heatmap <- function(
    avg_expr,
    filename,
    scale = "none",
    cluster_cols = TRUE
) {
    # Create directory for the target filename
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)

    # Call scMisc::pHeatmap
    # This will create a file named 'results/heatmap/hm_avg_expr.pdf' because
    # the variable passed is named 'avg_expr'
    scMisc::pHeatmap(
        avg_expr,
        scale = scale,
        cluster_cols = cluster_cols
    )

    # Move the generated file to the desired filename
    generated_file <- file.path("results", "heatmap", "hm_avg_expr.pdf")
    if (file.exists(generated_file)) {
        file.rename(generated_file, filename)
    } else {
        warning("Expected heatmap file not found: ", generated_file)
    }
}

# Function to process heatmaps for multiple datasets and markers
process_heatmaps <- function(
    datasets,
    marker_sets,
    group_by = "diagnosis",
    cluster_cols = TRUE,
    output_dir = file.path("results", "heatmap")
) {
    # Create output directory
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    for (dataset_name in names(datasets)) {
        for (marker_name in names(marker_sets)) {
            # Create average expression
            avg_expr <- create_expression_heatmap(
                datasets[[dataset_name]],
                marker_sets[[marker_name]],
                group_by = group_by
            )

            filename <- file.path(
                output_dir,
                paste0("hm_", dataset_name, "_", marker_name, "_avg.pdf")
            )

            create_save_heatmap(
                avg_expr = avg_expr,
                filename = filename,
                cluster_cols = FALSE,
                scale = "row"
            )
        }
    }
}
