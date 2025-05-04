##################################################
# differential gene expression analysis
# requires running annotate.R first
##################################################

# load libraries
library(tidyverse)
library(Seurat)
library(qs)
library(scMisc)
library(readxl)

# Helper functions ----
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
        pt.size = .1
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
create_violin_plots_for_object <- function(seurat_obj, base_suffix, feature_sets) {
    for (set in feature_sets) {
        suffix <- paste0(set$name, "_", base_suffix)
        create_save_violin_plot(
            seurat_obj = seurat_obj,
            features = set$features,
            filename_suffix = suffix
        )
    }
}

# Function to create expression heatmap with custom filename
create_expression_heatmap <- function(seurat_obj, features, group_by) {
    # Calculate average expression
    avg_expr <- AverageExpression(
        seurat_obj,
        features = features,
        group.by = group_by
    )$RNA
    return(avg_expr)
}

# load data ----
# read preprocessed data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

# Create subsets using direct subsetting
main_diagnoses <- c("CIDP", "GBS", "CIAP", "CTRL")
cd8_nk_clusters <- c(
    "CD8_naive",
    "CD8TEM_1",
    "CD8TEM_2",
    "CD8_NK",
    "NKCD56bright",
    "NKCD56dim"
)
cd8_nk_cluster <- c("CD8_NK")

sc_merge_main_diagnosis <- subset(
    sc_merge,
    subset = diagnosis %in% main_diagnoses
)

sc_merge_main_diagnosis_pbmc <- subset(
    sc_merge,
    subset = diagnosis %in% main_diagnoses & tissue == "PBMC"
)

sc_merge_main_diagnosis_csf <- subset(
    sc_merge,
    subset = diagnosis %in% main_diagnoses & tissue == "CSF"
)

sc_merge_main_diagnosis_pbmc_t_nk <- subset(
    sc_merge,
    subset = diagnosis %in%
        main_diagnoses &
        tissue == "PBMC" &
        cluster %in% cd8_nk_clusters
)

sc_merge_main_diagnosis_pbmc_cd8_nk <- subset(
    sc_merge,
    subset = diagnosis %in%
        main_diagnoses &
        tissue == "PBMC" &
        cluster %in% cd8_nk_cluster
)

# Load marker genes
markers <- read_csv(file.path("lookup", "markers.csv"))

# Define feature sets
feature_sets <- list(
    list(name = "bnb", features = markers$BNB),
    list(name = "bnb_ligands", features = markers$BNB_ligands),
    list(name = "chemokines", features = markers$chemokines)
)

# Create violin plots for each Seurat object
create_violin_plots_for_object(sc_merge_main_diagnosis, "main", feature_sets)
create_violin_plots_for_object(sc_merge_main_diagnosis_csf, "main_csf", feature_sets)
create_violin_plots_for_object(sc_merge_main_diagnosis_pbmc_cd8_nk, "pbmc_cd8_nk", feature_sets)
create_violin_plots_for_object(sc_merge_main_diagnosis_pbmc_t_nk, "pbmc_t_nk", feature_sets)

# Create heatmaps using the helper function with custom names
main_bnb_avg <- create_expression_heatmap(
    sc_merge_main_diagnosis,
    markers$BNB,
    group_by = "diagnosis"
)

main_pbmc_bnb_avg <- create_expression_heatmap(
    sc_merge_main_diagnosis_pbmc,
    markers$BNB,
    group_by = "diagnosis"
)

scMisc::pHeatmap(main_bnb_avg, scale = "row", )

scMisc::pHeatmap(main_pbmc_bnb_avg, scale = "row", )
