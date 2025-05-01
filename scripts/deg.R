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

    filename <- file.path(
        "results",
        "de",
        paste0("vln_bnb_", filename_suffix, ".png")
    )

    ggplot2::ggsave(
        filename = filename,
        plot = plot,
        width = 10,
        height = 20
    )

    return(plot)
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
bnb_markers <- read_csv(file.path("lookup", "markers.csv"))

# Create and save violin plots using the helper function
vln_bnb_plot <- create_save_violin_plot(
    sc_merge_main_diagnosis,
    bnb_markers$BNB,
    filename_suffix = "main"
)

vln_bnb_pbmc_cd8_nk_plot <- create_save_violin_plot(
    sc_merge_main_diagnosis_pbmc_cd8_nk,
    bnb_markers$BNB,
    filename_suffix = "pbmc_cd8_nk"
)

vln_bnb_pbmc_t_nk_plot <- create_save_violin_plot(
    sc_merge_main_diagnosis_pbmc_t_nk,
    bnb_markers$BNB,
    filename_suffix = "pbmc_t_nk"
)

# Create heatmaps using the helper function with custom names
main_bnb_avg <- create_expression_heatmap(
    sc_merge_main_diagnosis,
    bnb_markers$BNB,
    group_by = "diagnosis"
)

main_pbmc_bnb_avg <- create_expression_heatmap(
    sc_merge_main_diagnosis_pbmc,
    bnb_markers$BNB,
    group_by = "diagnosis"
)

scMisc::pHeatmap(main_bnb_avg, scale = "row", )

scMisc::pHeatmap(main_pbmc_bnb_avg, scale = "row", )
