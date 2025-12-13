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

# load helper functions ----
source(file.path("scripts", "deg_hypothesis_helper.R"))

# load data ----
# read preprocessed data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

# Create subsets using direct subsetting
main_diagnoses <- c("CIDP", "GBS", "CIAP", "CTRL")
tc_nk_clusters <- c(
    c(
        "CD4naive_1",
        "CD4TCM_1",
        "CD4TCM_2",
        "CD4TEM",
        "Treg",
        "MAIT",
        "gdT",
        "CD4CTL",
        "CD8naive",
        "CD8TCM",
        "CD8TEM_1",
        "CD8TEM_2",
        "CD8TEM_3",
        "NKCD56bright_1",
        "NKCD56bright_2",
        "NKCD56dim"
    )
)

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
        cluster %in% tc_nk_clusters
)

# Load marker genes
markers <- read_csv(file.path("lookup", "markers.csv"))

# Define feature sets
feature_sets <- list(
    list(name = "bnb", features = markers$BNB),
    list(name = "bnb_ligands", features = markers$BNB_ligands),
    list(name = "chemokines", features = markers$chemokines),
    list(name = "cytotoxicity", features = markers$cytotoxicity),
    list(name = "exhaustion", features = markers$exhaustion)
)

# Prepare datasets for analysis ----
datasets <- list(
    main = sc_merge_main_diagnosis,
    main_csf = sc_merge_main_diagnosis_csf,
    main_pbmc = sc_merge_main_diagnosis_pbmc,
    pbmc_t_nk = sc_merge_main_diagnosis_pbmc_t_nk
)

# Create violin plots for all datasets ----
for (dataset_name in names(datasets)) {
    create_violin_plots_for_object(
        datasets[[dataset_name]],
        dataset_name,
        feature_sets
    )
}

# Prepare marker sets for heatmaps ----
marker_sets <- list(
    bnb = markers$BNB,
    cytotoxicity = markers$cytotoxicity,
    exhaustion = markers$exhaustion
)

# Create all heatmaps ----
process_heatmaps(
    datasets = datasets,
    marker_sets = marker_sets,
    group_by = "diagnosis",
    cluster_cols = TRUE
)
