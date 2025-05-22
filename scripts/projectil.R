###############################################################################
# use ProjecTILs to project TILs from the reference map to the query
# requires running annotate.R first
###############################################################################

# load libraries ----
library(qs)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ProjecTILs)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)
cd8_nk <- subset(sc_merge, subset = cluster %in% c("CD8_NK"))

# Function to process and plot projections
process_diagnosis <- function(data, diagnosis = NULL, ref) {
    ref_name <- deparse(substitute(ref))
    # If diagnosis is provided, subset the data
    if (!is.null(diagnosis)) {
        data <- subset(data, subset = diagnosis %in% diagnosis)
        suffix <- paste0("_", tolower(paste(diagnosis, collapse = "_")))
    } else {
        suffix <- ""
    }

    # Create projection
    projected <- make.projection(
        data,
        ref = ref,
        ncores = 8,
        filter.cells = FALSE
    )

    # Create and save plot
    plot.projection(
        ref,
        projected,
        linesize = 0.5,
        pointsize = 0.5
    )

    # Save plot
    ggsave(
        file.path(
            "results",
            "projectil",
            paste0("cd8_nk", suffix, "_", ref_name, "_projectil.pdf")
        ),
        width = 8,
        height = 6
    )

    return(projected)
}

# Load reference map
ref <- load.reference.map()

ref_cd8_terekhova <- qread(file.path(
    "objects",
    "conventional_cd8_projectil_reference.qs"
))

str(ref@assays$RNA$counts)
str(ref_cd8_terekhova@assays$RNA$counts)
str(ref_cd8_terekhova@assays$RNA)

# Process all CD8_NK cells
all_projected <- process_diagnosis(cd8_nk, diagnosis = NULL, ref = ref)
all_projected <- process_diagnosis(cd8_nk, diagnosis = NULL, ref = ref_cd8_terekhova)

# Process each diagnosis
diagnoses <- list(
    "CIAP" = "CIAP",
    "CIDP" = "CIDP",
    "GBS" = "GBS"
)

# Process each diagnosis group
projections <- lapply(
    diagnoses,
    function(x) process_diagnosis(cd8_nk, diagnosis = x, ref = ref)
)

# Print summary of cells per diagnosis
diagnosis_counts <- table(cd8_nk$diagnosis)
print("Number of cells per diagnosis:")
print(diagnosis_counts)
