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
cd8tem_3 <- subset(sc_merge, subset = cluster %in% c("CD8TEM_3"))

# Function to process and plot projections
process_diagnosis <- function(data, diagnosis = NULL, ref, title) {
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
    ref_cd8_terekhova,
    all_projected,
    linesize = .5,
    pointsize = 1,
    ref.size = 0.1,
    cols = colors_terekhova,
  ) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(
        color = "black",
        linewidth = 1,
        fill = NA
      ),
      aspect.ratio = 0.7
    ) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    ggtitle(title)

  # Save plot
  ggsave(
    file.path(
      "results",
      "projectil",
      paste0("cd8tem_3", suffix, "_", ref_name, "_projectil.pdf")
    ),
    width = 7,
    height = 4
  )

  return(projected)
}

# Load reference map
ref_cd8_terekhova <- qread(file.path(
  "objects",
  "conventional_cd8_projectil_reference.qs"
))

str(ref_cd8_terekhova@assays$RNA$counts)
str(ref_cd8_terekhova@assays$RNA)

# try reproducing the colors from Terekhova et al.
colors_terekhova <- c(
  "Naive" = "#854c90",
  "Naive-IFN" = "#a3a5aa",
  "Trm" = "#bfa0d1",
  "Tmem KLRC2+" = "#4fa6a2",
  "Tem GZMK+" = "#406f9d",
  "Tem GZMB+" = "#e89e34",
  "NKT-like" = "#99ad6e",
  "Temra" = "#dd587d",
  "Proliferative" = "#fc8f7e",
  "HLA-DR+" = "#d22b90",
  "Tcm CCR4+" = "#eb8123",
  "Tcm CCR4-" = "#99ad6e"
)

# Process all CD8TEM_3 cells
all_projected <- process_diagnosis(
  cd8tem_3,
  diagnosis = NULL,
  ref = ref_cd8_terekhova,
  title = "CD8TEM_3 cells projected onto Terekhova et al. CD8"
)

# Process each diagnosis
diagnoses <- list(
  "CIAP" = "CIAP",
  "CIDP" = "CIDP",
  "GBS" = "GBS"
)

# Process each diagnosis group
projections <- lapply(
  diagnoses,
  function(x) {
    process_diagnosis(cd8tem_3, diagnosis = x, ref = ref_cd8_terekhova)
  }
)

# Print summary of cells per diagnosis
diagnosis_counts <- table(cd8tem_3$diagnosis)
print("Number of cells per diagnosis:")
print(diagnosis_counts)
