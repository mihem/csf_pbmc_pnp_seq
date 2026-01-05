##################################################
# differential gene expression analysis using miloDE
# requires running annotate.R first
##################################################

# Libraries and Setup ----
library(miloR)
library(miloDE)
library(SingleCellExperiment)
library(Seurat)
library(qs)
library(tidyverse)
library(BiocParallel)
library(pals)

# cidp vs ctrl csf
# gbs vs ctrl csf

# Load Preprocessed Data ----
# Load the preprocessed single-cell data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 6)

# Define colors
my_cols_25 <- pals::cols25()

sc_merge_main_csf <- subset(
  sc_merge,
  subset = diagnosis %in% c("CIDP", "GBS", "CIAP", "CTRL") & tissue == "CSF"
)

# Drop emptry levels
sc_merge_main_csf$diagnosis <- droplevels(sc_merge_main_csf$diagnosis)

# Remove unnecessary assays and dimensions
sc_diet <- DietSeurat(
  sc_merge_main_csf,
  counts = TRUE,
  data = TRUE,
  scale.data = TRUE,
  assays = "RNA",
  dimreducs = c("stacas.ss.all", "umap.stacas.ss.all")
)

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(sc_diet, assay = "RNA")

rm(sc_merge)
rm(sc_merge_main_csf)
rm(sc_diet)
scMisc::lss()
gc()

# Neighborhood Assignment ----
# Assign neighborhoods using Milo
set.seed(123)
milo_DE <- assign_neighbourhoods(
  sce,
  k = 20,
  prop = 0.1,
  d = 30,
  order = 2,
  filtering = TRUE,
  reducedDim_name = "STACAS.SS.ALL",
  verbose = TRUE
)

# Save Milo object
qs::qsave(milo_DE, file.path("objects", "milo_DE_csf.qs"))

# Neighborhood Annotation and Plotting ----
# Annotate neighborhoods and create plots
nhoods_sce <- miloR::nhoods(milo_DE)

nhood_stat_ct <- data.frame(
  Nhood = 1:ncol(nhoods_sce),
  Nhood_center = colnames(nhoods_sce)
)

nhood_stat_ct <- miloR::annotateNhoods(
  milo_DE,
  nhood_stat_ct,
  coldata_col = "cluster"
)

milo_hood_plot <- miloDE::plot_milo_by_single_metric(
  milo_DE,
  nhood_stat_ct,
  colour_by = "cluster",
  layout = "UMAP.STACAS.SS.ALL",
  size_range = c(1.5, 3),
  edge_width = c(0.01, 0.05)
) +
  scale_fill_manual(values = my_cols_25, name = "cluster")

ggsave(
  plot = milo_hood_plot,
  file.path("results", "miloDE", "milo_rna_nhood.pdf"),
  width = 10,
  height = 7
)

# Differential Expression Analysis ----
# Perform differential expression analysis using Milo
system.time(
  de_stat <- miloDE::de_test_neighbourhoods(
    milo_DE,
    sample_id = "sample",
    design = ~ 0 + diagnosis,
    covariates = c("diagnosis"),
    contrasts = c("diagnosisCIDP - diagnosisCTRL"),
    output_type = "SCE",
    plot_summary_stat = TRUE,
    layout = "UMAP.STACAS.SS.ALL",
    BPPARAM = NULL,
    verbose = TRUE,
    min_count = 3
  )
)

# Save Results ----
# Save the differential expression statistics
qs::qsave(de_stat, file.path("objects", "milo_de_stat_cidp_ctrl_csf.qs"))
# system("systemctl suspend")

stat_de_magnitude <- rank_neighbourhoods_by_DE_magnitude(de_stat)

p1 <- plot_milo_by_single_metric(
  milo_DE,
  stat_de_magnitude,
  colour_by = "n_DE_genes",
  layout = "UMAP.STACAS.SS.ALL",
  size_range = c(0.5, 5),
  edge_width = c(0.1, 1.0),
  edge_weight.thres = 10
) +
  viridis::scale_fill_viridis(name = "# DE genes", option = "inferno")

ggsave(
  plot = p1,
  filename = file.path("results", "miloDE", "milo_DE_CIDP_CTRL_CSF.pdf"),
  width = 6,
  height = 6,
  device = cairo_pdf
)
