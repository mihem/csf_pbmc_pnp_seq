##################################################
# annotate data
# requires running annotate.R first
##################################################
#libraries ---
library(Seurat)
library(tidyverse)
library(scMisc)
library(qs)
library(pheatmap)
library(speckle)

# detach(package:scMisc, unload = TRUE)
# remotes::install_github("mihem/scMisc")

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

sc_merge_csf <- subset(sc_merge, tissue %in% "CSF")
sc_merge_pbmc <- subset(sc_merge, tissue %in% "PBMC")

# meta data ----
lookup <-
  readxl::read_excel(file.path("lookup", "SEED_lookup_v6.xlsx")) |>
  janitor::clean_names() |>
  mutate(age = lubridate::time_length(difftime(date, birth_date), "years")) |>
  mutate(
    diagnosis = factor(
      diagnosis,
      levels = c(
        "CTRL",
        "CIAP",
        "CIDP",
        "GBS",
        "MAG",
        "MFS",
        "PNC",
        "CAN",
        "PPN"
      )
    )
  )

# abundance table
scMisc::abundanceTbl(sc_merge, "cluster", "sample")
scMisc::abundanceTbl(sc_merge, "cluster", "tissue_diagnosis")
scMisc::abundanceTbl(sc_merge, "cluster", "tissue_group")

# stacked plots --
scMisc::stackedPlot(
  object = sc_merge_csf,
  x_axis = "sample",
  y_axis = "cluster",
  x_order = unique(sc_merge_csf$sample),
  y_order = sc_merge_csf@misc$cluster_order,
  color = sc_merge_csf@misc$cluster_col,
  width = 7
)

scMisc::stackedPlot(
  object = sc_merge_pbmc,
  x_axis = "sample",
  y_axis = "cluster",
  x_order = unique(sc_merge_pbmc$sample),
  y_order = sc_merge_pbmc@misc$cluster_order,
  color = sc_merge_pbmc@misc$cluster_col,
  width = 7
)

scMisc::stackedPlot(
  object = sc_merge_csf,
  x_axis = "group",
  y_axis = "cluster",
  x_order = unique(sc_merge_csf$group),
  y_order = sc_merge_csf@misc$cluster_order,
  color = sc_merge_csf@misc$cluster_col,
  width = 3.8
)

scMisc::stackedPlot(
  object = sc_merge_pbmc,
  x_axis = "group",
  y_axis = "cluster",
  x_order = unique(sc_merge_pbmc$group),
  y_order = sc_merge_pbmc@misc$cluster_order,
  color = sc_merge_pbmc@misc$cluster_col,
  width = 3.8
)

scMisc::stackedPlot(
  object = sc_merge_csf,
  x_axis = "diagnosis",
  y_axis = "cluster",
  x_order = unique(sc_merge_csf$diagnosis),
  y_order = sc_merge_csf@misc$cluster_order,
  color = sc_merge_csf@misc$cluster_col,
  width = 5
)

scMisc::stackedPlot(
  object = sc_merge_pbmc,
  x_axis = "diagnosis",
  y_axis = "cluster",
  x_order = unique(sc_merge_pbmc$diagnosis),
  y_order = sc_merge_pbmc@misc$cluster_order,
  color = sc_merge_pbmc@misc$cluster_col,
  width = 5
)

# propeller abundance analysis ---
# PNP vs CTRL CSF
propeller_CSF_PNP_CTRL <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_csf,
    condition1 = "PNP",
    condition2 = "CTRL",
    cluster_col = "cluster",
    meta_col = "group",
    lookup = lookup,
    sample_col = "patient",
    formula = "~0 + group",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_CSF_PNP_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "CSF_PNP_CTRL",
  FDR = 0.1
)

# scMisc::dotplotPropeller(
#     data = propeller_PNP_CTRL,
#     color = sc_merge@misc$cluster_col,
#     filename = "PNP_CTRL",
# )

# # only plot logFC > 0.5
# propeller_PNP_CTRL |>
#   dplyr::filter(abs(log2ratio) > 0.5)  |>
#   scMisc::dotplotPropeller(
#     data = _,
#     color = sc_merge@misc$cluster_col,
#     filename = "PNP_CTRL_logFC_0.5",
#     width = 2.5,
#     height = 3
#   )

# PNP vs CTRL PBMC
propeller_PBMC_PNP_CTRL <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_pbmc,
    condition1 = "PNP",
    condition2 = "CTRL",
    cluster_col = "cluster",
    meta_col = "group",
    lookup = lookup,
    sample_col = "patient",
    formula = "~0 + group",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_PBMC_PNP_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "PBMC_PNP_CTRL",
  FDR = 0.1
)

# create all combinations with defined order
# filter out same conditions and reverse combinations
diagnosis <- factor(
  c("CIDP", "GBS", "CIAP", "CTRL"),
  levels = c("CIDP", "GBS", "CIAP", "CTRL")
)

seu_objects <- list(
  "CSF" = sc_merge_csf,
  "PBMC" = sc_merge_pbmc
)

combinations <- crossing(
  tissue = names(seu_objects),
  condition1 = diagnosis,
  condition2 = diagnosis
) |>
  dplyr::filter(
    condition1 != condition2,
    as.numeric(condition1) < as.numeric(condition2)
  ) |>
  dplyr::mutate(across(everything(), as.character))

# Calculate propeller results for all combinations
propeller_results <-
  pmap(
    combinations,
    function(tissue, condition1, condition2) {
      scMisc::propellerCalc(
        seu_obj1 = seu_objects[[tissue]],
        condition1 = condition1,
        condition2 = condition2,
        cluster_col = "cluster",
        meta_col = "diagnosis",
        lookup = lookup,
        sample_col = "patient",
        formula = "~0 + diagnosis + sex + age",
        min_cells = 30
      )
    }
  )

# Create filenames
filenames <- paste0(
  combinations$tissue, "_", 
  combinations$condition1, "_", 
  combinations$condition2
)

# Plot the results
map2(
  propeller_results,
  filenames,
  function(data, filename) {
    scMisc::plotPropeller(
      data = data,
      color = sc_merge@misc$cluster_col,
      filename = filename,
      FDR = 0.1
    )
  }
)

# abundance box plot ----
scMisc::abBoxPlot(
  object = sc_merge_csf,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = sc_merge_csf@misc$cluster_order,
  group_by = "diagnosis",
  group_order = sc_merge_csf@misc$diagnosis_order,
  color = sc_merge_csf@misc$diagnosis_col,
  number_of_tests = choose(9, 2),
  width = 14
)

scMisc::abBoxPlot(
  object = sc_merge_pbmc,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = sc_merge_pbmc@misc$cluster_order,
  group_by = "diagnosis",
  group_order = sc_merge_pbmc@misc$diagnosis_order,
  color = sc_merge_pbmc@misc$diagnosis_col,
  number_of_tests = choose(9, 2),
  width = 14
)
