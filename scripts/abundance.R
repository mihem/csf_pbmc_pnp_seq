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

# source functions ----
source(file.path("scripts", "flow_helper.R"))

# detach(package:scMisc, unload = TRUE)
# remotes::install_github("mihem/scMisc")

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

sc_merge_csf <- subset(sc_merge, tissue %in% "CSF")
sc_merge_pbmc <- subset(sc_merge, tissue %in% "PBMC")

# meta data ----
lookup <- qs::qread(file.path("objects", "lookup.qs"))

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
  object = sc_merge_csf,
  x_axis = "sample",
  y_axis = "cluster",
  x_order = unique(sc_merge_csf$sample),
  y_order = sc_merge_csf@misc$cluster_order,
  color = sc_merge_csf@misc$cluster_col,
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
  object = sc_merge_csf,
  x_axis = "group",
  y_axis = "cluster",
  x_order = unique(sc_merge_csf$group),
  y_order = sc_merge_csf@misc$cluster_order,
  color = sc_merge_csf@misc$cluster_col,
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
  object = sc_merge_csf,
  x_axis = "diagnosis",
  y_axis = "cluster",
  x_order = unique(sc_merge_csf$diagnosis),
  y_order = sc_merge_csf@misc$cluster_order,
  color = sc_merge_csf@misc$cluster_col,
  width = 5
)

# create all combinations with defined order
# filter out same conditions and reverse combinations
seu_objects <- list(
  "CSF" = sc_merge_csf,
  "PBMC" = sc_merge_pbmc
)

# Create combinations for different groupings
combinations_group <- createCombinations(
  conditions = c("PNP", "CTRL"),
  group_column_name = "group",
  tissues = c("CSF", "PBMC")
)

combinations_group2 <- createCombinations(
  conditions = c("IN", "NIN", "CTRL"),
  group_column_name = "group2",
  tissues = c("CSF", "PBMC")
)

combinations_diagnosis <- createCombinations(
  conditions = c("CIDP", "GBS", "CIAP", "CTRL"),
  group_column_name = "diagnosis",
  tissues = c("CSF", "PBMC")
)

# Create a configuration table for all volcano plot comparisons
abundance_configs <- bind_rows(
  combinations_group,
  combinations_group2,
  combinations_diagnosis
)


# Calculate propeller results for all combinations
propeller_results <- lapply(seq_len(nrow(abundance_configs)), function(i) {
  formula_str <- paste0(
    "~0 + ",
    abundance_configs$group_column[i],
    " + sex + age"
  )
  scMisc::propellerCalc(
    seu_obj1 = seu_objects[[abundance_configs$tissue[i]]],
    condition1 = abundance_configs$condition1[i],
    condition2 = abundance_configs$condition2[i],
    cluster_col = "cluster",
    meta_col = abundance_configs$group_column[i],
    lookup = lookup,
    sample_col = "patient",
    formula = as.formula(formula_str),
    min_cells = 9
  )
})

# Plot the results
lapply(seq_along(propeller_results), function(i) {
  filename <- paste0(
    abundance_configs$tissue[i],
    "_",
    abundance_configs$condition1[i],
    "_",
    abundance_configs$condition2[i]
  )
  scMisc::plotPropeller(
    data = propeller_results[[i]],
    color = sc_merge@misc$cluster_col,
    filename = filename,
    FDR = 0.1
  )
})


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

# investigate subgroup of CIDP/GBS which has NKCD56bright_1 cells
sc_merge_csf_subgroup <- sc_merge_csf

nkbright_high_group <- c("P16", "P17", "P18", "P21", "P23")
sc_merge_csf_subgroup$subgroup <- ifelse(
  sc_merge_csf_subgroup$patient %in% nkbright_high_group,
  "NKhi",
  as.character(sc_merge_csf_subgroup$diagnosis)
)

sc_merge_csf_subgroup@misc$diagnosis_order <- c(
  "CTRL",
  "NKhi",
  "CIAP",
  "GBS",
  "CIDP",
  "MAG",
  "MFS",
  "PNC",
  "CAN",
  "PPN"
)

scMisc::abBoxPlot(
  object = sc_merge_csf_subgroup,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = sc_merge_csf_subgroup@misc$cluster_order,
  group_by = "subgroup",
  group_order = sc_merge_csf_subgroup@misc$diagnosis_order,
  color = sc_merge_csf_subgroup@misc$diagnosis_col,
  number_of_tests = choose(9, 2),
  width = 16
)

lookup_subgroup <- lookup |>
  dplyr::mutate(
    subgroup = dplyr::case_when(
      patient %in% nkbright_high_group ~ "NKhi",
      diagnosis %in% c("CIDP", "GBS") ~ "IN",
      TRUE ~ as.character(diagnosis)
    )
  )

sc_merge_csf_subgroup$subgroup <- dplyr::case_when(
  sc_merge_csf_subgroup$patient %in% nkbright_high_group ~ "NKhi",
  sc_merge_csf_subgroup$diagnosis %in% c("CIDP", "GBS") ~ "IN",
  TRUE ~ as.character(sc_merge_csf_subgroup$diagnosis)
)

propeller_results_subgroup_csf <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_csf_subgroup,
    condition1 = "NKhi",
    condition2 = "IN",
    cluster_col = "cluster",
    meta_col = "subgroup",
    lookup = lookup_subgroup,
    sample_col = "patient",
    formula = "~0 + subgroup + sex + age",
    min_cells = 9
  )

scMisc::plotPropeller(
  data = propeller_results_subgroup_csf,
  color = sc_merge@misc$cluster_col,
  filename = "CSF_NKhi_IN",
  FDR = 0.1
)
