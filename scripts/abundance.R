##################################################
# abundancce of scRNA clusters
# requires running annotate.R first
##################################################
#libraries ---
library(Seurat)
library(tidyverse)
library(scMisc)
library(qs)
library(pheatmap)
library(speckle)
library(writexl)
library(permFDP)

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
    min_cells = 9,
    adjustment_method = "permFDP",
    fdr_threshold = 0.1,
    n_perms = 1000
  )
})


# Create named list of propeller results
names(propeller_results) <- sapply(seq_along(propeller_results), function(i) {
  sheet_name <- paste0(
    abundance_configs$tissue[i],
    "_",
    abundance_configs$condition1[i],
    "_vs_",
    abundance_configs$condition2[i]
  )
  # Truncate sheet name to 31 characters (Excel limit)
  substr(sheet_name, 1, 31)
})

# Save to Excel
writexl::write_xlsx(
  propeller_results,
  path = file.path("results", "abundance", "propeller_results_all.xlsx")
)

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
    use_permFDP = TRUE,
    dir_output = file.path("results", "abundance", "propeller_plots"),
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
