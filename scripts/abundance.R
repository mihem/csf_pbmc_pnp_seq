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
  mutate(diagnosis = factor(diagnosis, levels = c("CTRL", "CIAP", "CIDP", "GBS", "MAG", "MFS", "PNC", "CAN", "PPN")))

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

# CIDP vs CTRL CSF
propeller_CSF_CIDP_CTRL <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_csf,
    condition1 = "CIDP",
    condition2 = "CTRL",
    cluster_col = "cluster",
    meta_col = "diagnosis",
    lookup = lookup,
    sample_col = "patient",
    formula = "~0 + diagnosis",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_CSF_CIDP_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "CSF_CIDP_CTRL",
  FDR = 0.1
)

# CIDP vs CTRL PBMC
propeller_PBMC_CIDP_CTRL <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_pbmc,
    condition1 = "CIDP",
    condition2 = "CTRL",
    cluster_col = "cluster",
    meta_col = "diagnosis",
    lookup = lookup_PBMC,
    sample_col = "patient",
    formula = "~0 + diagnosis",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_PBMC_CIDP_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "PBMC_CIDP_CTRL",
  FDR = 0.1
)

# GBS vs CTRL CSF
propeller_CSF_GBS_CTRL <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_csf,
    condition1 = "GBS",
    condition2 = "CTRL",
    cluster_col = "cluster",
    meta_col = "diagnosis",
    lookup = lookup_CSF,
    sample_col = "patient",
    formula = "~0 + diagnosis",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_CSF_GBS_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "CSF_GBS_CTRL",
  FDR = 0.1
)

# GBS vs CTRL PBMC
propeller_PBMC_GBS_CTRL <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_pbmc,
    condition1 = "GBS",
    condition2 = "CTRL",
    cluster_col = "cluster",
    meta_col = "diagnosis",
    lookup = lookup_PBMC,
    sample_col = "patient",
    formula = "~0 + diagnosis",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_PBMC_GBS_CTRL,
  color = sc_merge@misc$cluster_col,
  filename = "PBMC_GBS_CTRL",
  FDR = 0.1
)

# CIDP vs CIAP CSF
propeller_CSF_CIDP_CIAP <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_csf,
    condition1 = "CIDP",
    condition2 = "CIAP",
    cluster_col = "cluster",
    meta_col = "diagnosis",
    lookup = lookup_CSF,
    sample_col = "patient",
    formula = "~0 + tissue_diagnosis",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_CSF_CIDP_CIAP,
  color = sc_merge@misc$cluster_col,
  filename = "CSF_CIDP_CIAP",
  FDR = 0.1
)

# CIDP vs CIAP PBMC
propeller_PBMC_CIDP_CIAP <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_pbmc,
    condition1 = "CIDP",
    condition2 = "CIAP",
    cluster_col = "cluster",
    meta_col = "diagnosis",
    lookup = lookup_PBMC,
    sample_col = "patient",
    formula = "~0 + diagnosis",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_PBMC_CIDP_CIAP,
  color = sc_merge@misc$cluster_col,
  filename = "PBMC_CIDP_CIAP",
  FDR = 0.1
)

# GBS vs CIAP CSF
propeller_CSF_GBS_CIAP <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_csf,
    condition1 = "GBS",
    condition2 = "CIAP",
    cluster_col = "cluster",
    meta_col = "diagnosis",
    lookup = lookup_CSF,
    sample_col = "patient",
    formula = "~0 + diagnosis",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_CSF_GBS_CIAP,
  color = sc_merge@misc$cluster_col,
  filename = "CSF_GBS_CIAP",
  FDR = 0.1
)

# GBS vs CIAP PBMC
propeller_PBMC_GBS_CIAP <-
  scMisc::propellerCalc(
    seu_obj1 = sc_merge_pbmc,
    condition1 = "GBS",
    condition2 = "CIAP",
    cluster_col = "cluster",
    meta_col = "diagnosis",
    lookup = lookup_PBMC,
    sample_col = "patient",
    formula = "~0 + diagnosis",
    min_cells = 30
  )

scMisc::plotPropeller(
  data = propeller_PBMC_GBS_CIAP,
  color = sc_merge@misc$cluster_col,
  filename = "PBMC_GBS_CIAP",
  FDR = 0.1
)


unique(sc_merge_csf$tissue_diagnosis)
unique(sc_merge_csf$tissue)

scMisc::abBoxPlot(
  object = sc_merge_csf,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = sc_merge_csf@misc$cluster_order,
  group_by = "diagnosis",
  group_order = sc_merge_csf@misc$diagnosis_order,
  color = sc_merge_csf@misc$diagnosis_col,
  number_of_tests = choose(9, 2))

# abundance box plot ----
scMisc::abBoxPlot(
  object = sc_merge_pbmc,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = sc_merge_pbmc@misc$cluster_order,
  group_by = "diagnosis",
  group_order = sc_merge_pbmc@misc$diagnosis_order,
  color = sc_merge_pbmc@misc$diagnosis_col,
  number_of_tests = choose(9, 2))


# abundance of mrVI groups ----
mrvi_lookup <- read_csv(file.path("lookup", "mrvi_lookup.csv"))

sc_merge@meta.data <-
    sc_merge@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(mrvi_lookup, by = "sample") |>
    tibble::column_to_rownames(var = "barcode")

vn_cidp_ciap_ctrl <- subset(sc_merge, level2 %in% c("VN", "CIDP", "CIAP", "CTRL"))

scMisc::abBoxPlot(
  object = vn_cidp_ciap_ctrl,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = vn_cidp_ciap_ctrl@misc$cluster_order,
  group_by =  "mrvi_cluster",
  group_order = paste0("cl", 1:5),
  color = pals::cols25()
)

phmap_mrvi_cluster <-
  table(vn_cidp_ciap_ctrl$cluster, vn_cidp_ciap_ctrl$mrvi_cluster) |>
  pheatmap(
    scale = "column",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = viridis::magma(100),
    cellwidth = 10,
    cellheight = 10,
    treeheight_row = 15,
    treeheight_col = 15,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    border_color = NA,
    cutree_rows = 5,
    main = "mrVI cluster"
  )

pdf(file.path("results", "abundance", "mrvi_heatmap_abundance_cluster.pdf"), width = 5, height = 7)
print(phmap_mrvi_cluster)
dev.off()


# abundance of mrVI groups in immune cells  ----
mrvi_lookup <- read_csv(file.path("lookup", "mrvi_lookup.csv"))  

ic@meta.data <-
    ic@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(mrvi_lookup, by = "sample") |>
    tibble::column_to_rownames(var = "barcode")

ic_vn_cidp_ciap_ctrl <- subset(ic, level2 %in% c("VN", "CIDP", "CIAP", "CTRL"))

str(ic_vn_cidp_ciap_ctrl@meta.data)

scMisc::abBoxPlot(
  object = ic_vn_cidp_ciap_ctrl,
  cluster_idents = "ic_cluster",
  sample = "sample",
  cluster_order = ic_vn_cidp_ciap_ctrl@misc$ic_cluster_order,
  group_by =  "mrvi_cluster",
  group_order = paste0("p-cl", 1:5),
  color = pals::cols25()
)

phmap_mrvi_cluster_ic <-
  table(ic_vn_cidp_ciap_ctrl$ic_cluster, ic_vn_cidp_ciap_ctrl$mrvi_cluster) |>
  pheatmap(
    scale = "column",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = viridis::magma(100),
    cellwidth = 10,
    cellheight = 10,
    treeheight_row = 15,
    treeheight_col = 15,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    border_color = NA,
    cutree_rows = 5,
    main = "mrVI cluster"
  )

pdf(file.path("results", "abundance", "mrvi_heatmap_abundance_ic_cluster.pdf"), width = 5, height = 7)
print(phmap_mrvi_cluster_ic)
dev.off()


# g ratio and axon diameter boxplots ----
g_ratio_axon_diameter_mrvi <-
  g_ratio |>
  group_by(sample) |>
  mutate(
    g_ratio = mean(g_ratio),
    axon_diameter = mean(axon_diameter),
  ) |>
  ungroup() |>
  distinct() |>
  left_join(mrvi_lookup, join_by(sample)) |>
  dplyr::filter(!is.na(mrvi_cluster))  

g_ratio_mrvi_stats <- scMisc:::compStat(x_var = "g_ratio", group = "mrvi_cluster", data = g_ratio_axon_diameter_mrvi, paired = FALSE)

g_ratio_mrvi_plot <-
  g_ratio_axon_diameter_mrvi |>
  ggplot(aes(x = mrvi_cluster, y = g_ratio)) +
  ggsignif::geom_signif(comparisons = g_ratio_mrvi_stats$comparisons, annotation = g_ratio_mrvi_stats$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7) +
  geom_boxplot(aes(fill = mrvi_cluster)) +
  geom_point() +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("g-ratio") +
  scale_fill_manual(values = pals::cols25()) +
  theme(legend.position = "none")
ggsave(file.path("results", "abundance", "boxplot_g_ratio_mrvi.pdf"), plot = g_ratio_mrvi_plot, width = 3, height = 3)

# axon diameter
axon_diameter_mrvi_stats <- scMisc:::compStat(x_var = "axon_diameter", group = "mrvi_cluster", data = g_ratio_axon_diameter_mrvi, paired = FALSE)

axon_diameter_mrvi_plot <-
  g_ratio_axon_diameter_mrvi |>
  ggplot(aes(x = mrvi_cluster, y = axon_diameter)) +
  ggsignif::geom_signif(comparisons = axon_diameter_mrvi_stats$comparisons, annotation = axon_diameter_mrvi_stats$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)  +
  geom_boxplot(aes(fill = mrvi_cluster)) +
  geom_point() +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("Axon diameter") +
  scale_fill_manual(values = pals::cols25()) +
  theme(legend.position = "none")
  
ggsave(file.path("results", "abundance", "boxplot_axon_diameter_mrvi.pdf"), plot = axon_diameter_mrvi_plot, width = 3, height = 3)

# axon counts
axon_count_mrvi <-
  axon_count_mean |>
  left_join(mrvi_lookup, join_by(sample)) |>
  dplyr::filter(!is.na(mrvi_cluster))

axon_count_mrvi_stats <- scMisc:::compStat(x_var = "log_axon_normal", group = "mrvi_cluster", data = axon_count_mrvi, paired = FALSE)

axon_count_mrvi_plot <-
  axon_count_mrvi |>
  ggplot(aes(x = mrvi_cluster, y = axon_normal)) +
  ggsignif::geom_signif(comparisons = axon_count_mrvi_stats$comparisons, annotation = axon_count_mrvi_stats$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)  +
  geom_boxplot(aes(fill = mrvi_cluster)) +
  geom_point() +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("Normal axon count") +
  scale_fill_manual(values = pals::cols25()) +
  theme(legend.position = "none")
  
ggsave(file.path("results", "abundance", "boxplot_axon_count_mrvi.pdf"), plot = axon_count_mrvi_plot, width = 3, height = 3)

# incat -----
incat_mrvi <-
  sample_lookup |>
  left_join(mrvi_lookup, join_by(sample)) |>
  dplyr::filter(!is.na(mrvi_cluster)) |>
  mutate(incat = as.numeric(incat))

incat_mrvi_stats <- scMisc:::compStat(x_var = "incat", group = "mrvi_cluster", data = incat_mrvi, paired = FALSE)

incat_mrvi_plot <-
  incat_mrvi |>
  ggplot(aes(x = mrvi_cluster, y = incat)) +
  ggsignif::geom_signif(comparisons = incat_mrvi_stats$comparisons, annotation = incat_mrvi_stats$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)  +
  geom_boxplot(aes(fill = mrvi_cluster)) +
  geom_jitter(height = 0, width = 0.1) +
  theme_bw() +
  xlab("") +
  ylab("") +
  ggtitle("INCAT score") +
  scale_fill_manual(values = pals::cols25()) +
  theme(legend.position = "none")
  
ggsave(file.path("results", "abundance", "boxplot_incat_mrvi.pdf"), plot = incat_mrvi_plot, width = 3, height = 3)

 renv::install(
  "stringi@1.8.4",
  rebuild = TRUE,
  repos = "https://packagemanager.rstudio.com/cran/latest"
)
