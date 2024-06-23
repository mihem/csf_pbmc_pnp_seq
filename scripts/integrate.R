##################################################
# integrate data
# requires running azimuth.R first
##################################################
library(tidyverse)
library(Seurat)
library(Azimuth)
library(qs)
library(scMisc)
library(scIntegrationMetrics)
library(STACAS)
library(SeuratWrappers)

# general settings  ----
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)

# load data ----
# read preprocessed data from preprocess.R
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

# integrate data ---
# split for Seurat integration
# sc_merge <- split(x = sc_merge, f = sc_merge$sample)
sc_merge <- split(x = sc_merge, f = sc_merge$batch)

# use different approaches for integration 
# 1st level: different integration approaches
# 2nd level: PCA vs NMF
# 3rd level: batch at sample level vs multiplexed batches

# PCA harmony
sc_merge <- IntegrateLayers(
    object = sc_merge,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = TRUE,
)

# NMF harmony
sc_merge <- IntegrateLayers(
    object = sc_merge,
    method = HarmonyIntegration,
    orig.reduction = "nmf",
    new.reduction = "nmf.harmony",
    verbose = TRUE,
)


qs::qsave(sc_merge, file.path("objects", "sc_merge.qs"))

# scvi PCA
sc_merge <- IntegrateLayers(
    object = sc_merge,
    method = scVIIntegration,
    orig.reduction = "pca",
    new.reduction = "scvi",
    conda_env = "~/miniconda3/envs/scvi",
    verbose = TRUE
)

# scvi NMF
sc_merge <- IntegrateLayers(
    object = sc_merge,
    method = scVIIntegration,
    orig.reduction = "nmf",
    new.reduction = "nmf.scvi",
    conda_env = "~/miniconda3/envs/scvi",
    verbose = TRUE
)

# RPCA PCA
sc_merge <- IntegrateLayers(
    object = sc_merge,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "rpca",
    verbose = TRUE,
    k.weight = 50,
)

# rejoin layers after integration 
sc_merge <- JoinLayers(sc_merge)

# stancas integration unsupervised ----
# stacas requires split objects
# sc_merge <- SplitObject(sc_merge, split.by = "sample")
sc_merge_stacas <- SplitObject(sc_merge, split.by = "batch")
sc_merge_stacas <- STACAS::Run.STACAS(sc_merge)
sc_merge$stacas <- sc_merge_stacas$pca

sc_merge <- qread(file.path("objects", "sc_merge.qs"))

# stancas integration semi-supervised based on azimuth annotation all
sc_merge_stacas_ss_all <- sc_merge
rm(sc_merge)
# sc_merge_stacas_ss_all <- SplitObject(sc_merge_stacas_ss_all, split.by = "sample")
sc_merge_stacas_ss_all <- SplitObject(sc_merge_stacas_ss_all, split.by = "batch")
sc_merge_stacas_ss_all <- STACAS::Run.STACAS(sc_merge_stacas_ss_all, cell.labels = "azimuth_pbmcref2", k.weight = 50)

# qsave(sc_merge_stacas_ss_all, file.path("objects", "sc_merge_stacas_ss_all.qs"))
qsave(sc_merge_stacas_ss_all, file.path("objects", "sc_merge_stacas_ss_all.qs"))

sc_merge$stacas.ss.all <- sc_merge_stacas_ss_all$pca

# umap ----
# integration_methods <- c("pca", "nmf", "harmony", "nmf.harmony", "rpca", "scvi", "nmf.scvi", "stacas", "stacas.ss.highconf", "stacas.ss.all")
integration_methods <- c("pca", "nmf", "harmony", "nmf.harmony", "rpca", "scvi", "nmf.scvi", "stacas", "stacas.ss.all")

for (integration_method in integration_methods) {
  umap_name <- paste0("umap.", integration_method)
  if (!umap_name %in% names(sc_merge@reductions)) {
    sc_merge <- RunUMAP(sc_merge, reduction = integration_method, reduction.name = umap_name, dims = 1:30)
  } else {
    message(paste("UMAP reduction", umap_name, "already exists. Skipping..."))
  }
}

# plot umaps per batch ----
my_cols_100 <- unname(Polychrome::createPalette(100, pals::cols25()))

# per tissue
umap_plots_tissue <-
    lapply(
        paste0("umap.", integration_methods),
        FUN = function(x) {
            DimPlot(sc_merge, reduction = x, pt.size = .5, raster = FALSE, alpha = 0.1, group.by = "tissue", cols = my_cols_100) +
                theme_rect() +
                xlab("UMAP1") +
                ylab("UMAP2")
        }
    )

names(umap_plots_tissue) <- integration_methods

for (name in names(umap_plots_tissue)) {
    ggsave(file.path("results", "umap", paste0(name, "_umap_tissue", ".png")), plot = umap_plots_tissue[[name]], width = 20, height = 8)
}

# per sample
umap_plots_sample <-
    lapply(
        paste0("umap.", integration_methods),
        FUN = function(x) {
            DimPlot(sc_merge, reduction = x, pt.size = .5, raster = FALSE, alpha = 0.1, group.by = "sample", cols = my_cols_100) +
                theme_rect() +
                xlab("UMAP1") +
                ylab("UMAP2")
        }
    )

names(umap_plots_sample) <- integration_methods

for (name in names(umap_plots_sample)) {
    ggsave(file.path("results", "umap", paste0(name, "_umap_sample", ".png")), plot = umap_plots_sample[[name]], width = 20, height = 8)
}

# plot predictions ----
pred_plots_pbmcref_level2 <-
    lapply(
        paste0("umap.", integration_methods),
        FUN = function(x) {
            DimPlot(
                sc_merge,
                reduction = x,
                group.by = "azimuth_pbmcref2",
                raster = FALSE,
                pt.size = .1,
                alpha = .1,
                cols = my_cols_100,
                label = TRUE
            ) +
                theme_rect()
        }
    )

names(pred_plots_pbmcref_level2) <- integration_methods

for (name in names(pred_plots_pbmcref_level2)) {
    ggsave(file.path("results", "map", paste0("map_pbmcref_level2_", name, ".png")), plot = pred_plots_pbmcref_level2[[name]], width = 20, height = 8)
}

# integration metrics ----

# first subset because integration metrics do not scale well
Idents(sc_merge) <- sc_merge$azimuth_pbmcref2
sc_merge_small <- subset(sc_merge, downsample = 1000)

# for each sample
metrics_sample <-
    lapply(
        integration_methods,
        FUN = function(x) {
            scIntegrationMetrics::getIntegrationMetrics(
                sc_merge_small,
                meta.label = "azimuth_pbmcref2",
                meta.batch = "sample",
                method.reduction = x,
                metrics = c("CiLISI", "iLISI", "celltype_ASW")
            )
        }
    )

integration_methods_underscore <- str_replace(integration_methods, "\\.", "_")
names(metrics_sample) <- str_replace_all(integration_methods_underscore, "\\.", "_")

as.data.frame(metrics_sample) |>
    pivot_longer(everything()) |>
    separate_wider_delim(name, names = c("method", "metric"), delim = ".") |>
    pivot_wider(names_from = metric, values_from = value)  |>
    arrange(desc(celltype_ASW)) |>
    write_csv(file = file.path("results", "qc", "integration_metrics_sample.csv"))

# for each tissue
metrics_tissue <-
    lapply(
        integration_methods,
        FUN = function(x) {
            scIntegrationMetrics::getIntegrationMetrics(
                sc_merge_small,
                meta.label = "azimuth_pbmcref2",
                meta.batch = "tissue",
                method.reduction = x,
                metrics = c("CiLISI", "iLISI", "celltype_ASW")
            )
        }
    )

names(metrics_tissue) <- str_replace_all(integration_methods, "\\.", "_")

as.data.frame(metrics_tissue) |>
    pivot_longer(everything()) |>
    separate_wider_delim(name, names = c("method", "metric"), delim = ".") |>
    pivot_wider(names_from = metric, values_from = value)  |>
    arrange(desc(celltype_ASW)) |>
    write_csv(file = file.path("results", "qc", "integration_metrics_tissue.csv"))
