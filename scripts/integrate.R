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
sc_merge <- split(x = sc_merge, f = sc_merge$sample)

sc_merge <- IntegrateLayers(
  object = sc_merge,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = TRUE,
)

sc_merge <- IntegrateLayers(
  object = sc_merge,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "rpca",
  verbose = TRUE,
  k.weight = 50,
)

set.seed(123)
sc_merge <- IntegrateLayers(
  object = sc_merge,
  method = scVIIntegration,
  orig.reduction = "pca",
  new.reduction = "scvi",
  conda_env = "~/miniconda3/envs/scvi",
  verbose = TRUE
)

# rejoin layers after integration 
sc_merge <- JoinLayers(sc_merge)

# stancas integration unsupervised ----
sc_merge <- SplitObject(sc_merge, split.by = "sample")
sc_merge_stacas <- STACAS::Run.STACAS(sc_merge)
sc_merge$stacas <- sc_merge_stacas$pca

# stancas integration semi-supervised based on azimuth annotation high confidence ---
sc_merge_stacas_ss_highconf <- sc_merge
pred_high_conf <-
    data.frame(pred = sc_merge$azimuth_pbmcref2, score = sc_merge$azimuth_pbmcref2_score) |>
    mutate(pred = ifelse(score < 0.99, "unknown", pred))
sc_merge_stacas_ss_highconf <- AddMetaData(sc_merge_stacas_ss_highconf, pred_high_conf)
sc_merge_stacas_ss_highconf <- SplitObject(sc_merge_stacas_ss_highconf, split.by = "sample")
sc_merge_stacas_ss_highconf <- STACAS::Run.STACAS(sc_merge_stacas_ss_highconf, cell.labels = "pred")

qsave(sc_merge_stacas_ss_highconf, file.path("objects", "sc_merge_stacas_ss_highconf_highconf.qs"))
sc_merge$stacas.ss.highconf <- sc_merge_stacas_ss_highconf$pca

# stancas integration semi-supervised based on azimuth annotation all ---
# use only one core to save memory
future::plan("multicore", workers = 1) 
options(future.globals.maxSize = 16000 * 1024^2)
sc_merge_stacas_ss_all <- sc_merge
sc_merge_stacas_ss_all <- SplitObject(sc_merge_stacas_ss_all, split.by = "sample")
sc_merge_stacas_ss_all <- STACAS::Run.STACAS(sc_merge_stacas_ss_all, cell.labels = "azimuth_pbmcref2", k.weight = 50)

qsave(sc_merge_stacas_ss_all, file.path("objects", "sc_merge_stacas_ss_all.qs"))
sc_merge$stacas.ss.all <- sc_merge_stacas_ss_all$pca

# umap ----
sc_merge <- RunUMAP(sc_merge, reduction = "harmony", reduction.name = "umap.harmony", dims = 1:30)
sc_merge <- RunUMAP(sc_merge, reduction = "rpca", reduction.name = "umap.rpca", dims = 1:30)
sc_merge <- RunUMAP(sc_merge, reduction = "integrated.scvi", reduction.name = "umap.scvi", dims = 1:30)
sc_merge <- RunUMAP(sc_merge, reduction = "stacas", reduction.name = "umap.stacas", dims = 1:30)
sc_merge <- RunUMAP(sc_merge, reduction = "stacas.ss.highconf", reduction.name = "umap.stacas.ss.highconf", dims = 1:30)
sc_merge <- RunUMAP(sc_merge, reduction = "stacas.ss.all", reduction.name = "umap.stacas.ss.all", dims = 1:30)

qsave(sc_merge, file.path("objects", "sc_merge.qs"))

# plot umaps per batch ----
my_cols_100 <- unname(Polychrome::createPalette(100, pals::cols25()))
integration_methods <- c("harmony", "rpca", "scvi", "stacas", "stacas.ss.highconf", "stacas.ss.all")
umap_methods <- paste0("umap.", integration_methods)

# per tissue
umap_plots_tissue <-
    lapply(
        umap_methods,
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
        umap_methods,
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
integration_methods <- c("harmony", "rpca", "scvi", "stacas", "stacas.ss.highconf", "stacas.ss.all")
umap_methods <- paste0("umap.", integration_methods)

pred_plots_pbmcref_level2 <-
    lapply(
        umap_methods,
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

integration_methods <- str_replace(integration_methods, "\\.", "_")
names(metrics_sample) <- str_replace_all(integration_methods, "\\.", "_")

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