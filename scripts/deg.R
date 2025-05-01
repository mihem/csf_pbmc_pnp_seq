##################################################
# differential gene expression analysis
# requires running annotate.R first
##################################################

# load libraries
library(tidyverse)
library(Seurat)
library(qs)
library(scMisc)
library(presto)
library(clustree)
library(readxl)
library(writexl)
library(cerebroAppLite)
library(Polychrome)

# load data ----
# read preprocessed data
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)


# subset seurat object ---
sc_merge_main_diagnosis <- subset(
    sc_merge,
    diagnosis %in% c("CIDP", "GBS", "CIAP", "CTRL")
)

sc_merge_main_diagnosis_pbmc <- subset(
    sc_merge_main_diagnosis,
    tissue %in% "PBMC"
)

sc_merge_main_diagnosis_pbmc_t_nk <- subset(
    sc_merge_main_diagnosis_pbmc,
    cluster %in% c("CD8_naive", "CD8TEM_1", "CD8TEM_2", "CD8_NK", "NKCD56bright", "NKCD56dim")
)

sc_merge_main_diagnosis_pbmc_cd8_nk <- subset(
    sc_merge_main_diagnosis_pbmc,
    cluster %in% c("CD8_NK")
)


levels(sc_merge$cluster)

bnb_markers <- read_csv(file.path("lookup", "markers.csv"))

#sc merge main
vln_bnb_plot <-
    VlnPlot(
        sc_merge_main_diagnosis,
        features = bnb_markers$BNB,
        group.by = "diagnosis",
        pt.size = .1
    ) +
    NoLegend() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(
    filename = file.path("results", "de", "vln_bnb_main.png"),
    plot = vln_bnb_plot,
    width = 10,
    height = 20
)

#pbmc cd8 nk
vln_bnb_pbmc_cd8_nk_plot <-
    VlnPlot(
        sc_merge_main_diagnosis_pbmc_cd8_nk,
        features = bnb_markers$BNB,
        group.by = "diagnosis",
        pt.size = .1
    ) +
    NoLegend() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(
    filename = file.path("results", "de", "vln_bnb_pbmc_cd8_nk.png"),
    plot = vln_bnb_pbmc_cd8_nk_plot,
    width = 10,
    height = 20
)

#pbmc cd8 nk
vln_bnb_pbmc_t_nk_plot <-
    VlnPlot(
        sc_merge_main_diagnosis_pbmc_t_nk,
        features = bnb_markers$BNB,
        group.by = "diagnosis",
        pt.size = .1
    ) +
    NoLegend() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(
    filename = file.path("results", "de", "vln_bnb_pbmc_t_nk.png"),
    plot = vln_bnb_pbmc_t_nk_plot,
    width = 10,
    height = 20
)

# sc merge main pbmc

# heatmap of BNB genes ----
main_bnb_avg <-
    AverageExpression(
        sc_merge_main_diagnosis,
        features = bnb_markers$BNB,
        group.by = "diagnosis"
    )$RNA

main_pbmc_bnb_avg <-
    AverageExpression(
        sc_merge_main_diagnosis_pbmc,
        features = bnb_markers$BNB,
        group.by = "diagnosis"
    )$RNA

scMisc::pHeatmap(main_bnb_avg, scale = "row")
scMisc::pHeatmap(main_pbmc_bnb_avg, scale = "row")
