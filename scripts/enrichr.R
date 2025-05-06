###############################################################################
# perform gene set enrichment annalysis of CD8_NK population
# requires running annotate.R first
###############################################################################

# load libraries ----
library(enrichR)
library(qs)
library(Seurat)
library(tidyverse)
library(patchwork)
library(readxl)
library(writexl)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

# get enrichment dbs ----
dbs <- enrichR::listEnrichrDbs()
dbs <- c(
    "TF_Perturbations_Followed_by_Expression",
    "Transcription_Factor_PPIs",
    "WikiPathways_2024_Human",
    "KEGG_2021_Human",
    "Reactome_2022",
    "Panther_2016",
    "NCI-Nature_2016",
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023"
)

# get top markers ----
cluster_sel <- c("CD8_NK")
top <-
    lapply(c(cluster_sel), FUN = function(x) {
        read_excel(file.path("results", "de", "topmarkers.xlsx"), sheet = x) |>
            dplyr::filter(
                p_val_adj < 0.05
            ) |>
            dplyr::filter(
                avg_log2FC > 2
            ) |>
            dplyr::slice_max(
                order_by = avg_log2FC,
                n = 100,
                with_ties = FALSE
            )
    }) |>
    setNames(c(cluster_sel))

# get de markers ----
de_cidp_ctrl_csf_clusters <- c("CD8_NK", "CD4TCM_2")

de_cidp_ctrl_csf_pos <-
    lapply(c(de_cidp_ctrl_csf_clusters), FUN = function(x) {
        readxl::read_excel(
            file.path("results", "de", "de_cidp_ctrl_csf.xlsx"),
            sheet = x
        ) |>
            dplyr::filter(
                p_val_adj < 0.05
            ) |>
            dplyr::filter(
                avg_log2FC > 2
            ) |>
            dplyr::slice_max(
                order_by = avg_log2FC,
                n = 100,
                with_ties = FALSE
            )
    }) |>
    setNames(c(de_cidp_ctrl_csf_clusters))

# perform enrichment with selected clusters of topmarkers ----
top_enrichr <-
    lapply(c(cluster_sel), FUN = function(x) {
        enrichr(top[[x]]$gene, dbs)
    }) |>
    setNames(c(cluster_sel))

lapply(
    c(cluster_sel),
    FUN = function(x) {
        write_xlsx(
            top_enrichr[[x]],
            file.path("results", "enrichr", paste0("enrichr_", x, ".xlsx"))
        )
    }
)

# perform enrichment with selected clusters of de_cidp_ctrl_csf ----
de_cidp_ctrl_csf_enrichr_pos <-
    lapply(c(de_cidp_ctrl_csf_clusters), FUN = function(x) {
        enrichr(de_cidp_ctrl_csf_pos[[x]]$gene, dbs)
    }) |>
    setNames(c(de_cidp_ctrl_csf_clusters))

lapply(
    de_cidp_ctrl_csf_clusters,
    FUN = function(x) {
        write_xlsx(
            de_cidp_ctrl_csf_enrichr_pos[[x]],
            file.path(
                "results",
                "enrichr",
                paste0("enrichr_de_cidp_ctrl_csf_pos_", x, ".xlsx")
            )
        )
    }
)

# function to plot enrichment results ----
plotEnrichrFun <- function(filename, sheet, width, height) {
    colors <- RColorBrewer::brewer.pal(5, "Set2")
    enrichr <- readxl::read_excel(
        file.path("results", "enrichr", glue::glue("enrichr_{filename}.xlsx")),
        sheet = sheet
    ) |>
        dplyr::slice_min(
            order_by = Adjusted.P.value,
            n = 10,
            with_ties = FALSE
        ) |>
        tidyr::separate(Overlap, into = c("overlap1", "overlap2")) |> # separate overlap in two columns
        dplyr::mutate(
            Term = gsub(x = Term, pattern = "\\s\\(.+\\)", replacement = "")
        ) |>
        dplyr::mutate(overlap = as.numeric(overlap1) / as.numeric(overlap2)) |> # calculcate overlap
        ggplot(aes(
            y = reorder(Term, -log10(Adjusted.P.value)),
            x = -log10(Adjusted.P.value)
        )) +
        geom_col(fill = scales::hue_pal()(5)[1]) +
        labs(
            x = "-Log10 Adjusted P value",
            y = ""
        ) +
        theme_classic() +
        theme(legend.position = "none")
    ggsave(
        file.path(
            "results",
            "enrichr",
            glue::glue("barplot_enrichr_{filename}_{sheet}.pdf")
        ),
        width = width,
        height = height
    )
}

# plot enrichment of topmarkers ----
lapply(
    c(cluster_sel),
    FUN = function(x) {
        plotEnrichrFun(
            x,
            sheet = "GO_Biological_Process_2023",
            width = 6,
            height = 2
        )
    }
)

# plot enrichment of de_cidp_ctrl_csf ----
enrichr_filenames_de_cidp_csf_pos <- c(
    "de_cidp_ctrl_csf_pos_CD8_NK",
    "de_cidp_ctrl_csf_pos_CD4TCM_2"
)

lapply(
    c(enrichr_filenames_de_cidp_csf_pos),
    FUN = function(x) {
        plotEnrichrFun(
            x,
            sheet = "GO_Biological_Process_2023",
            width = 6,
            height = 2
        )
    }
)
