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

# get de of macro ----
cluster_sel <- c("CD8_NK")

top <-
    lapply(c(cluster_sel), FUN = function(x) {
        read_excel(file.path("results", "de", "topmarkers.xlsx"), sheet = x) |>
            dplyr::slice_max(
                order_by = avg_log2FC,
                n = 100,
                with_ties = FALSE
            )
    }) |>
    setNames(c(cluster_sel))

# perform enrichment with macro ----
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

# plot enrichment of macro ----
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
