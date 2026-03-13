###############################################################################
# perform gene set enrichment analysis of CD8_NK population
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
library(rlang)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

# get enrichment dbs ----
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

# Helper functions ----
get_top_genes <- function(
  file_path,
  sheet_name,
  n = 100,
  min_log2fc = 2,
  max_padj = 0.05
) {
  read_excel(file_path, sheet = sheet_name) |>
    dplyr::filter(
      p_val_adj < max_padj,
      avg_log2FC > min_log2fc
    ) |>
    dplyr::slice_min(
      order_by = p_val_adj,
      n = n,
      with_ties = FALSE
    )
}

perform_enrichment <- function(clusters, file_path, prefix = "") {
  # Get enrichment results
  enrichr_results <- lapply(clusters, function(x) {
    genes <- get_top_genes(file_path, x)$gene
    enrichr_res <- enrichr(genes, dbs)

    # Write results to Excel
    write_xlsx(
      enrichr_res,
      file.path(
        "results",
        "enrichr",
        paste0("enrichr_", prefix, x, ".xlsx")
      )
    )

    enrichr_res
  }) |>
    setNames(clusters)

  # Plot results
  plot_enrichment(clusters, prefix)

  return(enrichr_results)
}

plotEnrichrFun <- function(filename, sheet, width, height) {
  plot_data <- readxl::read_excel(
    file.path("results", "enrichr", glue::glue("enrichr_{filename}.xlsx")),
    sheet = sheet
  ) |>
    dplyr::slice_min(
      order_by = Adjusted.P.value,
      n = 10,
      with_ties = FALSE
    ) |>
    tidyr::separate(Overlap, into = c("overlap1", "overlap2")) |>
    dplyr::mutate(
      Term = gsub(
        x = Term,
        pattern = "\\s\\(.+\\)",
        replacement = ""
      ),
      overlap = as.numeric(overlap1) / as.numeric(overlap2)
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      y = reorder(Term, -log10(Adjusted.P.value)),
      x = -log10(Adjusted.P.value)
    )
  ) +
    ggplot2::geom_col(fill = scales::hue_pal()(5)[1]) +
    ggplot2::labs(
      x = "-Log10 Adjusted P value",
      y = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none")

  ggplot2::ggsave(
    file.path(
      "results",
      "enrichr",
      glue::glue("barplot_enrichr_{filename}_{sheet}.pdf")
    ),
    plot = p,
    width = width,
    height = height
  )
}

plot_enrichment <- function(
  clusters,
  prefix = "",
  sheet = "GO_Biological_Process_2023",
  width = 6,
  height = 2
) {
  filenames <- if (prefix == "") clusters else paste0(prefix, clusters)
  lapply(filenames, function(x) {
    plotEnrichrFun(x, sheet = sheet, width = width, height = height)
  })
}

# Main analysis ----
# Analysis parameters
cluster_sel <- c("CD8_NK")
de_cidp_ctrl_csf_clusters <- c("CD8_NK", "CD4TCM_2")

# Perform enrichment analysis for top markers
top_enrichr <- perform_enrichment(
  clusters = cluster_sel,
  file_path = file.path("results", "de", "topmarkers.xlsx"),
  prefix = "topmarkers_"
)

# Perform enrichment analysis for DE genes
de_cidp_ctrl_csf_enrichr_pos <- perform_enrichment(
  clusters = de_cidp_ctrl_csf_clusters,
  file_path = file.path("results", "de", "de_cidp_ctrl_csf.xlsx"),
  prefix = "de_cidp_ctrl_csf_pos_"
)
