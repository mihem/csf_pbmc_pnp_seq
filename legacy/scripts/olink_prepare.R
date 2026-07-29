##################################################
# Find genes of interest in Olink panel
# requires running annotate.R first
##################################################

library(tidyverse)
library(readxl)

# compare with olink data
olink_flex_list <- read_csv2(file.path("lookup", "Olink Flex - 2025-05-13.csv"))

read_all_excel_sheets <- function(file_path) {
  # Get all sheet names
  sheets <- excel_sheets(file_path)

  # Read all sheets into a named list
  results <- lapply(sheets, function(sheet) {
    read_xlsx(file_path, sheet = sheet)
  })
  names(results) <- sheets

  return(results)
}

# Get the path to the Excel file
gsea_path <- file.path(
  "results",
  "enrich",
  "de_combined",
  "gsea",
  "de_combined_gsea_results.xlsx"
)

# Read all sheets using the function
gsea_results <- read_all_excel_sheets(gsea_path)

# get genes of the top 10 enriched pathways and split into individual genes
gsea_genes <- gsea_results$gbs_ctrl_pbmc$core_enrichment[1:10] |>
  strsplit(split = "/") |>
  unlist() |>
  unique() |>
  sort()

write_csv(
  data.frame(gene = intersect(gsea_genes, olink_flex_list$Gene)),
  file.path("results", "enrich", "de_combined_gsea_olink_intersect.csv")
)

# intersect of combined DE genes with olink genes
olink_intersect <- function(file_name) {
  result <- read_xlsx(file.path("results", "de", file_name)) |>
    dplyr::filter(
      p_val_adj < 0.05,
      abs(avg_log2FC) > 1
    ) |>
    dplyr::filter(gene %in% olink_flex_list$Gene)
  if (nrow(result) != 0) {
    write_xlsx(
      result,
      file.path("results", "de", paste0("olink_", file_name))
    )
  } else {
    message("No significant genes found in ", file_name)
  }
}

file_names <- c(
  "de_cidp_ctrl_csf_combined.xlsx",
  "de_gbs_ctrl_csf_combined.xlsx",
  "de_cidp_ctrl_pbmc_combined.xlsx",
  "de_gbs_ctrl_pbmc_combined.xlsx"
)

lapply(
  file_names,
  olink_intersect
)
