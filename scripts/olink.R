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
excel_path <- file.path("results", "enrich", "de_combined_gsea_results.xlsx")

# Read all sheets using the function
gsea_results <- read_all_excel_sheets(excel_path)

# get genes of the top 10 enriched pathways and split into individual genes
genes <- gsea_results$gbs_ctrl_pbm$core_enrichment[1:10] |>
    strsplit(split = "/") |>
    unlist() |>
    unique() |>
    sort()

write_csv(
    data.frame(gene = intersect(genes, olink_flex_list$Gene)),
    file.path("results", "enrich", "de_combined_gsea_olink_intersect.csv")
)
