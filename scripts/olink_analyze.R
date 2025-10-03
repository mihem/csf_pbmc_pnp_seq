##################################################
# Analyze Olink data
##################################################

# Load necessary libraries
library(tidyverse)
library(qs)
library(readxl)
library(tidyxl)
library(broom)
library(emmeans)
library(writexl)
library(patchwork)

# source helper files
source(file.path("scripts", "olink_analyze_helper.R"))

# Read the metadata file
olink_metadata_file <- file.path("lookup", "olink_lookup.xlsx")
olink_metadata <- read_xlsx(olink_metadata_file)

# colors
color_olink_diagnosis <- setNames(
    pals::cols25(3),
    c("CTRL", "GBS", "CIDP")
)

# Read the Olink quantified data
olink_quant_file <- file.path(
    "raw",
    "olink",
    "P049E085_Müller-Miny_QUANT_long.xlsx"
)

olink_npx_file <- file.path(
    "raw",
    "olink",
    "P049E085_Müller-Miny_NPX_long.xlsx"
)

olink_quant <- read_xlsx(olink_quant_file) |>
    mutate(Quantified_value = as.numeric(Quantified_value))

olink_npx <- read_xlsx(olink_npx_file) |>
    mutate(NPX = as.numeric(NPX))


# Filter red cells for both datasets
olink_quant_filtered <- filterRedCells(
    olink_quant,
    olink_quant_file,
    "Quantified_value"
)

olink_npx_filtered <- filterRedCells(
    olink_npx,
    olink_npx_file,
    "NPX"
)

# Analyze QUANT data filtered
results_quant <- processOlinkData(
    data = olink_quant_filtered,
    value_col = "Quantified_value",
    unit_label = "pg/ml"
)

qs::qsave(results_quant, file.path("objects", "olink_quant.qs"))

write_xlsx(
    results_quant$stats,
    file.path("results", "olink", "olink_quant_stats.xlsx")
)

ggsave(
    file.path("results", "olink", "olink_quant_boxplots.pdf"),
    plot = results_quant_filtered$boxplots,
    width = 10,
    height = 12
)

##################################################
# Analyze NPX data
##################################################

results_npx <- processOlinkData(
    data = olink_npx_filtered,
    value_col = "NPX",
    unit_label = "NPX"
)

# Save NPX results
qs::qsave(results_npx$data_wide, file.path("objects", "olink_npx_wide.qs"))

write_xlsx(
    results_npx$stats,
    file.path("results", "olink", "olink_npx_stats_all.xlsx")
)

ggsave(
    file.path("results", "olink", "olink_npx_boxplots.pdf"),
    plot = results_npx$boxplots,
    width = 10,
    height = 12
)
