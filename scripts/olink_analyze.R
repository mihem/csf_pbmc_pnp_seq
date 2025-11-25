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
library(lme4)

# source helper files
source(file.path("scripts", "olink_analyze_helper.R"))

# Read the metadata file
olink_metadata_file <- file.path("lookup", "olink_flow_lookup.xlsx")
olink_metadata <- read_xlsx(olink_metadata_file)

# Define group configurations
group_configs <- list(
    diagnosis = list(
        levels = c("CTRL", "GBS", "CIDP"),
        colors = setNames(pals::cols25(3), c("CTRL", "GBS", "CIDP")),
        suffix = ""
    ),
    group2 = list(
        levels = c("CTRL", "IN"),
        colors = setNames(
            pals::cols25(2),
            c("CTRL", "IN")
        ), 
        suffix = "_meta"
    )
)

# Read the Olink quantified data
olink_quant_file <- file.path(
    "raw",
    "olink",
    "olink_quant_long_filtered.xlsx"
)

olink_npx_file <- file.path(
    "raw",
    "olink",
    "olink_npx_long_filtered.xlsx"
)

olink_quant <- read_xlsx(olink_quant_file) |>
    mutate(Quantified_value = as.numeric(Quantified_value))

olink_npx <- read_xlsx(olink_npx_file) |>
    mutate(NPX = as.numeric(NPX))

# Loop through both group and meta_group
for (group_name in names(group_configs)) {
    config <- group_configs[[group_name]]

    # Analyze QUANT data
    results_quant <- processOlinkData(
        data = olink_quant,
        value_col = "Quantified_value",
        unit_label = "pg/ml",
        group_var = group_name,
        group_levels = config$levels,
        colors = config$colors
    )

    qs::qsave(
        results_quant,
        file.path("objects", paste0("olink_quant", config$suffix, ".qs"))
    )

    write_xlsx(
        results_quant$stats,
        file.path(
            "results",
            "olink",
            paste0("olink_quant_stats_lme4", config$suffix, ".xlsx")
        )
    )

    ggsave(
        file.path(
            "results",
            "olink",
            paste0("olink_quant_boxplot_lme4", config$suffix, ".pdf")
        ),
        plot = results_quant$boxplots,
        width = 8,
        height = 12
    )

    # Analyze NPX data
    results_npx <- processOlinkData(
        data = olink_npx,
        value_col = "NPX",
        unit_label = "NPX",
        group_var = group_name,
        group_levels = config$levels,
        colors = config$colors
    )

    qs::qsave(
        results_npx$data_wide,
        file.path("objects", paste0("olink_npx_wide", config$suffix, ".qs"))
    )

    write_xlsx(
        results_npx$stats,
        file.path(
            "results",
            "olink",
            paste0("olink_npx_stats_all", config$suffix, ".xlsx")
        )
    )

    ggsave(
        file.path(
            "results",
            "olink",
            paste0("olink_npx_boxplots", config$suffix, ".pdf")
        ),
        plot = results_npx$boxplots,
        width = 8,
        height = 12
    )
}
