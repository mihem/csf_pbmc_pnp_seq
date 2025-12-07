##################################################
# analyze spectral flow cytometry data with Aurora
##################################################

# Load libraries -----
library(tidyverse)
library(readxl)
library(pals)

# source functions -----
source(file.path("scripts", "flow_helper.R"))

# read flow data -----
aurora_flow_pre <-
    read_excel(file.path("raw", "aurora", "final_data.xlsx")) |>
    janitor::clean_names() |>
    dplyr::select(file_id, starts_with("del5")) |>
    rename_with(function(x) gsub("del5_", "", x)) |>
    mutate(file_id = gsub("^\\S+\\s+", "", file_id))

# Read the metadata file
aurora_metadata <- read_xlsx(file.path("lookup", "olink_flow_lookup.xlsx")) |>
    rename(file_id = SampleID)

# Join flow data with metadata
aurora_flow <- aurora_flow_pre |>
    inner_join(aurora_metadata, by = join_by(file_id))

# sanity checks
identical(nrow(aurora_flow_pre), nrow(aurora_flow))
aurora_flow_pre |>
    anti_join(aurora_metadata, by = join_by(file_id))

# Extract flow variables -----
flow_vars <-
    aurora_flow |>
    select(cd4_percent_t_cells:n_kmem_tn_fa) |>
    names()

# abundance volcano plot ----
# set colors for flow variables (use colorRampPalette for many variables)
flow_vars_cols <- setNames(
    colorRampPalette(pals::cols25())(length(flow_vars)),
    flow_vars
)

flow_vars_cols <- setNames(
    rep("black", length(flow_vars)),
    flow_vars
)

# Create combinations for different groupings
combinations_group <- createCombinations(
    conditions = c("PNP", "CTRL"),
    group_column_name = "group",
    tissue = "CSF"
)

combinations_diagnosis <- createCombinations(
    conditions = c("CIDP", "GBS", "CTRL"),
    group_column_name = "diagnosis",
    tissue = "CSF"
)

# Create a configuration table for all volcano plot comparisons
volcano_configs <- bind_rows(
    combinations_group,
    combinations_diagnosis
)

# Generate all volcano plots
volcano_results_list <- pmap(
    volcano_configs,
    function(tissue, condition1, condition2, group_column) {
        createVolcanoPlot(
            data = aurora_flow,
            group_column = group_column,
            group1 = condition1,
            group2 = condition2,
            tissue = tissue,
            output_dir = "aurora",
            width = 10,
            height = 10,
            top_n = 100
        )
    }
)
