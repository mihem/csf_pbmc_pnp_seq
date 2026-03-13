##################################################
# NKT flow analysis
##################################################

# load libraries ------
library(tidyverse)
library(readxl)
library(patchwork)
library(pals)

# source functions ----
source(file.path("scripts", "flow_helper.R"))

# define and colors
group2_colors <- setNames(pals::cols25(2), c("CTRL", "IN"))
diagnosis_colors <- setNames(pals::cols25(3), c("CTRL", "GBS", "CIDP"))

# load flow and lookup file ----
# load flow data, keep only first measurement if multiple are available
nkt_flow_pre <-
  read_excel(file.path("raw", "flow_nkt", "nkt_result_v1.xlsx")) |>
  mutate(sample = gsub("(\\d+-V\\d).*", "\\1", Sample))

flow_metadata <-
  read_excel(file.path("lookup", "olink_flow_lookup.xlsx")) |>
  mutate(tissue = ifelse(grepl("C", SampleID), "PBMC", "CSF")) |>
  dplyr::filter(tissue == "PBMC") |>
  mutate(sample = gsub("(\\d+)_(\\d).*", "\\1-V\\2", SampleID)) |>
  mutate(group2 = factor(group2, levels = c("CTRL", "IN"))) |>
  mutate(diagnosis = factor(diagnosis, levels = c("CTRL", "GBS", "CIDP")))

# sanity checks
anyDuplicated(nkt_flow_pre$sample)

nkt_flow_pre |>
  anti_join(flow_metadata, join_by(sample))

anyDuplicated(flow_metadata$sample)

flow_metadata |>
  group_by(sample) |>
  filter(n() > 1) |>
  select(orbis_id, SampleID)

# join flow and lookup
nkt_flow <-
  nkt_flow_pre |>
  inner_join(flow_metadata, join_by(sample))

# Sanity check 2
nrow(nkt_flow) == nrow(nkt_flow_pre)

# Extract flow variables -----
flow_vars <-
  nkt_flow |>
  select(NKT:NK) |>
  names()

# Create all boxplots -----
createAndSaveBoxplots(
  data = nkt_flow,
  flow_vars = flow_vars,
  group = "group2",
  cols = group2_colors,
  output_file = "nkt_plots_group2.pdf",
  width = 10,
  height = 10
)

createAndSaveBoxplots(
  data = nkt_flow,
  flow_vars = flow_vars,
  group = "diagnosis",
  cols = diagnosis_colors,
  output_file = "nkt_plots_diagnosis.pdf",
  width = 10,
  height = 10
)

# abundance volcano plot ----
# set colors for flow variables
flow_vars_cols <- setNames(pals::cols25(length(flow_vars)), flow_vars)

# Create combinations for different groupings
combinations_group2 <- createCombinations(
  conditions = c("IN", "CTRL"),
  group_column_name = "group2",
  tissues = "PBMC"
)

combinations_diagnosis <- createCombinations(
  conditions = c("CIDP", "GBS", "CTRL"),
  group_column_name = "diagnosis",
  tissues = "PBMC"
)

# Create a configuration table for all volcano plot comparisons
volcano_configs <- bind_rows(
  combinations_group2,
  combinations_diagnosis
)

# Generate all volcano plots
volcano_results_list <- lapply(seq_len(nrow(volcano_configs)), function(i) {
  createVolcanoPlot(
    data = nkt_flow,
    group_column = volcano_configs$group_column[i],
    group1 = volcano_configs$condition1[i],
    group2 = volcano_configs$condition2[i],
    tissue = volcano_configs$tissue[i],
    output_dir = file.path("flow", "nkt")
  )
})
