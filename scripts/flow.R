#####################################
# basic flow analysis
#####################################

# load libraries ------
library(tidyverse)
library(readxl)
library(patchwork)
library(broom)
library(pals)
library(janitor)

# source functions ----
source(file.path("scripts", "flow_helper.R"))

# define order and colors
diagnosis_order <- c(
  "CTRL",
  "CIAP",
  "CIDP",
  "GBS",
  "MAG",
  "MFS",
  "PNC",
  "CAN",
  "PPN"
)

diagnosis_color <- setNames(
  pals::cols25(length(diagnosis_order)),
  diagnosis_order
)

group_order <- c("CTRL", "PNP")
group_color <- setNames(pals::cols25(length(group_order)), group_order)

# load flow and lookup file ----
# load flow data, keep only first measurement if multiple are available
flow_pre <-
  read_excel(file.path("raw", "flow", "flowbasic_v3.xlsx")) |>
  mutate(date = as_date(date)) |>
  group_by(last_name, first_name, tissue) |>
  filter(date == min(date)) |>
  ungroup()

lookup <-
  read_excel(file.path("lookup", "SEED_lookup_v9.xlsx")) |>
  janitor::clean_names() |>
  mutate(age = lubridate::time_length(difftime(date, birth_date), "years")) |>
  mutate(diagnosis = factor(diagnosis, levels = diagnosis_order)) |>
  mutate(group = factor(group, levels = group_order))

# sanity checks
lookup |>
  anti_join(flow_pre, join_by(last_name, first_name)) |>
  select(last_name, first_name, birth_date, date)

flow_pre |>
  anti_join(lookup, join_by(last_name, first_name, date)) |>
  print(n = Inf)

# join flow and lookup
flow <-
  flow_pre |>
  inner_join(lookup, join_by(last_name, first_name)) |>
  (function(df) split(df, df$tissue))()

# Extract flow variables -----
flow_vars <-
  flow$CSF |>
  select(Gran:intMono) |>
  names()

# Create all boxplots -----
# CSF plots
createAndSaveBoxplots(
  data = flow$CSF,
  flow_vars = flow_vars,
  group = "group",
  cols = group_color,
  output_file = "csf_con_plots_group.pdf",
  width = 5,
  height = 15
)

createAndSaveBoxplots(
  data = flow$CSF,
  flow_vars = flow_vars,
  group = "diagnosis",
  cols = diagnosis_color,
  output_file = "csf_con_plots_dx.pdf",
  width = 10,
  height = 15
)

# Blood plots
createAndSaveBoxplots(
  data = flow$blood,
  flow_vars = flow_vars,
  group = "group",
  cols = group_color,
  output_file = "blood_con_plots_group.pdf",
  width = 5,
  height = 15
)

createAndSaveBoxplots(
  data = flow$blood,
  flow_vars = flow_vars,
  group = "diagnosis",
  cols = diagnosis_color,
  output_file = "blood_con_plots_dx.pdf",
  width = 10,
  height = 15
)

# abundance volcano plot ----
# set colors for flow variables
flow_vars_cols <- setNames(pals::cols25(length(flow_vars)), flow_vars)

# Create combinations for different groupings
combinations_group <- createCombinations(
  conditions = c("PNP", "CTRL"),
  group_column_name = "group"
)

combinations_group2 <- createCombinations(
  conditions = c("IN", "NIN", "CTRL"),
  group_column_name = "group2"
)

combinations_diagnosis <- createCombinations(
  conditions = c("CIDP", "GBS", "CIAP", "CTRL"),
  group_column_name = "diagnosis"
)

# Create a configuration table for all volcano plot comparisons
volcano_configs <- bind_rows(
  combinations_group,
  combinations_group2,
  combinations_diagnosis
)

# Generate all volcano plots
volcano_results_list <- pmap(
  volcano_configs,
  function(tissue, condition1, condition2, group_column) {
    createVolcanoPlot(
      data = flow[[tissue]],
      group_column = group_column,
      group1 = condition1,
      group2 = condition2,
      tissue = tissue
    )
  }
)
