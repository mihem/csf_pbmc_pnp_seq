###########################################
# basic flow analysis
# requires running flow_pre.R first
###########################################

# load libraries ------
library(tidyverse)
library(readxl)
library(patchwork)
library(broom)
library(pals)
library(janitor)
library(qs)

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

# load preprocessed flow data from this study and  Heming et al . Front Imm Paper ----
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.00515/full
flow <- qread(file.path("objects", "flow_pre.qs")) 

# Extract flow variables -----
flow_vars <-
  flow$CSF |>
  select(Lymph:intMono) |>
  names()

# Create all boxplots -----
# CSF plots
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
volcano_configs <- createCombinations(
  conditions = c("CIDP", "GBS", "CTRL"),
  group_column_name = "diagnosis"
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
      tissue = tissue,
      output_dir = "flow"
    )
  }
)

# Batch comparison plots for all flow variables ----
flow_batch_plots <- lapply(
  flow_vars,
  function(var) {
    plot <-
      flow$CSF |>
      mutate(
        batch = ifelse(
          grepl(patient, pattern = "P\\d{2}"),
          "scRNAseq",
          "Frontiers"
        )
      ) |>
      select(patient, diagnosis, batch, all_of(var)) |>
      dplyr::filter(diagnosis %in% c("CTRL", "CIDP", "GBS")) |>
      ggplot(aes(x = batch, y = .data[[var]])) +
      geom_boxplot() +
      facet_wrap(~diagnosis) +
      labs(y = var, x = "Batch")

    ggsave(
      file.path("results", "flow", "batch", paste0(var, "_csf.pdf")),
      plot = plot,
      width = 5,
      height = 3
    )

    return(plot)
  }
)
