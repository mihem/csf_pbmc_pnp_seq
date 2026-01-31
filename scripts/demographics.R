#===============================================================================
# Demographics Analysis Script
#===============================================================================
# Purpose: Analyze and visualize patient characteristics across disease groups:
#===============================================================================

# libraries ---
library(tidyverse)
library(qs)
library(pals)
library(readxl)

# source helper functions
source(file.path("scripts", "demographics_helper.R"))

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

lookup <-
  qread(file.path("objects", "lookup.qs")) |>
  dplyr::filter(grepl("scRNA", cohort))

# Define plot configurations ----
boxplot_configs <- list(
  age = list(
    y_var = "age",
    title = "age",
    geom_type = "point",
    width = 3.5,
    height = 3
  ),
  ncv_tibial_motoric = list(
    y_var = "ncv_tibial_motoric",
    title = "motoric NCV tibial nerve (m/s)",
    geom_type = "jitter",
    width = 3,
    height = 3
  ),
  incat = list(
    y_var = "incat_at_lumbar_puncture",
    title = "INCAT score",
    geom_type = "jitter",
    width = 3,
    height = 3
  ),
  incat_progress = list(
    y_var = "incat_progress",
    title = "INCAT score progress",
    geom_type = "jitter",
    width = 5,
    height = 5
  ),
  disease_duration = list(
    y_var = "disease_duration_in_months",
    title = "disease duration (months)",
    geom_type = "jitter",
    width = 3.5,
    height = 3
  ),
  csf_protein = list(
    y_var = "csf_protein",
    title = "CSF protein (mg/L)",
    geom_type = "jitter",
    width = 3.8,
    height = 3
  )
)

barplot_configs <- list(
  disease = list(
    fill_var = "diagnosis",
    title = "diagnosis",
    color_palette = pals::cols25(9),
    width = 4.5,
    height = 3
  ),
  sex = list(
    fill_var = "sex",
    title = "sex",
    color_palette = rev(pals::cols25(2)),
    width = 4.5,
    height = 3
  ),
  therapy_status = list(
    fill_var = "therapy",
    title = "therapy status",
    color_palette = pals::cols25(2),
    width = 4.5,
    height = 3
  )
)

# Generate and save boxplots ----
create_and_save_boxplots(
  data = lookup,
  configs = boxplot_configs,
  x_var = "diagnosis",
  group_var = "diagnosis",
  color_palette = sc_merge@misc$diagnosis_col,
  output_dir = file.path("results", "demographics")
)

# Generate and save barplots ----
create_and_save_barplots(
  data = lookup,
  configs = barplot_configs,
  x_var = "diagnosis",
  output_dir = file.path("results", "demographics")
)

# Olink cohort analysis ----
olink_quant <- read_xlsx(file.path("raw", "olink", "olink_quant_long_filtered.xlsx"))
olink_metadata <- read_xlsx(file.path("lookup", "olink_flow_lookup.xlsx"))

olink_patients <- olink_quant |>
  left_join(olink_metadata, by = "SampleID") |>
  dplyr::select(SampleID, orbis_id, age, sex, diagnosis) |>
  distinct(SampleID, orbis_id, .keep_all = TRUE)

overview_table_olink <- create_overview_table(olink_patients)

writexl::write_xlsx(
  overview_table_olink,
  file.path("results", "table", "overview_table_olink.xlsx")
)
