#===============================================================================
# Demographics Analysis Script
#===============================================================================
# Purpose: Analyze and visualize patient characteristics across disease groups:
#===============================================================================

# libraries ---
library(tidyverse)
library(qs)
library(pals)

# source helper functions
source(file.path("scripts", "demographics_helper.R"))

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

lookup <-
  qread(file.path("objects", "lookup.qs")) |>
  dplyr::filter(grepl("scRNA", cohort))

# Number of patients per group ----
disease_plot <- create_barplot(
  data = lookup,
  x_var = "diagnosis",
  fill_var = "diagnosis",
  title = "diagnosis",
  color_palette = pals::cols25(9)
)

ggsave(
  file.path("results", "demographics", "barplot_disease.pdf"),
  plot = disease_plot,
  width = 4.5,
  height = 3
)


# Age Analysis ----
# Compare age distribution between disease groups using boxplots
age_plot <- create_boxplot(
  data = lookup,
  x_var = "diagnosis",
  y_var = "age",
  group_var = "diagnosis",
  title = "age",
  color_palette = sc_merge@misc$diagnosis_col
)

ggsave(
  file.path("results", "demographics", "boxplot_age.pdf"),
  plot = age_plot,
  width = 3.5,
  height = 3
)

# NCV Analysis ----
# Compare tibial nerve motor conduction velocity between groups
ncv_tibial_motoric_plot <- create_boxplot(
  data = lookup,
  x_var = "diagnosis",
  y_var = "ncv_tibial_motoric",
  group_var = "diagnosis",
  title = "motoric NCV tibial nerve (m/s)",
  color_palette = sc_merge@misc$diagnosis_col,
  geom_type = "jitter"
)

ggsave(
  file.path("results", "demographics", "boxplot_ncv_tibial_motoric.pdf"),
  plot = ncv_tibial_motoric_plot,
  width = 3,
  height = 3
)

# INCAT Score Analysis ----
incat_plot <- create_boxplot(
  data = lookup,
  x_var = "diagnosis",
  y_var = "incat_at_lumbar_puncture",
  group_var = "diagnosis",
  title = "INCAT score",
  color_palette = sc_merge@misc$diagnosis_col,
  geom_type = "jitter"
)

ggsave(
  file.path("results", "demographics", "boxplot_incat.pdf"),
  plot = incat_plot,
  width = 3,
  height = 3
)

# INCAT Progress Analysis ----
incat_followup_plot <- create_boxplot(
  data = lookup,
  x_var = "diagnosis",
  y_var = "incat_progress",
  group_var = "diagnosis",
  title = "INCAT score progress",
  color_palette = sc_merge@misc$diagnosis_col,
  geom_type = "jitter"
)

ggsave(
  file.path("results", "demographics", "boxplot_incat_progress.pdf"),
  plot = incat_followup_plot,
  width = 5,
  height = 5
)

# Disease Duration Analysis ----
disease_duration_plot <- create_boxplot(
  data = lookup,
  x_var = "diagnosis",
  y_var = "disease_duration_in_months",
  group_var = "diagnosis",
  title = "disease duration (months)",
  color_palette = sc_merge@misc$diagnosis_col,
  geom_type = "jitter"
)

ggsave(
  file.path("results", "demographics", "boxplot_disease_duration.pdf"),
  plot = disease_duration_plot,
  width = 3.5,
  height = 3
)

# CSF Protein Analysis ----
csf_protein_plot <- create_boxplot(
  data = lookup,
  x_var = "diagnosis",
  y_var = "csf_protein",
  group_var = "diagnosis",
  title = "CSF protein (mg/L)",
  color_palette = sc_merge@misc$diagnosis_col,
  geom_type = "jitter"
)

ggsave(
  file.path("results", "demographics", "boxplot_csf_protein.pdf"),
  plot = csf_protein_plot,
  width = 3.8,
  height = 3
)

# Sex Distribution ----
sex_plot <- create_barplot(
  data = lookup,
  x_var = "diagnosis",
  fill_var = "sex",
  color_palette = rev(pals::cols25(2)),
  title = "sex"
)

ggsave(
  file.path("results", "demographics", "barplot_sex.pdf"),
  plot = sex_plot,
  width = 4.5,
  height = 3
)

# Therapy Status Distribution ----
therapy_plot <- create_barplot(
  data = lookup,
  x_var = "diagnosis",
  fill_var = "therapy",
  title = "therapy status",
  color_palette = pals::cols25(2)
)

ggsave(
  file.path("results", "demographics", "barplot_therapy_status.pdf"),
  plot = therapy_plot,
  width = 4.5,
  height = 3
)

# Olink cohort ----
olink_quant_file <- file.path(
  "raw",
  "olink",
  "olink_quant_long_filtered.xlsx"
)

olink_quant <- read_xlsx(olink_quant_file)

olink_metadata_file <- file.path("lookup", "olink_flow_lookup.xlsx")
olink_metadata <- read_xlsx(olink_metadata_file)

olink_patients <- olink_quant |>
  left_join(olink_metadata, by = "SampleID") |>
  dplyr::select(SampleID, orbis_id, age, sex, diagnosis) |>
  distinct(SampleID, orbis_id, .keep_all = TRUE)

overview_table_olink <-
  olink_patients |>
  dplyr::mutate(sex_cat = if_else(sex == "male", 1, 0)) |>
  dplyr::group_by(diagnosis) |>
  dplyr::summarize(
    samples = n(),
    patients = n_distinct(orbis_id),
    age = mean(age, na.rm = TRUE),
    female = (1 - mean(sex_cat, na.rm = TRUE)) * 100
  )

writexl::write_xlsx(
  overview_table_olink,
  file.path("results", "table", "overview_table_olink.xlsx")
)
