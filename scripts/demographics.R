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

lookup$ncv_tibial_motoric


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
    width = 5,
    height = 5
)

# NCV Analysis ----
# Compare tibial nerve motor conduction velocity between groups
ncv_tibial_motoric_plot <- create_boxplot(
    data = lookup,
    x_var = "diagnosis",
    y_var = "ncv_tibial_motoric",
    group_var = "diagnosis",
    title = "motoric NCV tibial nerve (m/s)",
    color_palette = sc_merge@misc$diagnosis_col
)
ggsave(
    file.path("results", "demographics", "boxplot_ncv_tibial_motoric.pdf"),
    plot = ncv_tibial_motoric_plot,
    width = 5,
    height = 5
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
    width = 5,
    height = 5
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
    width = 5,
    height = 5
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
    width = 5,
    height = 5
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
    width = 5,
    height = 5
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
    width = 5,
    height = 5
)

