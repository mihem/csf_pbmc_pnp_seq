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
# Compare disability scores between disease groups
# Scale fixed from 0-6 for consistent visualization
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

# Sex Distribution ----
# Visualize male/female distribution across disease groups

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

therapy_plot <- create_barplot(
    data = lookup,
    x_var = "diagnosis",
    fill_var = "therapy_naive",
    title = "therapy status",
    color_palette = rev(pals::cols25(2))
)

ggsave(
    file.path("results", "demographics", "barplot_therapy_status.pdf"),
    plot = therapy_plot,
    width = 5,
    height = 5
)
