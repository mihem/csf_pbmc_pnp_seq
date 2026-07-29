library(targets)
library(tarchetypes)

tar_source("R")

source("pipeline/targets_inputs.R")
source("pipeline/targets_preprocess.R")
source("pipeline/targets_mapping.R")
source("pipeline/targets_integration.R")
source("pipeline/targets_annotation.R")
source("pipeline/targets_abundance.R")
source("pipeline/targets_deg.R")
source("pipeline/targets_enrichment.R")
source("pipeline/targets_liana.R")
source("pipeline/targets_tcr.R")
source("pipeline/targets_sukenikova.R")
source("pipeline/targets_flow.R")
source("pipeline/targets_olink.R")
source("pipeline/targets_correlation.R")
source("pipeline/targets_demographics.R")
source("pipeline/targets_projectil.R")
source("pipeline/targets_cerebro.R")
source("pipeline/targets_tables.R")

tar_option_set(
  packages = c(
    "Seurat",
    "Azimuth",
    "broom",
    "cerebroAppLite",
    "dplyr",
    "edgeR",
    "emmeans",
    "factoextra",
    "fgsea",
    "ggplot2",
    "ggrepel",
    "ggsignif",
    "glue",
    "janitor",
    "liana",
    "Libra",
    "limma",
    "lme4",
    "lubridate",
    "msigdbr",
    "patchwork",
    "permFDP",
    "Polychrome",
    "pals",
    "purrr",
    "qs",
    "rcna",
    "readr",
    "readxl",
    "scRepertoire",
    "scMisc",
    "SeuratObject",
    "singlet",
    "speckle",
    "STACAS",
    "stringr",
    "tibble",
    "tidyr",
    "viridis",
    "withr",
    "writexl",
    "yaml"
  ),
  format = "qs",
  error = "abridge",
  memory = "transient",
  garbage_collection = TRUE
)

pipeline <- c(
  targets_inputs,
  targets_preprocess,
  targets_mapping,
  targets_integration,
  targets_annotation,
  targets_abundance,
  targets_deg,
  targets_enrichment,
  targets_liana,
  targets_tcr,
  targets_sukenikova,
  targets_flow,
  targets_olink,
  targets_correlation,
  targets_demographics,
  targets_projectil,
  targets_tables
)

if (include_heavy_targets()) {
  pipeline <- c(pipeline, targets_cerebro)
}

pipeline
