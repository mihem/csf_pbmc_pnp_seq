##################################################
# analyze spectral flow cytometry data with Aurora
##################################################

# Load libraries -----
library(tidyverse)
library(readxl)
library(pals)
library(permFDP)

# source functions -----
source(file.path("scripts", "flow_helper.R"))

names(aurora_flow_pre)
flow_vars[637]

# read flow data -----
aurora_flow_pre <-
  read_excel(file.path("raw", "aurora", "final_data.xlsx")) |>
  rename_with(tolower) |>
  janitor::clean_names() |>
  dplyr::select(file_id, starts_with("del5")) |>
  rename_with(function(x) gsub("del5_", "", x, fixed = TRUE)) |>
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
    select(cd4_percent_t_cells:n_kmem_tnfa) |>
#   select(cd4_percent_t_cells:tfh_tnfa) |> # tcell panel only
  names()

# olink data read -----
olink_quant_file <- file.path(
  "raw",
  "olink",
  "olink_quant_long_filtered.xlsx"
)

olink_quant <- read_xlsx(olink_quant_file) |>
  mutate(Quantified_value = as.numeric(Quantified_value))

olink_vars <- tolower(unique(olink_quant$Assay))
olink_vars[olink_vars == "gzma"] <- "gr_a"

# define CD8 TEM markers to filter flow variables -----
cd8tem_3 <- c(
  "CD194",
  "CXCR6",
  "SLAMF6",
  "FASLG",
  "M-CSF",
  "IL-10RA",
  "SIRT2",
  "TNFRSF9"
)

# filter flow variables to those matching olink variables -----
flow_vars <- flow_vars[grepl(
  paste(olink_vars, collapse = "|"),
  # paste(cd8tem_3, collapse = "|"),
  flow_vars,
  ignore.case = TRUE
)]

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
  tissue = "PBMC"
)

combinations_diagnosis <- createCombinations(
  conditions = c("CIDP", "GBS", "CTRL"),
  group_column_name = "diagnosis",
  tissue = "PBMC"
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
      output_dir = file.path("aurora", "olink"),
      width = 10,
      height = 10,
      top_n = NULL
    )
  }
)

# Export volcano plot data to Excel (separate sheets per comparison)
volcano_data_export <- lapply(seq_along(volcano_results_list), function(i) {
  volcano_results_list[[i]]$data |>
    arrange(p.value) |>
    select(
      var,
      log2_ratio,
      p.value,
      p.adj.threshold,
      neg_log10_p,
      significant
    )
})

names(volcano_data_export) <- paste0(
  volcano_configs$condition1,
  "_vs_",
  volcano_configs$condition2,
  "_",
  volcano_configs$tissue
)

writexl::write_xlsx(
  volcano_data_export,
  file.path("results", "aurora", "olink", "volcano_results.xlsx")
)
