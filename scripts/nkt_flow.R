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

# load flow and lookup file ----
# load flow data, keep only first measurement if multiple are available
nkt_flow_pre <-
  read_excel(file.path("raw", "flow_nkt", "nkt_result_v1.xlsx")) |>
  mutate(sample = gsub("(\\d+-V\\d).*", "\\1", Sample))

flow_metadata <-
  read_excel(file.path("lookup", "olink_flow_lookup.xlsx")) |>
  mutate(tissue = ifelse(grepl("C", SampleID), "PBMC", "CSF")) |>
  dplyr::filter(tissue == "PBMC") |>
  mutate(sample = gsub("(\\d+)_(\\d).*", "\\1-V\\2", SampleID))

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
  nkt_flow$CSF |>
  select(Gran:intMono) |>
  names()
