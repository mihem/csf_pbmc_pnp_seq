##################################################
# analyze spectral flow cytometry data with Aurora
##################################################

# Load libraries -----
library(tidyverse)
library(readxl)
library(pals)

# read flow data -----
aurora_flow_pre <-
    read_excel(file.path("raw", "aurora", "final_data.xlsx")) |>
    janitor::clean_names() |>
    dplyr::select(file_id, starts_with("del5")) |>
    rename_with(function(x) gsub("del5_", "", x)) |>
    mutate(file_id = gsub("^\\S+\\s+", "", file_id))

# Read the metadata file
aurora_metadata <- read_xlsx(file.path("lookup", "olink_flow_lookup.xlsx"))

# Define group configurations
group_configs <- list(
    diagnosis = list(
        levels = c("CTRL", "GBS", "CIDP"),
        colors = setNames(pals::cols25(3), c("CTRL", "GBS", "CIDP")),
        suffix = ""
    ),
    group2 = list(
        levels = c("CTRL", "IN"),
        colors = setNames(
            pals::cols25(2),
            c("CTRL", "IN")
        ),
        suffix = "_meta"
    )
)
