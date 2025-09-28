##################################################
# Analyze Olink data
##################################################

# Load necessary libraries
library(tidyverse)
library(qs)
library(readxl)
library(tidyxl)

# Read the metadata file
olink_metadata_file <- file.path("lookup", "olink_lookup.xlsx")
olink_metadata <- read_xlsx(olink_metadata_file)

# Read the Olink quantified data
olink_quant_file <- file.path(
    "raw",
    "olink",
    "P049E085_Müller-Miny_QUANT_long.xlsx"
)

olink_quant <- read_xlsx(olink_quant_file) |>
    mutate(Quantified_value = as.numeric(Quantified_value))

# Retrieve formatting information from the Excel file
cells <- tidyxl::xlsx_cells(olink_quant_file)
formats <- tidyxl::xlsx_formats(olink_quant_file)

# Find the column number for 'Quantified_value'
quantified_col <- cells |>
    filter(row == 1, character == "Quantified_value") |>
    pull(col)

# Filter cells to only include the 'Quantified_value' column (excluding header)
quantified_cells <- cells |>
    filter(col == quantified_col, row > 1)

# Identify red-marked cells in the 'Quantified_value' column
red_cells <- quantified_cells |>
    filter(!is.na(local_format_id)) |>
    mutate(
        is_red = map_lgl(
            local_format_id,
            function(x) {
                # Use local_format_id as index into the format arrays
                if (
                    x > 0 &&
                        x <= length(formats$local$fill$patternFill$fgColor$rgb)
                ) {
                    # Check foreground color for red highlighting
                    fg_rgb <- formats$local$fill$patternFill$fgColor$rgb[x]
                    if (!is.na(fg_rgb)) {
                        # Check for light red/pink colors (FFF5AA9C is the red highlighting)
                        return(
                            fg_rgb %in% c("FFF5AA9C", "FFFF0000", "FFFF6B6B")
                        )
                    }
                }
                return(FALSE)
            }
        )
    ) |>
    select(row, col, is_red)

# Join red cell information with the quantified data
quantified_data <- quantified_cells |>
    filter(!is_blank) |>
    left_join(red_cells, by = c("row", "col")) |>
    select(row, content, is_red) |>
    rename(Quantified_value = content, uncertain = is_red) |>
    mutate(
        row_index = row - 1, # Convert to 0-indexed to match data frame rows
        uncertain = ifelse(is.na(uncertain), FALSE, uncertain)
    ) |>
    select(row_index, uncertain)

# Combine with the original olink_quant data
olink_combined <- olink_quant |>
    mutate(row_index = row_number()) |> # Add row index to match with Excel rows
    left_join(
        quantified_data,
        by = "row_index"
    ) |>
    select(-row_index) |> # Remove the helper column
    dplyr::filter(uncertain == FALSE | is.na(uncertain)) # Keep only certain values

olink_combined_metadata <-
    olink_combined |>
    left_join(olink_metadata, by = "SampleID") |>
    relocate(Assay, Quantified_value, group, age, sex)

# Save the cleaned data
qs::qsave(
    olink_combined_metadata,
    file.path("objects", "olink_combined_metadata.qs")
)
