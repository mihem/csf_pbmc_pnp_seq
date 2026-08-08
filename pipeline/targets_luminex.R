targets_luminex <- list(
  targets::tar_target(
    luminex_file, "raw/luminex/Luminex_CIDPvsGBS.xlsx", format = "file"
  ),
  targets::tar_target(luminex_input, read_luminex_input(luminex_file))
)
