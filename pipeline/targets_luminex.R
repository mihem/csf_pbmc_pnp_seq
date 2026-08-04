targets_luminex <- list(
  targets::tar_target(
    luminex_file, "raw/luminex/Luminex_CIDP_GBS.xlsx", format = "file"
  ),
  targets::tar_target(luminex_input, read_luminex_input(luminex_file)),
  targets::tar_target(
    luminex_analysis, analyze_luminex_data(luminex_input, 400L)
  ),
  targets::tar_target(
    luminex_files, write_luminex_artifacts(luminex_analysis), format = "file"
  ),
  targets::tar_target(
    luminex_volcano_results, luminex_volcano_data(luminex_analysis, 500L)
  ),
  targets::tar_target(
    luminex_volcano_file,
    write_luminex_volcano(luminex_volcano_results, 500L),
    format = "file"
  )
)
