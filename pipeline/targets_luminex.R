targets_luminex <- c(list(
  targets::tar_target(
    luminex_file, "raw/luminex/Luminex_CIDP_GBS.xlsx", format = "file"
  ),
  targets::tar_target(
    luminex_control_file, "raw/luminex/Luminex_NIRAE.xlsx", format = "file"
  ),
  targets::tar_target(luminex_disease_input, read_luminex_input(luminex_file)),
  targets::tar_target(
    luminex_control_input, read_luminex_controls(luminex_control_file)
  ),
  targets::tar_target(
    luminex_input,
    combine_luminex_inputs(luminex_disease_input, luminex_control_input)
  ),
  targets::tar_target(
    luminex_analysis, analyze_luminex_data(luminex_input, 400L)
  ),
  targets::tar_target(
    luminex_files, write_luminex_artifacts(luminex_analysis), format = "file"
  )
  ), tarchetypes::tar_map(
    values = luminex_contrasts(), names = comparison, unlist = TRUE,
    targets::tar_target(
      luminex_volcano_results,
      luminex_volcano_data(
        luminex_analysis, comparison, group1, group2, seed
      )
    ),
    targets::tar_target(
      luminex_volcano_file,
      write_luminex_volcano(
        luminex_volcano_results, comparison, group1, group2, seed
      ),
      format = "file"
    )
  ))
