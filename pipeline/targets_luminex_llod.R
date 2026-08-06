targets_luminex_llod <- c(list(
  targets::tar_target(
    luminex_llod_analysis, analyze_luminex_llod(luminex_input)
  ),
  targets::tar_target(
    luminex_llod_files,
    write_luminex_llod_artifacts(luminex_llod_analysis),
    format = "file"
  )
  ), tarchetypes::tar_map(
    values = luminex_contrasts(), names = comparison, unlist = TRUE,
    targets::tar_target(
      luminex_llod_volcano_file,
      write_luminex_llod_volcano(
        luminex_llod_analysis, comparison, group1, group2, seed
      ),
      format = "file"
    )
  ))
