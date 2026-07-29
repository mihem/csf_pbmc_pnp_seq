targets_olink <- c(list(
  targets::tar_target(
    olink_quant_file, "raw/olink/olink_quant_long_filtered.xlsx", format = "file"
  ),
  targets::tar_target(
    olink_npx_file, "raw/olink/olink_npx_long_filtered.xlsx", format = "file"
  ),
  targets::tar_target(
    olink_metadata_file, "lookup/olink_flow_lookup_anonymized.xlsx", format = "file"
  ),
  targets::tar_target(olink_metadata, read_olink_metadata(olink_metadata_file)),
  targets::tar_target(
    olink_quant_input, read_olink_input(olink_quant_file, "Quantified_value")
  ),
  targets::tar_target(olink_npx_input, read_olink_input(olink_npx_file, "NPX")),
  targets::tar_target(
    olink_quant_diagnosis,
    analyze_olink_data(
      olink_quant_input, olink_metadata, "Quantified_value", "diagnosis", "pg/ml", 123L
    )
  ),
  targets::tar_target(
    olink_quant_group2,
    analyze_olink_data(
      olink_quant_input, olink_metadata, "Quantified_value", "group2", "pg/ml", 124L
    )
  ),
  targets::tar_target(
    olink_npx_diagnosis,
    analyze_olink_data(olink_npx_input, olink_metadata, "NPX", "diagnosis", "NPX", 125L)
  ),
  targets::tar_target(
    olink_npx_group2,
    analyze_olink_data(olink_npx_input, olink_metadata, "NPX", "group2", "NPX", 126L)
  ),
  targets::tar_target(
    olink_quant_diagnosis_files,
    write_olink_artifacts(olink_quant_diagnosis, "quant"), format = "file"
  ),
  targets::tar_target(
    olink_quant_group2_files,
    write_olink_artifacts(olink_quant_group2, "quant"), format = "file"
  ),
  targets::tar_target(
    olink_npx_diagnosis_files,
    write_olink_artifacts(olink_npx_diagnosis, "npx"), format = "file"
  ),
  targets::tar_target(
    olink_npx_group2_files,
    write_olink_artifacts(olink_npx_group2, "npx"), format = "file"
  )
  ), tarchetypes::tar_map(
    values = olink_contrasts(), names = comparison, unlist = TRUE,
    targets::tar_target(
      olink_quant_volcano_data,
      olink_volcano_data(olink_quant_diagnosis, group1, group2, 200L)
    ),
    targets::tar_target(
      olink_quant_volcano_file,
      write_olink_volcano(olink_quant_volcano_data, comparison, "quant", 200L),
      format = "file"
    ),
    targets::tar_target(
      olink_npx_volcano_data,
      olink_volcano_data(olink_npx_diagnosis, group1, group2, 300L)
    ),
    targets::tar_target(
      olink_npx_volcano_file,
      write_olink_volcano(olink_npx_volcano_data, comparison, "npx", 300L),
      format = "file"
    )
  ))
