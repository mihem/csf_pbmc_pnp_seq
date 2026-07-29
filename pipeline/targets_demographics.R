targets_demographics <- list(
  targets::tar_target(
    demographics_olink_quant_file,
    "raw/olink/olink_quant_long_filtered.xlsx",
    format = "file"
  ),
  targets::tar_target(
    demographics_olink_metadata_file,
    "lookup/olink_flow_lookup_anonymized.xlsx",
    format = "file"
  ),
  targets::tar_target(
    demographics_olink_data,
    prepare_olink_demographics(
      demographics_olink_quant_file, demographics_olink_metadata_file
    )
  ),
  targets::tar_target(
    demographics_flow_data,
    prepare_flow_demographics(flow_data)
  ),
  targets::tar_target(
    demographics_scrna_files,
    write_scrna_demographics(lookup, sc_annotated@misc$diagnosis_col, 123L),
    format = "file"
  ),
  targets::tar_target(
    demographics_olink_files,
    write_modality_demographics(
      demographics_olink_data, "olink", "orbis_id",
      sc_annotated@misc$diagnosis_col, 223L
    ),
    format = "file"
  ),
  targets::tar_target(
    demographics_flow_files,
    write_modality_demographics(
      demographics_flow_data, "flow", "patient",
      sc_annotated@misc$diagnosis_col, 323L
    ),
    format = "file"
  ),
  targets::tar_target(
    demographics_pdf_check,
    check_pdf_artifacts(c(
      demographics_scrna_files,
      demographics_olink_files[grepl("[.]pdf$", demographics_olink_files)],
      demographics_flow_files[grepl("[.]pdf$", demographics_flow_files)]
    ))
  )
)
