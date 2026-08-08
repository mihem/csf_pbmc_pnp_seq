targets_olink <- list(
  targets::tar_target(
    olink_quant_file, "raw/olink/olink_quant_long_filtered.xlsx",
    format = "file"
  ),
  targets::tar_target(
    olink_metadata_file, "lookup/olink_flow_lookup_anonymized.xlsx",
    format = "file"
  ),
  targets::tar_target(
    olink_metadata, read_olink_metadata(olink_metadata_file)
  ),
  targets::tar_target(
    olink_quant_input,
    read_olink_input(olink_quant_file, "Quantified_value")
  ),
  targets::tar_target(
    olink_quant_data,
    prepare_olink_quantified(olink_quant_input, olink_metadata)
  )
)
