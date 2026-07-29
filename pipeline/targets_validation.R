targets_validation <- list(
  tar_target(legacy_sc_list_file, "objects/sc_list.qs", format = "file"),
  tar_target(legacy_sc_merge_pre_file, "objects/sc_merge_pre.qs", format = "file"),
  tar_target(legacy_sc_azimuth_file, "objects/sc_merge_azimuth.qs", format = "file"),
  tar_target(legacy_sc_integrated_file, "objects/sc_merge_batch.qs", format = "file"),
  tar_target(legacy_sc_annotated_file, "objects/sc_merge.qs", format = "file"),
  tar_target(validate_sc_list_result, validate_sc_list(sc_list, legacy_sc_list_file)),
  tar_target(
    validate_sc_merge_pre_result,
    validate_preprocessed_merge(sc_merge_pre, legacy_sc_merge_pre_file)
  ),
  tar_target(
    validate_sc_preprocessed_result,
    validate_preprocessed_reductions(sc_preprocessed, legacy_sc_integrated_file)
  ),
  tar_target(
    validate_sc_azimuth_result,
    validate_azimuth(sc_azimuth, legacy_sc_azimuth_file)
  ),
  tar_target(
    validate_sc_integrated_result,
    validate_integration(sc_integrated, legacy_sc_integrated_file)
  ),
  tar_target(
    validate_sc_annotated_result,
    validate_annotation(sc_annotated, legacy_sc_annotated_file)
  )
)
