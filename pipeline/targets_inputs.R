targets_inputs <- list(
  tar_target(analysis_config_file, "config/analysis.yml", format = "file"),
  tar_target(analysis_config, read_analysis_config(analysis_config_file)),
  tar_target(trust4_config, read_trust4_config(analysis_config_file)),
  tar_target(
    lookup_file,
    analysis_config$paths$lookup,
    format = "file"
  ),
  tar_target(
    qc_thresholds_file,
    analysis_config$paths$qc_thresholds,
    format = "file"
  ),
  tar_target(rna_manifest, discover_rna_inputs(analysis_config$paths$raw_rna)),
  tar_target(rna_files, rna_input_files(rna_manifest), format = "file"),
  tar_target(manual_files, manual_annotation_files(analysis_config), format = "file"),
  tar_target(lookup, clean_lookup(lookup_file)),
  tar_target(donor_assignments, read_donor_assignments(rna_manifest, rna_files))
)
