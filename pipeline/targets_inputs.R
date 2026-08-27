targets_inputs <- list(
  tar_target(analysis_config_file, "config/analysis.yml", format = "file"),
  tar_target(analysis_config, read_analysis_config(analysis_config_file)),
  tar_target(trust4_config, read_trust4_config(analysis_config_file)),
  tar_target(config_isolation_check, validate_config_isolation(analysis_config)),
  tar_target(rna_creation_config, rna_creation_settings(analysis_config)),
  tar_target(rna_filter_config, rna_filter_settings(analysis_config)),
  tar_target(
    rna_normalization_config, rna_normalization_settings(analysis_config)
  ),
  tar_target(azimuth_config, azimuth_settings(analysis_config)),
  tar_target(clustering_config, clustering_settings(analysis_config)),
  tar_target(annotation_config, annotation_settings(analysis_config)),
  tar_target(
    trust4_discovery_config, trust4_discovery_settings(trust4_config)
  ),
  tar_target(trust4_main_config, trust4_main_settings(trust4_config)),
  tar_target(trust4_ic_config, trust4_ic_settings(trust4_config)),
  tar_target(trust4_patient_config, trust4_patient_settings(trust4_config)),
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
  tar_target(raw_rna_dir, analysis_config$paths$raw_rna),
  tar_target(rna_manifest, discover_rna_inputs(raw_rna_dir)),
  tar_target(rna_files, rna_input_files(rna_manifest), format = "file"),
  tar_target(
    annotation_file, analysis_config$paths$annotations, format = "file"
  ),
  tar_target(
    annotation_marker_file, analysis_config$paths$markers, format = "file"
  ),
  tar_target(
    nk_overrides_file, analysis_config$paths$nk_overrides, format = "file"
  ),
  tar_target(
    annotation_outlier_files,
    file.path(
      "lookup", paste0(c("cl9", "cl12", "cl23", "cl24"), "_outliers.csv")
    ),
    format = "file"
  ),
  tar_target(lookup, clean_lookup(lookup_file)),
  tar_target(
    lookup_identity,
    project_lookup(lookup, c("patient", "pseudonym"))
  ),
  tar_target(
    lookup_preprocess,
    project_lookup(
      lookup,
      c(
        "patient", "sex", "group", "diagnosis",
        "incat_at_lumbar_puncture", "incat_follow_up",
        "onls_at_lumbar_puncture", "onls_follow_up",
        "mrc_sum_score_60_at_lumbar_puncture",
        "mrc_sum_score_60_follow_up", "icu", "age"
      )
    )
  ),
  tar_target(
    lookup_batch,
    project_lookup(lookup, c("patient", "batch"))
  ),
  tar_target(
    lookup_model,
    project_lookup(
      lookup,
      c("patient", "sex", "age", "group", "group2", "diagnosis")
    )
  ),
  tar_target(
    lookup_treatment,
    project_lookup(
      lookup,
      c(
        "patient", "sex", "age", "group", "group2", "diagnosis",
        "therapy", "immunomodulators"
      )
    )
  ),
  tar_target(
    lookup_correlation,
    project_lookup(
      lookup,
      unique(c("patient", "sex", "age", correlation_scores()))
    )
  ),
  tar_target(
    lookup_severity,
    project_lookup(
      lookup, c("patient", "diagnosis", "incat_at_lumbar_puncture")
    )
  ),
  tar_target(
    lookup_tcr_shared,
    project_lookup(
      lookup,
      c("patient", "disease_duration_in_months", "csf_protein")
    )
  ),
  tar_target(
    lookup_albumin,
    project_lookup(
      lookup, c("patient", "albumin_quotient")
    )
  ),
  tar_target(
    lookup_csf_cell_count,
    project_lookup(lookup, c("patient", "cell_count"))
  ),
  tar_target(
    lookup_demographics,
    project_lookup(
      lookup,
      c(
        "patient", "diagnosis", "sex", "age", "ncv_tibial_motoric",
        "incat_at_lumbar_puncture", "incat_progress",
        "disease_duration_in_months", "csf_protein", "therapy", "bbbd"
      )
    )
  ),
  tar_target(donor_assignments, read_donor_assignments(rna_manifest, rna_files))
)
