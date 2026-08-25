targets_correlation <- c(
  list(
    targets::tar_target(
      severity_subgroup_plot_file,
      write_severity_subgroup_plot(lookup_severity, 123L),
      format = "file"
    )
  ),
  tarchetypes::tar_map(
    values = correlation_jobs(), names = job, unlist = TRUE,
    targets::tar_target(
      cross_modality_correlation,
      run_abundance_clinical_correlation(
        sc_annotated, lookup_correlation, score, tissue
      )
    ),
    targets::tar_target(
      cross_modality_correlation_files,
      write_correlation_result(cross_modality_correlation),
      format = "file"
    )
  )
)
