targets_correlation <- c(
  list(
    targets::tar_target(
      severity_subgroup_plot_file,
      write_severity_subgroup_plot(lookup, 123L),
      format = "file"
    )
  ),
  tarchetypes::tar_map(
    values = correlation_jobs(), names = job, unlist = TRUE,
    targets::tar_target(
      cross_modality_correlation,
      run_abundance_clinical_correlation(sc_annotated, lookup, score, tissue)
    ),
    targets::tar_target(
      cross_modality_correlation_files,
      write_correlation_result(cross_modality_correlation),
      format = "file"
    ),
    targets::tar_target(
      legacy_cross_modality_correlation_file,
      file.path(
        "results", "correlation",
        paste0("correlation_", score, "_", tissue, ".csv")
      ),
      format = "file"
    ),
    targets::tar_target(
      cross_modality_correlation_validation,
      validate_correlation_result(
        cross_modality_correlation_files,
        legacy_cross_modality_correlation_file
      )
    )
  ),
  tarchetypes::tar_map(
    values = tibble::tibble(
      score = correlation_scores(),
      rcna_job = paste0("rcna_", correlation_scores()),
      seed = 1000L + seq_along(correlation_scores())
    ),
    names = rcna_job,
    unlist = TRUE,
    targets::tar_target(
      rcna_association,
      run_rcna_association(sc_annotated, lookup, score, seed)
    ),
    targets::tar_target(
      rcna_plot_file,
      write_rcna_plot(rcna_association, score),
      format = "file"
    )
  )
)
