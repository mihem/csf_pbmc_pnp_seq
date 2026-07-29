targets_tcr <- list(
  tar_target(tcr_manifest, discover_tcr_inputs(file.path("raw", "tcr"))),
  tar_target(tcr_files, tcr_input_files(tcr_manifest), format = "file"),
  tar_target(
    tcr_contigs,
    import_tcr_contigs(tcr_manifest, tcr_files, donor_assignments, lookup)
  ),
  tar_target(combined_tcr, combine_tcr_contigs(tcr_contigs)),
  tar_target(sc_tcr, annotate_tcr_cells(sc_annotated, combined_tcr)),
  tar_target(
    legacy_combined_tcr_file,
    file.path("objects", "combined_tcr.qs"),
    format = "file"
  ),
  tar_target(
    legacy_sc_tcr_file,
    file.path("objects", "sc_tcr.qs"),
    format = "file"
  ),
  tar_target(
    validate_combined_tcr_result,
    validate_combined_tcr_exact(combined_tcr, legacy_combined_tcr_file)
  ),
  tar_target(
    validate_sc_tcr_result,
    validate_sc_tcr_exact(sc_tcr, legacy_sc_tcr_file)
  ),
  tar_target(tcr_report_tables, make_tcr_report_tables(sc_tcr)),
  tar_target(
    tcr_report_table_files,
    write_tcr_report_tables(tcr_report_tables),
    format = "file"
  ),
  tar_target(
    tcr_basic_plot_files,
    write_tcr_basic_plots(combined_tcr, sc_tcr, seed = 42L),
    format = "file"
  ),
  tar_target(
    tcr_alluvial_plot_files,
    write_tcr_alluvial_plots(sc_tcr, seed = 42L),
    format = "file"
  ),
  tar_target(
    tcr_abundance_table_files,
    write_tcr_abundance_tables(sc_tcr),
    format = "file"
  ),
  tar_target(
    tcr_abundance_plot_files,
    write_tcr_abundance_plots(sc_tcr),
    format = "file"
  ),
  tar_target(
    tcr_invariant_plot_files,
    write_tcr_invariant_plots(sc_tcr, seed = 42L),
    format = "file"
  ),
  tar_target(
    legacy_tcr_report_table_files,
    tcr_legacy_table_paths(),
    format = "file"
  ),
  tar_target(
    validate_tcr_report_tables_result,
    validate_tcr_report_tables(tcr_report_tables, legacy_tcr_report_table_files)
  ),
  tar_target(tcr_comparison_manifest, tcr_comparison_input_manifest()),
  tar_target(tcr_comparison_input_files, tcr_comparison_manifest, format = "file"),
  tar_target(
    tcr_comparison,
    prepare_tcr_comparison(sc_tcr, tcr_comparison_input_files, seed = 42L)
  ),
  tar_target(
    tcr_comparison_plot_files,
    write_tcr_comparison_plots(tcr_comparison, seed = 42L),
    format = "file"
  )
)
