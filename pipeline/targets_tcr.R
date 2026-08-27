targets_tcr <- list(
  tar_target(tcr_manifest, discover_tcr_inputs(file.path("raw", "tcr"))),
  tar_target(tcr_files, tcr_input_files(tcr_manifest), format = "file"),
  tar_target(
    tcr_contigs,
    import_tcr_contigs(
      tcr_manifest, tcr_files, donor_assignments, lookup_identity
    )
  ),
  tar_target(combined_tcr, combine_tcr_contigs(tcr_contigs)),
  tar_target(sc_tcr, annotate_tcr_cells(sc_annotated, combined_tcr)),
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
    tcr_shared_clone_plot_files,
    write_tcr_shared_clone_plots(
      sc_tcr, lookup_tcr_shared, tcr_report_tables
    ),
    format = "file"
  ),
  tar_target(
    tcr_albumin_quotient_plot_file,
    write_tcr_albumin_quotient_plot(
      sc_tcr, lookup_albumin, tcr_report_tables
    ),
    format = "file"
  ),
  tar_target(
    tcr_csf_cell_count_plot_file,
    write_tcr_csf_cell_count_plot(
      sc_tcr, lookup_csf_cell_count, tcr_report_tables
    ),
    format = "file"
  ),
  tar_target(
    tcr_barrier_adjustment,
    prepare_tcr_barrier_adjustment(
      tcr_report_tables, lookup_tcr_barrier
    )
  ),
  tar_target(
    tcr_barrier_adjustment_table_file,
    write_tcr_barrier_adjustment_table(tcr_barrier_adjustment),
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
