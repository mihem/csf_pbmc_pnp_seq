targets_abundance <- list(
  tar_target(
    abundance_object,
    prepare_abundance_object(sc_annotated, lookup)
  ),
  tar_target(abundance_tables, make_abundance_tables(abundance_object)),
  tar_target(
    abundance_table_sample_file,
    write_abundance_table(
      abundance_tables$sample,
      file.path(abundance_result_dir(), "abundance_tbl_sc_annotated_sample.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    abundance_table_tissue_diagnosis_file,
    write_abundance_table(
      abundance_tables$tissue_diagnosis,
      file.path(
        abundance_result_dir(),
        "abundance_tbl_sc_annotated_tissue_diagnosis.xlsx"
      )
    ),
    format = "file"
  ),
  tar_target(
    abundance_table_tissue_group_file,
    write_abundance_table(
      abundance_tables$tissue_group,
      file.path(
        abundance_result_dir(), "abundance_tbl_sc_annotated_tissue_group.xlsx"
      )
    ),
    format = "file"
  ),
  tar_target(
    abundance_stacked_csf_sample_file,
    write_stacked_abundance_plot(
      abundance_object, "CSF", "sample",
      file.path(abundance_result_dir(), "stacked_barplot_csf_sample.pdf"), 7
    ),
    format = "file"
  ),
  tar_target(
    abundance_stacked_pbmc_sample_file,
    write_stacked_abundance_plot(
      abundance_object, "PBMC", "sample",
      file.path(abundance_result_dir(), "stacked_barplot_pbmc_sample.pdf"), 7
    ),
    format = "file"
  ),
  tar_target(
    abundance_stacked_csf_diagnosis_file,
    write_stacked_abundance_plot(
      abundance_object, "CSF", "diagnosis",
      file.path(abundance_result_dir(), "stacked_barplot_csf_diagnosis.pdf"), 5
    ),
    format = "file"
  ),
  tar_target(
    abundance_stacked_pbmc_diagnosis_file,
    write_stacked_abundance_plot(
      abundance_object, "PBMC", "diagnosis",
      file.path(abundance_result_dir(), "stacked_barplot_pbmc_diagnosis.pdf"), 5
    ),
    format = "file"
  ),
  tar_target(
    abundance_boxplot_csf_diagnosis_file,
    write_abundance_boxplot(
      abundance_object,
      "CSF",
      file.path(abundance_result_dir(), "boxplot_cluster_csf_diagnosis.pdf")
    ),
    format = "file"
  ),
  tar_target(
    abundance_boxplot_pbmc_diagnosis_file,
    write_abundance_boxplot(
      abundance_object,
      "PBMC",
      file.path(abundance_result_dir(), "boxplot_cluster_pbmc_diagnosis.pdf")
    ),
    format = "file"
  ),
  tar_target(abundance_comparison_config, abundance_comparisons()),
  tar_target(
    abundance_propeller_results,
    run_abundance_propeller(
      abundance_object,
      lookup,
      abundance_comparison_config,
      seed = 123L,
      n_perms = 1000L
    )
  ),
  tar_target(
    abundance_propeller_results_file,
    write_propeller_results(
      abundance_propeller_results,
      file.path(abundance_result_dir(), "propeller_results_all.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    abundance_propeller_plot_files,
    write_propeller_plots(
      abundance_propeller_results,
      abundance_object@misc$cluster_col,
      file.path(abundance_result_dir(), "propeller_plots")
    ),
    format = "file"
  ),
  tar_target(abundance_pca_results, make_abundance_pcas(abundance_object)),
  tar_target(
    abundance_pca_tables_file,
    write_abundance_pca_tables(
      abundance_pca_results,
      file.path(abundance_result_dir(), "abundance_pca.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    abundance_pca_plot_files,
    write_abundance_pca_plots(
      abundance_pca_results,
      abundance_object@misc$diagnosis_col,
      abundance_result_dir(),
      seed = 123L
    ),
    format = "file"
  ),
  tar_target(
    abundance_legacy_table_files,
    legacy_abundance_table_paths(),
    format = "file"
  ),
  tar_target(
    validate_abundance_tables_result,
    validate_abundance_tables(abundance_tables, abundance_legacy_table_files)
  ),
  tar_target(
    abundance_legacy_propeller_file,
    legacy_propeller_path(),
    format = "file"
  ),
  tar_target(
    validate_abundance_propeller_result,
    validate_propeller_results(
      abundance_propeller_results, abundance_legacy_propeller_file
    )
  )
)
