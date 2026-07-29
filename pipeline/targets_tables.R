targets_tables <- list(
  tar_target(
    scrna_overview_table_file,
    write_scrna_overview_table(
      lookup,
      file.path(table_result_dir(), "overview_table.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    publication_table_files,
    collect_publication_tables(c(
      scrna_overview_table_file,
      demographics_olink_files,
      demographics_flow_files,
      tcr_report_table_files,
      tcr_abundance_table_files,
      sukenikova_reactive_table_file,
      abundance_table_sample_file,
      abundance_table_tissue_diagnosis_file,
      abundance_table_tissue_group_file,
      abundance_table_csf_sample_file
    )),
    format = "file"
  )
)
