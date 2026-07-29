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
    olink_overview_table_file,
    write_table_workbook(
      demographics_overview(demographics_olink_data, "orbis_id"),
      file.path(table_result_dir(), "overview_table_olink.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    flow_overview_table_file,
    write_table_workbook(
      demographics_overview(demographics_flow_data, "patient"),
      file.path(table_result_dir(), "overview_table_flow.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    publication_tcr_table_files,
    write_named_table_workbooks(tcr_report_tables),
    format = "file"
  ),
  tar_target(
    publication_sukenikova_table_file,
    write_table_workbook(
      sukenikova_reactive_table,
      file.path(table_result_dir(), "sukenikova_reactive_table.xlsx")
    ),
    format = "file"
  )
)
