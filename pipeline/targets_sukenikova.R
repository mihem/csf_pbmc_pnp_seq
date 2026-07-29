targets_sukenikova <- list(
  tar_target(
    sukenikova_reactive_clone_file,
    file.path("lookup", "sukenikova_nature_gbs_supp_table_2.csv"),
    format = "file"
  ),
  tar_target(
    sukenikova_reactive_clones,
    read_sukenikova_reactive_clones(sukenikova_reactive_clone_file)
  ),
  tar_target(sukenikova_study_clones, extract_study_trb_clones(sc_tcr)),
  tar_target(
    sukenikova_shared_clones,
    match_sukenikova_reactive_clones(
      sukenikova_study_clones, sukenikova_reactive_clones
    )
  ),
  tar_target(
    sukenikova_reactive_table,
    format_sukenikova_reactive_table(sukenikova_shared_clones)
  ),
  tar_target(
    sukenikova_shared_clones_file,
    write_sukenikova_workbook(
      sukenikova_shared_clones,
      file.path(sukenikova_result_dir(), "sukenikova_reactive_shared.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    sukenikova_reactive_table_file,
    write_sukenikova_workbook(
      sukenikova_reactive_table,
      file.path(sukenikova_result_dir(), "sukenikova_reactive_table.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    sukenikova_reactive_table_pdf,
    write_sukenikova_table_pdf(
      sukenikova_reactive_table,
      file.path(sukenikova_result_dir(), "sukenikova_reactive_table.pdf")
    ),
    format = "file"
  )
)
