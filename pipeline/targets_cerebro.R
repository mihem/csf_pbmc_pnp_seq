targets_cerebro <- list(
  tar_target(cerebro_prepared_object, prepare_cerebro_export(sc_annotated)),
  tar_target(cerebro_marked_object, add_cerebro_marker_genes(cerebro_prepared_object)),
  tar_target(
    cerebro_object_file,
    write_cerebro_object(
      cerebro_marked_object,
      file.path(cerebro_result_dir(), "sc_annotated_cerebro.qs")
    ),
    format = "file"
  ),
  tar_target(
    cerebro_export_file,
    export_cerebro_artifact(
      cerebro_marked_object,
      file.path(cerebro_result_dir(), "sc_annotated_cerebro.crb")
    ),
    format = "file"
  ),
  tar_target(
    cerebro_export_status_file,
    write_cerebro_export_status(
      cerebro_object_file,
      cerebro_export_file,
      file.path(cerebro_result_dir(), "export_status.csv")
    ),
    format = "file"
  )
)
