targets_projectil <- list(
  tar_target(
    projectil_reference_file,
    file.path("objects", "conventional_cd8_projectil_reference.qs"),
    format = "file"
  ),
  tar_target(
    projectil_readiness,
    assess_projectil_readiness(sc_annotated, projectil_reference_file)
  ),
  tar_target(
    projectil_readiness_file,
    write_projectil_readiness(
      projectil_readiness,
      file.path(projectil_result_dir(), "readiness.csv")
    ),
    format = "file"
  )
)
