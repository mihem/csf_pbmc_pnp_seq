write_targets_dag <- function(path = "results/pipeline_dag.html") {
  stopifnot(
    length(path) == 1L,
    !is.na(path),
    nzchar(path),
    requireNamespace("htmlwidgets", quietly = TRUE),
    requireNamespace("visNetwork", quietly = TRUE)
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  asset_directory <- paste0(tools::file_path_sans_ext(path), "_files")
  unlink(asset_directory, recursive = TRUE, force = TRUE)
  graph <- targets::tar_visnetwork(
    targets_only = TRUE,
    outdated = TRUE,
    physics = FALSE
  )
  htmlwidgets::saveWidget(
    graph,
    file = path,
    selfcontained = TRUE,
    title = "CSF/PBMC targets pipeline DAG"
  )
  unlink(asset_directory, recursive = TRUE, force = TRUE)
  path
}
