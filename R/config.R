include_heavy_targets <- function() {
  value <- Sys.getenv("TAR_INCLUDE_HEAVY", unset = "false")
  tolower(value) %in% c("1", "true", "yes")
}

read_analysis_config <- function(path) {
  config <- yaml::read_yaml(path)
  config$trust4 <- NULL
  config$paths[c(
    "sural_ic_metadata_umap", "sural_metadata_umap", "sural_trust4"
  )] <- NULL
  stopifnot(
    is.list(config),
    is.list(config$preprocess),
    is.list(config$integration),
    is.list(config$annotation),
    is.list(config$paths)
  )
  config
}

read_trust4_config <- function(path) {
  config <- yaml::read_yaml(path)
  stopifnot(
    is.list(config$trust4),
    is.list(config$trust4$sample_map),
    is.list(config$trust4$patient_map),
    all(c(
      "sural_ic_metadata_umap", "sural_metadata_umap", "sural_trust4"
    ) %in% names(config$paths))
  )
  list(
    paths = config$paths[c(
      "sural_ic_metadata_umap", "sural_metadata_umap", "sural_trust4"
    )],
    trust4 = config$trust4
  )
}

configure_runtime <- function(config) {
  workers <- as.integer(config$project$workers)
  future::plan("multicore", workers = workers)
  options(
    future.globals.maxSize =
      as.numeric(config$project$future_globals_gib) * 1024^3
  )
  invisible(workers)
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
