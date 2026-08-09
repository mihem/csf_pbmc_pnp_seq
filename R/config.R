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
    is.list(config$annotation),
    is.list(config$proteomics),
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

rna_creation_settings <- function(config) {
  config$preprocess[c("min_cells", "min_features")]
}

rna_filter_settings <- function(config) {
  config$preprocess["max_features_default"]
}

rna_normalization_settings <- function(config) {
  config$preprocess[c(
    "normalization_method", "scale_factor", "variable_features",
    "nmf_factors", "nmf_seed"
  )]
}

azimuth_settings <- function(config) {
  config$azimuth[c("reference", "score_threshold")]
}

clustering_settings <- function(config) {
  config$annotation[c("dimensions", "resolutions", "cluster_seed")]
}

annotation_settings <- function(config) {
  resolutions <- as.numeric(unlist(config$annotation$resolutions))
  final_resolution <- as.numeric(config$annotation$final_resolution)
  stopifnot(
    length(final_resolution) == 1L, is.finite(final_resolution),
    any(abs(resolutions - final_resolution) < .Machine$double.eps^0.5)
  )
  list(
    cluster_column = paste0(
      "RNA_snn_res.",
      format(final_resolution, scientific = FALSE, trim = TRUE)
    )
  )
}

trust4_discovery_settings <- function(config) {
  list(sample_map = config$trust4$sample_map)
}

trust4_main_settings <- function(config) {
  config$trust4[c("reduction", "cluster_column")]
}

trust4_ic_settings <- function(config) {
  config$trust4[c("ic_reduction", "ic_cluster_column")]
}

trust4_patient_settings <- function(config) {
  list(patient_map = config$trust4$patient_map)
}

scientific_config_projections <- function(config) {
  list(
    rna_creation = rna_creation_settings(config),
    rna_filter = rna_filter_settings(config),
    rna_normalization = rna_normalization_settings(config),
    azimuth = azimuth_settings(config),
    clustering = clustering_settings(config),
    annotation = annotation_settings(config),
    proteomics = proteomics_settings(config$proteomics$fdr)
  )
}

validate_config_isolation <- function(config) {
  baseline <- scientific_config_projections(config)
  checks <- list(
    rna_creation = function(x) { x$preprocess$min_cells <- x$preprocess$min_cells + 1; x },
    rna_filter = function(x) {
      x$preprocess$max_features_default <- x$preprocess$max_features_default + 1
      x
    },
    rna_normalization = function(x) {
      x$preprocess$variable_features <- x$preprocess$variable_features + 1
      x
    },
    azimuth = function(x) { x$azimuth$score_threshold <- x$azimuth$score_threshold / 2; x },
    clustering = function(x) { x$annotation$cluster_seed <- x$annotation$cluster_seed + 1; x },
    annotation = function(x) {
      resolutions <- as.numeric(unlist(x$annotation$resolutions))
      current <- as.numeric(x$annotation$final_resolution)
      x$annotation$final_resolution <- resolutions[resolutions != current][[1]]
      x
    },
    proteomics = function(x) { x$proteomics$fdr <- x$proteomics$fdr / 2; x }
  )
  for (name in names(checks)) {
    changed <- scientific_config_projections(checks[[name]](config))
    changed_names <- names(baseline)[vapply(
      names(baseline),
      function(projection) {
        !identical(baseline[[projection]], changed[[projection]])
      },
      logical(1)
    )]
    stopifnot(identical(changed_names, name))
  }
  unrelated <- config
  unrelated$reporting <- list(note = "dependency isolation check")
  stopifnot(identical(baseline, scientific_config_projections(unrelated)))
  TRUE
}

validate_pipeline_config_dependencies <- function(pipeline) {
  target_names <- vapply(
    pipeline, function(target) target$settings$name, character(1)
  )
  direct_dependencies <- function(symbol) {
    target_names[vapply(
      pipeline,
      function(target) symbol %in% all.names(target$command$expr),
      logical(1)
    )]
  }
  analysis_allowlist <- c(
    "config_isolation_check", "rna_creation_config", "rna_filter_config",
    "rna_normalization_config", "azimuth_config", "clustering_config",
    "annotation_config",
    "lookup_file",
    "qc_thresholds_file", "raw_rna_dir", "annotation_file",
    "annotation_marker_file", "nk_overrides_file", "proteomics_parameters"
  )
  trust4_allowlist <- c(
    "trust4_discovery_config", "trust4_main_config", "trust4_ic_config",
    "trust4_patient_config", "sural_metadata_umap_file",
    "sural_ic_metadata_umap_file", "sural_trust4_manifest"
  )
  unexpected_analysis <- setdiff(
    direct_dependencies("analysis_config"), analysis_allowlist
  )
  unexpected_trust4 <- setdiff(
    direct_dependencies("trust4_config"), trust4_allowlist
  )
  if (length(unexpected_analysis) || length(unexpected_trust4)) {
    stop(
      "Broad config dependency detected: ",
      paste(c(unexpected_analysis, unexpected_trust4), collapse = ", ")
    )
  }
  invisible(pipeline)
}

configure_runtime <- function() {
  workers <- as.integer(Sys.getenv("ANALYSIS_WORKERS", unset = "6"))
  future_globals_gib <- as.numeric(
    Sys.getenv("FUTURE_GLOBALS_GIB", unset = "16")
  )
  stopifnot(
    length(workers) == 1L, is.finite(workers), workers > 0L,
    length(future_globals_gib) == 1L, is.finite(future_globals_gib),
    future_globals_gib > 0
  )
  future::plan("multicore", workers = workers)
  options(
    future.globals.maxSize = future_globals_gib * 1024^3
  )
  invisible(workers)
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}
