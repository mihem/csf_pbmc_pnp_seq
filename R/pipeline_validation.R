compare_pipeline_target_hashes <- function(
  baseline_path, exclusions = c("analysis_config", "analysis_config_file")
) {
  baseline <- readRDS(baseline_path) |>
    dplyr::select("name", baseline_hash = "data")
  current <- targets::tar_meta(fields = c("name", "parent", "data"))
  active_names <- targets::tar_manifest(fields = "name")$name
  repeat {
    expanded <- union(
      active_names,
      current$name[!is.na(current$parent) & current$parent %in% active_names]
    )
    if (identical(sort(active_names), sort(expanded))) break
    active_names <- expanded
  }
  current <- current |>
    dplyr::filter(.data$name %in% .env$active_names) |>
    dplyr::select(-"parent") |>
    dplyr::rename(current_hash = "data")
  dplyr::full_join(baseline, current, by = "name") |>
    dplyr::filter(!.data$name %in% .env$exclusions) |>
    dplyr::mutate(
      identical = !is.na(.data$baseline_hash) & !is.na(.data$current_hash) &
        .data$baseline_hash == .data$current_hash
    ) |>
    dplyr::arrange(.data$identical, .data$name)
}

validate_pipeline_equivalence <- function(baseline_path, exclusions = NULL) {
  default_exclusions <- c("analysis_config", "analysis_config_file")
  comparison <- compare_pipeline_target_hashes(
    baseline_path,
    exclusions = unique(c(default_exclusions, exclusions))
  )
  changed <- dplyr::filter(comparison, !.data$identical)
  if (nrow(changed)) {
    stop(
      "Scientific target hashes changed: ",
      paste(changed$name, collapse = ", ")
    )
  }
  comparison
}
