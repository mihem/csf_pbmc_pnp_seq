discover_rna_inputs <- function(raw_rna_dir) {
  h5 <- list.files(
    raw_rna_dir,
    pattern = "raw_feature_bc_matrix_filtered\\.h5$",
    recursive = TRUE,
    full.names = TRUE
  )
  vireo <- list.files(
    raw_rna_dir,
    pattern = "donor_ids\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  metrics <- list.files(
    raw_rna_dir,
    pattern = "metrics_summary\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  h5_library <- basename(dirname(h5))
  vireo_library <- basename(dirname(dirname(vireo)))
  metrics_library <- basename(dirname(dirname(metrics)))

  manifest <- tibble::tibble(
    library_id = h5_library,
    h5 = h5,
    tissue = dplyr::case_when(
      stringr::str_detect(h5_library, "^CSF") ~ "CSF",
      stringr::str_detect(h5_library, "^PBMC") ~ "PBMC",
      h5_library == "PNP38" ~ "CSF",
      TRUE ~ NA_character_
    ),
    pool = dplyr::case_when(
      h5_library == "PNP38" ~ 10L,
      TRUE ~ as.integer(stringr::str_extract(h5_library, "[0-9]+$"))
    )
  ) |>
    dplyr::left_join(
      tibble::tibble(library_id = vireo_library, vireo = vireo),
      by = "library_id"
    ) |>
    dplyr::left_join(
      tibble::tibble(library_id = metrics_library, metrics = metrics),
      by = "library_id"
    ) |>
    dplyr::arrange(factor(tissue, levels = c("CSF", "PBMC")), pool)

  stopifnot(
    nrow(manifest) == 19L,
    !anyDuplicated(manifest$library_id),
    sum(!is.na(manifest$vireo)) == 18L,
    all(!is.na(manifest$metrics)),
    identical(manifest$library_id[manifest$pool == 10L], "PNP38")
  )
  manifest
}

rna_input_files <- function(manifest) {
  files <- c(manifest$h5, stats::na.omit(manifest$vireo), manifest$metrics)
  stopifnot(all(file.exists(files)))
  unname(files)
}

manual_annotation_files <- function(config) {
  paths <- c(
    config$paths$annotations,
    config$paths$markers,
    config$paths$nk_overrides,
    file.path("lookup", paste0(c("cl9", "cl12", "cl23", "cl24"), "_outliers.csv"))
  )
  stopifnot(all(file.exists(paths)))
  paths
}
