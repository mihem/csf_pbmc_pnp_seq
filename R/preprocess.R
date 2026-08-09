clean_lookup <- function(path) {
  readxl::read_excel(path) |>
    janitor::clean_names() |>
    dplyr::filter(cohort %in% c("scRNA", "scRNA_flow")) |>
    dplyr::mutate(
      incat_progress =
        (incat_follow_up - incat_at_lumbar_puncture) / follow_up,
      onls_progress =
        (onls_follow_up - onls_at_lumbar_puncture) / follow_up,
      mrc_sum_score_60_progress =
        (mrc_sum_score_60_follow_up -
          mrc_sum_score_60_at_lumbar_puncture) / follow_up,
      group = factor(group, levels = c("CTRL", "PNP")),
      group2 = factor(group2, levels = c("CTRL", "NIN", "IN", "other")),
      diagnosis = factor(
        diagnosis,
        levels = c(
          "CTRL", "CIAP", "CIDP", "GBS", "MAG",
          "MFS", "PNC", "CAN", "PPN"
        )
      ),
      dplyr::across(dml_median_motoric:ncv_sural_sensory, as.numeric)
    )
}

read_donor_assignments <- function(manifest, input_files) {
  stopifnot(length(input_files) > 0L)
  multiplexed <- dplyr::filter(manifest, !is.na(vireo))
  assignments <- purrr::map(
    multiplexed$vireo,
    ~ readr::read_tsv(.x, show_col_types = FALSE)
  )
  names(assignments) <- multiplexed$library_id

  recode <- c(
    "1B" = "PNP02",
    "2B" = "PNP03",
    "3B" = "PNP06",
    "4B" = "PNP09",
    "doublet" = "doublet",
    "unassigned" = "unassigned"
  )
  for (library_id in c("CSF_pool_1", "PBMC_pool_1")) {
    assignments[[library_id]]$donor_id <- unname(
      recode[assignments[[library_id]]$donor_id]
    )
  }
  stopifnot(all(vapply(assignments, function(x) !anyDuplicated(x$cell), logical(1))))
  assignments
}

create_sample_objects <- function(
  manifest,
  lookup,
  donor_assignments,
  settings,
  input_files
) {
  stopifnot(length(input_files) > 0L)
  configure_runtime()
  pool_objects <- purrr::map2(
    manifest$h5,
    manifest$library_id,
    function(path, library_id) {
      counts <- scMisc::ReadCellBender_h5(path)
      Seurat::CreateSeuratObject(
        counts = counts,
        min.cells = as.integer(settings$min_cells),
        min.features = as.integer(settings$min_features),
        project = "PNP"
      )
    }
  )
  names(pool_objects) <- manifest$library_id

  pseudonym_lookup <- lookup |>
    dplyr::select(pseudonym, patient) |>
    dplyr::add_row(pseudonym = "unassigned", patient = "unassigned") |>
    dplyr::add_row(pseudonym = "doublet", patient = "doublet")

  sample_objects <- list()
  for (library_id in setdiff(manifest$library_id, "PNP38")) {
    object <- pool_objects[[library_id]]
    assignments <- donor_assignments[[library_id]]
    assignment_index <- match(colnames(object), assignments$cell)
    object$pseudonym <- assignments$donor_id[assignment_index]
    object@meta.data <- object@meta.data |>
      tibble::rownames_to_column("barcode") |>
      dplyr::left_join(pseudonym_lookup, by = "pseudonym") |>
      tibble::column_to_rownames("barcode")

    split_objects <- Seurat::SplitObject(object, split.by = "patient")
    split_objects <- split_objects[
      !names(split_objects) %in% c("doublet", "unassigned")
    ]
    split_objects <- split_objects[!is.na(names(split_objects))]
    tissue <- manifest$tissue[match(library_id, manifest$library_id)]
    names(split_objects) <- paste(tissue, names(split_objects), sep = "_")
    sample_objects <- c(sample_objects, split_objects)
  }

  sample_objects$CSF_P14 <- pool_objects$PNP38
  sample_objects <- sample_objects[order(names(sample_objects))]
  stopifnot(
    length(sample_objects) == 61L,
    !anyDuplicated(names(sample_objects)),
    identical(names(sample_objects), sort(names(sample_objects)))
  )

  purrr::map(sample_objects, function(object) {
    object[["percent_mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^MT")
    object
  })
}

filter_sample_objects <- function(sample_objects, path, settings) {
  thresholds <- readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::rename(sample = patient, max_features = rna, max_percent_mt = mt)
  stopifnot(
    !anyDuplicated(thresholds$sample),
    setequal(names(sample_objects), thresholds$sample)
  )

  filtered <- purrr::imap(sample_objects, function(object, sample) {
    threshold <- dplyr::filter(thresholds, .data$sample == .env$sample)
    keep <- object$nFeature_RNA > settings$max_features_default &
      object$nFeature_RNA < threshold$max_features[[1]] &
      object$percent_mt < threshold$max_percent_mt[[1]]
    object[, keep]
  })
  stopifnot(sum(vapply(filtered, ncol, numeric(1))) == 119123L)
  filtered
}

merge_sample_objects <- function(filtered, lookup) {
  object <- merge(
    x = filtered[[1]],
    y = filtered[-1],
    merge.data = TRUE,
    add.cell.ids = names(filtered)
  )
  object$tissue <- stringr::str_extract(colnames(object), "(PBMC|CSF)")
  object$patient <- stringr::str_extract(colnames(object), "P\\d+")
  object$sample <- stringr::str_extract(colnames(object), "(PBMC|CSF)_P\\d+")
  object <- SeuratObject::JoinLayers(object)
  object <- split(x = object, f = object$sample)

  metadata <- lookup |>
    dplyr::select(
      patient,
      sex,
      group,
      diagnosis,
      incat_at_lumbar_puncture,
      incat_follow_up,
      onls_at_lumbar_puncture,
      onls_follow_up,
      mrc_sum_score_60_at_lumbar_puncture,
      mrc_sum_score_60_follow_up,
      icu,
      age
    )
  object@meta.data <- object@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(metadata, by = "patient") |>
    tibble::column_to_rownames("barcode")
  stopifnot(identical(as.numeric(dim(object)), c(27006, 119123)))
  object
}

normalize_and_reduce <- function(object, settings) {
  configure_runtime()
  object <- Seurat::NormalizeData(
    object,
    normalization.method = settings$normalization_method,
    scale.factor = settings$scale_factor,
    verbose = TRUE
  )
  object <- Seurat::FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = as.integer(settings$variable_features)
  )
  object <- Seurat::ScaleData(object)
  object <- Seurat::RunPCA(object, seed.use = 42L)

  nmf_object <- object
  nmf_object[["RNA3"]] <- methods::as(nmf_object[["RNA"]], "Assay")
  Seurat::DefaultAssay(nmf_object) <- "RNA3"
  nmf_object[["RNA"]] <- NULL
  nmf_object <- SeuratObject::RenameAssays(nmf_object, RNA3 = "RNA")
  nmf_object <- withr::with_seed(
    as.integer(settings$nmf_seed),
    singlet::RunNMF(
      nmf_object,
      k = as.integer(settings$nmf_factors)
    )
  )
  object[["nmf"]] <- nmf_object[["nmf"]]
  object
}

add_batch_metadata <- function(object, lookup) {
  batches <- lookup |>
    dplyr::select(patient, batch)
  object@meta.data <- object@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::select(-dplyr::any_of("batch")) |>
    dplyr::left_join(batches, by = "patient") |>
    tibble::column_to_rownames("barcode")
  object$batch <- paste0(object$batch, "_", object$tissue)
  stopifnot(length(unique(object$batch)) == 19L, !anyNA(object$batch))
  object
}

summarize_cell_counts <- function(sample_objects, merged_object) {
  tibble::tibble(
    sample = names(sample_objects),
    before = vapply(sample_objects, ncol, numeric(1))
  ) |>
    dplyr::left_join(
      dplyr::count(merged_object@meta.data, sample, name = "after"),
      by = "sample"
    )
}

summarize_gene_counts <- function(merged_object) {
  tibble::tibble(
    feature = merged_object$nFeature_RNA,
    sample = merged_object$sample
  ) |>
    dplyr::group_by(sample) |>
    dplyr::summarize(median_genes_after = stats::median(feature), .groups = "drop")
}

summarize_cellranger_metrics <- function(manifest) {
  purrr::map2_dfr(
    manifest$metrics,
    manifest$library_id,
    ~ readr::read_csv(.x, show_col_types = FALSE) |>
      dplyr::mutate(sample = .y, .before = 1)
  ) |>
    dplyr::arrange(sample)
}

write_csv_result <- function(data, path) {
  ensure_parent_dir(path)
  readr::write_csv(data, path)
  path
}
