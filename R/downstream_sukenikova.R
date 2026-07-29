sukenikova_result_dir <- function() {
  file.path("results", "targets", "sukenikova")
}

read_sukenikova_reactive_clones <- function(path) {
  clones <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot(
    "cdr3b_aa" %in% names(clones),
    nrow(clones) > 0L,
    !anyNA(clones$cdr3b_aa),
    all(nzchar(clones$cdr3b_aa))
  )
  clones |>
    dplyr::rename(CTaa_TRB = cdr3b_aa) |>
    dplyr::mutate(diagnosis = "GBS") |>
    dplyr::rename_with(
      function(column) paste0("sukenikova_", column),
      .cols = -CTaa_TRB
    )
}

extract_study_trb_clones <- function(sc_tcr) {
  required <- c("patient", "cluster", "CTaa", "CTgene", "CTnt", "tissue", "diagnosis")
  stopifnot(
    inherits(sc_tcr, "Seurat"),
    all(required %in% colnames(sc_tcr[[]]))
  )

  sc_tcr[[]] |>
    dplyr::select(dplyr::all_of(required)) |>
    tidyr::separate_wider_delim(CTaa, names = c("TRA", "TRB"), delim = "_") |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::rename(CTaa_TRB = TRB, CTaa_TRA = TRA) |>
    dplyr::rename_with(
      function(column) paste0("heming_", column),
      .cols = -CTaa_TRB
    )
}

match_sukenikova_reactive_clones <- function(study_clones, reactive_clones) {
  stopifnot(
    "CTaa_TRB" %in% names(study_clones),
    "CTaa_TRB" %in% names(reactive_clones)
  )
  result <- study_clones |>
    dplyr::inner_join(reactive_clones, by = "CTaa_TRB") |>
    dplyr::relocate(CTaa_TRB, .before = dplyr::everything())
  stopifnot(nrow(result) > 0L, !anyNA(result$CTaa_TRB))
  result
}

format_sukenikova_reactive_table <- function(shared_clones) {
  required <- c(
    "CTaa_TRB", "heming_patient", "heming_tissue", "heming_diagnosis",
    "sukenikova_pt", "sukenikova_diagnosis", "sukenikova_source",
    "sukenikova_specificity", "sukenikova_clone_id"
  )
  stopifnot(all(required %in% names(shared_clones)))

  shared_clones |>
    dplyr::arrange(heming_tissue, heming_diagnosis, CTaa_TRB) |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::rename(
      "CDR3 beta" = CTaa_TRB,
      "Patient" = heming_patient,
      "Tissue" = heming_tissue,
      "Diagnosis" = heming_diagnosis,
      "Sukenikova Patient" = sukenikova_pt,
      "Sukenikova Diagnosis" = sukenikova_diagnosis,
      "Sukenikova Tissue" = sukenikova_source,
      "Sukenikova Specificity" = sukenikova_specificity,
      "Sukenikova Clone ID" = sukenikova_clone_id
    )
}

write_sukenikova_workbook <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(data, path)
  stopifnot(file.exists(path), file.info(path)$size > 0L)
  path
}

write_sukenikova_table_pdf <- function(data, path) {
  stopifnot(nrow(data) > 0L)
  table_plot <- gridExtra::tableGrob(
    data,
    rows = NULL,
    theme = gridExtra::ttheme_default(
      core = list(
        fg_params = list(cex = 0.8),
        bg_params = list(fill = c("white", "lightgray"), alpha = 0.5)
      ),
      colhead = list(
        fg_params = list(cex = 0.9, fontface = "bold"),
        bg_params = list(fill = "darkgray", alpha = 0.8)
      )
    )
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = 12, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.draw(table_plot)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  stopifnot(file.exists(path), file.info(path)$size > 0L)
  path
}
