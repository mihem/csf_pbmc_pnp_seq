olink_assays <- function() {
  c(
    "CCL2", "CCL3", "CXCL8", "GZMA", "KLRD1", "SIRT2", "TNFSF14",
    "TREM2", "VEGFA", "IFNG", "CXCL12", "CXCL10", "TNFRSF9", "CCL7",
    "TNFSF12", "IL1B", "IL18"
  )
}

read_olink_input <- function(path, value_column, assays = olink_assays()) {
  result <- readxl::read_xlsx(path)
  stopifnot(all(c("SampleID", "Assay", value_column) %in% names(result)))
  result[[value_column]] <- as.numeric(result[[value_column]])
  result <- dplyr::filter(result, .data$Assay %in% .env$assays)
  stopifnot(nrow(result) > 0L, all(unique(result$Assay) %in% assays))
  result
}

read_olink_metadata <- function(path) {
  result <- suppressMessages(readxl::read_xlsx(path))
  required <- c("SampleID", "orbis_id", "diagnosis", "age", "sex")
  stopifnot(all(required %in% names(result)), !anyDuplicated(result$SampleID))
  result
}

prepare_olink_quantified <- function(data, metadata) {
  joined <- dplyr::left_join(data, metadata, by = "SampleID")
  stopifnot(!anyNA(joined$orbis_id), !anyNA(joined$diagnosis))
  wide <- tidyr::pivot_wider(
    joined,
    id_cols = tidyselect::all_of(c(
      "SampleID", "orbis_id", "diagnosis", "age", "sex"
    )),
    names_from = "Assay", values_from = "Quantified_value"
  )
  wide$diagnosis <- factor(
    wide$diagnosis, levels = c("CTRL", "GBS", "CIDP")
  )
  assays <- sort(intersect(unique(joined$Assay), names(wide)))
  stopifnot(length(assays) > 0L)
  list(data_wide = wide, assays = assays)
}
