run_azimuth <- function(object, config) {
  configure_runtime(config)
  # The pinned Azimuth reference calls the SeuratObject Key setter by name.
  suppressPackageStartupMessages(library("Seurat", character.only = TRUE))
  query <- SeuratObject::JoinLayers(object)
  mapped <- Azimuth::RunAzimuth(
    query = query,
    reference = config$azimuth$reference,
    assay = "RNA"
  )
  threshold <- as.numeric(config$azimuth$score_threshold)
  predictions <- data.frame(row.names = colnames(mapped))
  for (level in 1:3) {
    source_label <- paste0("predicted.celltype.l", level)
    source_score <- paste0(source_label, ".score")
    target_label <- paste0("azimuth_pbmcref", level)
    target_score <- paste0(target_label, "_score")
    labels <- as.character(mapped[[source_label, drop = TRUE]])
    scores <- mapped[[source_score, drop = TRUE]]
    labels[scores < threshold] <- "unknown"
    predictions[[target_label]] <- labels
    predictions[[target_score]] <- scores
  }
  Seurat::AddMetaData(query, predictions)
}
