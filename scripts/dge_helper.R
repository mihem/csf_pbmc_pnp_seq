########################################################
# differential gene expression analysis helper functions
########################################################

runLimma <- function(seurat_object, cluster, lookup, condition1, condition2) {
  pseudobulk_data <- Libra::to_pseudobulk(
    seurat_object,
    cell_type_col = "cluster",
    label_col = "diagnosis",
    replicate_col = "patient"
  )

  dge <- edgeR::DGEList(
    counts = pseudobulk_data[[cluster]],
    group = colnames(pseudobulk_data[[cluster]])
  )
  count_check <- edgeR::cpm(dge) > 1
  keep <- which(rowSums(count_check) > 2)
  dge <- dge[keep, ]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  meta_limma <-
    data.frame(
      patient = gsub(x = colnames(dge), pattern = ":.+", replacement = "")
    ) |>
    dplyr::left_join(lookup, by = "patient")
  designMat <- model.matrix(~ 0 + diagnosis + sex + age, data = meta_limma)
  my_contrasts <- glue::glue("diagnosis{condition1}-diagnosis{condition2}")
  my_args <- list(my_contrasts, levels = designMat)
  my_contrasts <- do.call(limma::makeContrasts, my_args)
  dge_voom <- limma::voomWithQualityWeights(dge, designMat, plot = FALSE)
  dge_voom <- dge_voom |>
    limma::lmFit(design = designMat, block = NULL) |>
    limma::contrasts.fit(my_contrasts) |>
    limma::eBayes(robust = TRUE)

  # topgenes_sig <- limma::topTable(dge_voom, n = Inf, adjust.method = "BH") |>
  #     dplyr::filter(adj.P.Val < 0.05) |>
  #     tibble::rownames_to_column("gene") |>
  #     tibble::tibble() |>
  #     dplyr::arrange(desc(logFC)) |>
  #     dplyr::rename(avg_log2FC = logFC, p_val_adj = adj.P.Val)

  topgenes_all <- limma::topTable(dge_voom, n = Inf, adjust.method = "BH") |>
    tibble::rownames_to_column("gene") |>
    tibble::tibble() |>
    dplyr::arrange(desc(logFC)) |>
    dplyr::rename(avg_log2FC = logFC, p_val_adj = adj.P.Val)

  return(topgenes_all)
}

# Create error-safe version of runLimma using purrr
safe_runLimma <- purrr::possibly(runLimma, otherwise = NULL)

# plot number of DEG per cluster ---
plotDE <- function(name, title) {
  sheets <- readxl::excel_sheets(
    path = file.path("results", "de", paste0(name, "_cluster.xlsx"))
  )
  cl_sig <-
    lapply(
      sheets,
      function(sheet) {
        read_xlsx(
          path = file.path("results", "de", paste0(name, "_cluster.xlsx")),
          sheet = sheet
        ) |>
          dplyr::filter(p_val_adj < 0.05) |>
          nrow()
      }
    )
  result <- tibble(
    cluster = sheets,
    n = unlist(cl_sig)
  )

  plot <-
    result |>
    mutate(cluster = fct_reorder(cluster, n)) |>
    ggplot(aes(x = cluster, y = n, fill = cluster)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = sc_merge@misc$cluster_col) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(
      x = "",
      y = "",
      title = title
    )
  ggsave(
    plot = plot,
    filename = file.path("results", "de", paste0(name, ".pdf")),
    width = 3,
    height = 3
  )
}

# Function to run limma on combined clusters
runLimmaCombined <- function(seurat_object, lookup, condition1, condition2) {
  pseudobulk_data <- Libra::to_pseudobulk(
    seurat_object,
    cell_type_col = "orig.ident",
    label_col = "diagnosis",
    replicate_col = "patient"
  )

  dge <- edgeR::DGEList(
    counts = pseudobulk_data$PNP,
    group = colnames(pseudobulk_data$PNP)
  )

  # Filter low expressed genes
  count_check <- edgeR::cpm(dge) > 1
  keep <- which(rowSums(count_check) > 2)
  dge <- dge[keep, ]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")

  # Create metadata for limma
  meta_limma <-
    data.frame(
      patient = gsub(x = colnames(dge), pattern = ":.+", replacement = "")
    ) |>
    dplyr::left_join(lookup, by = "patient")

  # Create design matrix
  designMat <- model.matrix(~ 0 + diagnosis + sex + age, data = meta_limma)

  # Create contrast matrix
  my_contrast <- paste0("diagnosis", condition1, "-diagnosis", condition2)
  contrast_matrix <- limma::makeContrasts(
    contrasts = my_contrast,
    levels = designMat
  )

  # Voom transformation and fitting
  dge_voom <- limma::voomWithQualityWeights(dge, designMat, plot = FALSE)

  dge_voom <- dge_voom |>
    limma::lmFit(design = designMat, block = NULL) |>
    limma::contrasts.fit(contrast_matrix) |>
    limma::eBayes(robust = TRUE)

  # Get all results
  topgenes_df <- limma::topTable(dge_voom, n = Inf, adjust.method = "BH")
  topgenes_all <- tibble::rownames_to_column(topgenes_df, "gene")
  topgenes_all <- dplyr::arrange(topgenes_all, -logFC)
  topgenes_all <- dplyr::rename(
    topgenes_all,
    avg_log2FC = logFC,
    p_val_adj = adj.P.Val
  )

  return(topgenes_all)
}

# Create a general function for differential expression analysis
performDEAnalysis <- function(
  seurat_object,
  condition1,
  condition2,
  tissue_type
) {
  # Create descriptive name for output files
  comparison_name <- paste0(
    "de_",
    tolower(condition1),
    "_",
    tolower(condition2),
    "_",
    tolower(tissue_type)
  )

  # Subset data for the specific conditions and tissue
  sc_subset <- subset(
    seurat_object,
    subset = diagnosis %in%
      c(condition1, condition2) &
      tissue == tissue_type
  )
  sc_subset$diagnosis <- droplevels(sc_subset$diagnosis)

  # Verify tissue filtering is correct
  unique_tissue <- unique(sc_subset$tissue)
  stopifnot(unique_tissue == tissue_type)

  # Verify condition filtering is correct
  unique_diagnoses <- unique(as.character(sc_subset$diagnosis))
  stopifnot(all(c(condition1, condition2) %in% unique_diagnoses))

  # Sanity check
  message(paste(
    "Performing DE analysis:",
    condition1,
    "vs",
    condition2,
    "in",
    tissue_type
  ))

  # Run DE analysis for each cluster
  de_results <- purrr::map(
    levels(sc_subset),
    function(cluster) {
      result <- safe_runLimma(
        cluster = cluster,
        seurat_object = sc_subset,
        condition1 = condition1,
        condition2 = condition2,
        lookup = lookup
      )
      if (is.null(result)) {
        message(paste("Failed to run DE for cluster:", cluster))
      }
      return(result)
    }
  )
  names(de_results) <- levels(sc_subset)

  # Remove NULL results
  de_results <- purrr::compact(de_results)

  # Save results to Excel
  writexl::write_xlsx(
    de_results,
    file.path("results", "de", paste0(comparison_name, "_cluster.xlsx"))
  )

  # Add combined cluster analysis
  combined_result <- runLimmaCombined(
    seurat_object = sc_subset,
    lookup = lookup,
    condition1 = condition1,
    condition2 = condition2
  )
  combined_result <- list(
    combined = combined_result
  )

  writexl::write_xlsx(
    combined_result,
    file.path("results", "de", paste0(comparison_name, "_combined.xlsx"))
  )

  # Plot number of DEGs per cluster
  plotDE(
    comparison_name,
    title = paste(condition1, "vs", condition2, tissue_type)
  )

  return(de_results)
}
