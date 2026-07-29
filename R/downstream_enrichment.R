enrichment_cluster_parameters <- function() {
  list(
    cidp_ctrl_csf = c("CD4TCM_2", "CD8TEM_3"),
    gbs_ctrl_csf = "pDC",
    cidp_ctrl_pbmc = c("CD4TCM_2", "NKCD56dim", "CD8TEM_3"),
    gbs_ctrl_pbmc = c("Plasma", "CD16Mono")
  )
}

collect_enrichment_deg_results <- function(
  combined_results,
  cluster_results,
  cluster_parameters = enrichment_cluster_parameters()
) {
  comparisons <- deg_comparisons()$comparison
  stopifnot(
    identical(names(combined_results), comparisons),
    identical(names(cluster_results), comparisons),
    identical(names(cluster_parameters), comparisons),
    all(vapply(combined_results, is.data.frame, logical(1))),
    all(vapply(cluster_results, is.list, logical(1)))
  )

  result <- stats::setNames(vector("list", length(comparisons)), comparisons)
  for (comparison in comparisons) {
    clusters <- cluster_parameters[[comparison]]
    stopifnot(all(clusters %in% names(cluster_results[[comparison]])))
    result[[comparison]] <- list(
      combined = combined_results[[comparison]],
      cluster = cluster_results[[comparison]][clusters]
    )
  }
  result
}

map_enrichment_entrez <- function(genes) {
  stopifnot(is.character(genes), !anyNA(genes))
  if (!length(genes)) {
    return(character())
  }
  database <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
  unname(AnnotationDbi::mapIds(
    database,
    keys = genes,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  ))
}

enrichment_background_genes <- function(object) {
  stopifnot(inherits(object, "Seurat"), nrow(object) > 0L)
  genes <- map_enrichment_entrez(rownames(object))
  sort(unique(genes[!is.na(genes)]))
}

find_cd8tem3_markers <- function(object) {
  stopifnot(inherits(object, "Seurat"), "cluster" %in% colnames(object[[]]))
  Seurat::Idents(object) <- object$cluster
  Seurat::FindMarkers(
    object = object,
    ident.1 = "CD8TEM_3",
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    assay = "RNA"
  ) |>
    tibble::rownames_to_column("gene") |>
    dplyr::filter(.data$p_val_adj < 0.05) |>
    dplyr::arrange(dplyr::desc(.data$avg_log2FC))
}

run_cd8tem3_marker_ora <- function(markers, background_genes) {
  selected <- markers |>
    dplyr::mutate(entrez_id = map_enrichment_entrez(.data$gene)) |>
    dplyr::filter(!is.na(.data$entrez_id), .data$avg_log2FC > 1) |>
    dplyr::slice_min(order_by = .data$p_val_adj, n = 200)
  database <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
  clusterProfiler::enrichGO(
    gene = selected$entrez_id,
    universe = background_genes,
    OrgDb = database,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

write_cd8tem3_marker_ora <- function(
  result,
  root = file.path(enrichment_result_dir(), "cd8tem_3")
) {
  stopifnot(nrow(as.data.frame(result)) > 0L)
  workbook <- file.path(root, "cd8tem_3_go_ora_results.xlsx")
  dotplot <- file.path(root, "cd8tem_3_go_ora_dotplot.pdf")
  ensure_parent_dir(workbook)
  writexl::write_xlsx(as.data.frame(result), workbook)
  save_enrichment_plot(
    enrichplot::dotplot(result, showCategory = 10),
    dotplot,
    width = 6,
    height = 6
  )
  c(workbook, dotplot)
}

make_enrichment_analyses <- function(deg_results) {
  analyses <- list()
  for (comparison in names(deg_results)) {
    analyses[[comparison]] <- list(
      comparison = comparison,
      cluster = NA_character_,
      result = deg_results[[comparison]]$combined
    )
    for (cluster in names(deg_results[[comparison]]$cluster)) {
      name <- paste(comparison, cluster, sep = "__")
      analyses[[name]] <- list(
        comparison = comparison,
        cluster = cluster,
        result = deg_results[[comparison]]$cluster[[cluster]]
      )
    }
  }
  stopifnot(length(analyses) > 0L, !anyDuplicated(names(analyses)))
  analyses
}

prepare_enrichment_rank <- function(result) {
  stopifnot(
    is.data.frame(result),
    all(c("gene", "avg_log2FC") %in% names(result)),
    !anyNA(result$gene),
    !anyNA(result$avg_log2FC)
  )
  entrez <- map_enrichment_entrez(result$gene)
  keep <- !is.na(entrez)
  ranked <- result$avg_log2FC[keep]
  names(ranked) <- entrez[keep]
  sort(ranked, decreasing = TRUE)
}

prepare_enrichment_ora_genes <- function(result, n = 100L) {
  stopifnot(
    is.data.frame(result),
    all(c("gene", "avg_log2FC", "p_val_adj") %in% names(result)),
    length(n) == 1L,
    n > 0L
  )
  filtered <- result[
    !is.na(result$p_val_adj) & result$p_val_adj < 0.05 &
      !is.na(result$avg_log2FC) & abs(result$avg_log2FC) > 1,
    ,
    drop = FALSE
  ]
  filtered$entrez_id <- map_enrichment_entrez(filtered$gene)
  filtered <- filtered[!is.na(filtered$entrez_id), , drop = FALSE]

  select_direction <- function(direction) {
    selected <- filtered[sign(filtered$avg_log2FC) == direction, , drop = FALSE]
    selected <- selected[order(selected$p_val_adj, seq_len(nrow(selected))), , drop = FALSE]
    utils::head(selected, n)
  }
  list(down = select_direction(-1), up = select_direction(1))
}

run_enrichment_ora <- function(
  deg_results,
  background_genes,
  minimum_genes = 10L
) {
  analyses <- make_enrichment_analyses(deg_results)
  database <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
  results <- list()
  for (analysis_name in names(analyses)) {
    genes_by_direction <- prepare_enrichment_ora_genes(
      analyses[[analysis_name]]$result
    )
    for (direction in names(genes_by_direction)) {
      genes <- genes_by_direction[[direction]]$entrez_id
      result <- NULL
      if (length(genes) >= minimum_genes) {
        result <- clusterProfiler::enrichGO(
          gene = genes,
          universe = background_genes,
          OrgDb = database,
          ont = "BP",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.01,
          qvalueCutoff = 0.05,
          readable = TRUE
        )
      }
      name <- paste(analysis_name, direction, sep = "__")
      results[[name]] <- c(
        analyses[[analysis_name]][c("comparison", "cluster")],
        list(direction = direction, result = result)
      )
    }
  }
  results
}

run_enrichment_gsea <- function(deg_results, seed = 123L) {
  stopifnot(length(seed) == 1L, !is.na(seed))
  analyses <- make_enrichment_analyses(deg_results)
  database <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
  lapply(analyses, function(analysis) {
    ranked <- prepare_enrichment_rank(analysis$result)
    result <- withr::with_seed(seed, clusterProfiler::gseGO(
      geneList = ranked,
      OrgDb = database,
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      verbose = FALSE,
      seed = TRUE
    ))
    result <- clusterProfiler::setReadable(
      result,
      OrgDb = database,
      keyType = "ENTREZID"
    )
    c(
      analysis[c("comparison", "cluster")],
      list(result = result, ranked = ranked)
    )
  })
}

enrichment_result_dir <- function(root = "results/targets/enrichment") {
  root
}

enrichment_scope <- function(cluster) {
  if (is.na(cluster)) "de_combined" else "de_cluster"
}

enrichment_prefix <- function(comparison, cluster = NA_character_) {
  if (is.na(cluster)) comparison else paste(comparison, cluster, sep = "_")
}

save_enrichment_plot <- function(plot, path, width, height) {
  ensure_parent_dir(path)
  ggplot2::ggsave(path, plot = plot, width = width, height = height, dpi = 300)
  path
}

write_enrichment_plots <- function(result, ranked, prefix, directory) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0L) {
    return(character())
  }
  dotplot <- enrichplot::dotplot(result, showCategory = 10)
  heatplot <- suppressMessages(
    enrichplot::heatplot(
      result,
      showCategory = 10,
      foldChange = ranked
    ) + ggplot2::scale_fill_viridis_c()
  )
  c(
    save_enrichment_plot(
      dotplot,
      file.path(directory, paste0(prefix, "_dotplot.pdf")),
      7,
      6
    ),
    save_enrichment_plot(
      heatplot,
      file.path(directory, paste0(prefix, "_heatplot.pdf")),
      12,
      8
    )
  )
}

write_enrichment_ora_outputs <- function(
  results,
  root = enrichment_result_dir()
) {
  paths <- character()
  for (item in results) {
    if (is.null(item$result) || nrow(as.data.frame(item$result)) == 0L) {
      next
    }
    directory <- file.path(
      root,
      enrichment_scope(item$cluster),
      "ora"
    )
    prefix <- paste0(
      enrichment_prefix(item$comparison, item$cluster),
      "_",
      item$direction,
      "_go_ora"
    )
    workbook <- file.path(directory, paste0(prefix, ".xlsx"))
    ensure_parent_dir(workbook)
    writexl::write_xlsx(as.data.frame(item$result), workbook)
    paths <- c(
      paths,
      workbook,
      write_enrichment_plots(item$result, NULL, prefix, directory)
    )
  }
  stopifnot(length(paths) == 0L || all(file.exists(paths)))
  paths
}

write_enrichment_gsea_outputs <- function(
  results,
  root = enrichment_result_dir()
) {
  paths <- character()
  combined <- results[vapply(results, function(x) is.na(x$cluster), logical(1))]
  combined_file <- file.path(
    root,
    "de_combined",
    "gsea",
    "de_combined_gsea_results.xlsx"
  )
  ensure_parent_dir(combined_file)
  writexl::write_xlsx(
    stats::setNames(
      lapply(combined, function(x) as.data.frame(x$result)),
      vapply(combined, `[[`, character(1), "comparison")
    ),
    combined_file
  )
  paths <- c(paths, combined_file)

  cluster <- results[!vapply(results, function(x) is.na(x$cluster), logical(1))]
  for (comparison in unique(vapply(cluster, `[[`, character(1), "comparison"))) {
    selected <- cluster[vapply(
      cluster,
      function(x) identical(x$comparison, comparison),
      logical(1)
    )]
    workbook <- file.path(
      root,
      "de_cluster",
      "gsea",
      paste0("de_", comparison, "_cluster_gsea_results.xlsx")
    )
    ensure_parent_dir(workbook)
    writexl::write_xlsx(
      stats::setNames(
        lapply(selected, function(x) as.data.frame(x$result)),
        vapply(selected, `[[`, character(1), "cluster")
      ),
      workbook
    )
    paths <- c(paths, workbook)
  }

  for (item in results) {
    directory <- file.path(
      root,
      enrichment_scope(item$cluster),
      "gsea"
    )
    prefix <- paste0(
      enrichment_prefix(item$comparison, item$cluster),
      "_go_gsea"
    )
    paths <- c(
      paths,
      write_enrichment_plots(item$result, item$ranked, prefix, directory)
    )
  }
  stopifnot(length(paths) > 0L, all(file.exists(paths)), !anyDuplicated(paths))
  paths
}
