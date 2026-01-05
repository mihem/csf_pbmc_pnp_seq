# Enrichment helper functions ----
save_plot <- function(plot, filename, width, height) {
  dir.create(
    file.path("results", "enrich"),
    showWarnings = FALSE,
    recursive = TRUE
  )
  ggsave(
    file.path("results", "enrich", filename),
    plot,
    width = width,
    height = height,
    dpi = 300
  )
}

plot_enrichment_results <- function(
  enrichment_result,
  fold_change,
  n_categories = 10,
  prefix = "",
  width_dp = 12,
  height_dp = 8,
  width_hm = 12,
  height_hm = 8
) {
  # Create dotplot
  dp <- dotplot(enrichment_result, showCategory = n_categories)
  if (prefix != "") {
    save_plot(
      dp,
      paste0(prefix, "_dotplot.pdf"),
      width = width_dp,
      height = height_dp
    )
  }

  # Create heatplot
  hp <- heatplot(
    enrichment_result,
    showCategory = n_categories,
    foldChange = fold_change
  ) +
    viridis::scale_fill_viridis()
  if (prefix != "") {
    save_plot(
      hp,
      paste0(prefix, "_heatplot.pdf"),
      width = width_hm,
      height = height_hm
    )
  }
}

map_to_entrez <- function(genes, from_type = "SYMBOL") {
  entrez_ids <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = genes,
    keytype = from_type,
    column = "ENTREZID"
  )
  return(entrez_ids)
}

# Split DE genes into up and down regulated using base::split
read_de_combined_top <- function(condition) {
  # Read and filter for significant genes
  de <- readxl::read_xlsx(
    file.path("results", "de", paste0("de_", condition, "_combined.xlsx"))
  ) |>
    dplyr::filter(
      p_val_adj < 0.05,
      abs(avg_log2FC) > 1
    )

  if (nrow(de) == 0) {
    return(NULL)
  }

  # Add entrez IDs and split by direction
  de_split <- de |>
    dplyr::mutate(entrez_id = map_to_entrez(gene)) |>
    dplyr::filter(!is.na(entrez_id)) |>
    (function(x) split(x, f = sign(x$avg_log2FC)))() |>
    setNames(c("down", "up"))

  # Get top 100 for each direction by p-value
  de_split$up <- de_split$up |>
    dplyr::slice_min(order_by = p_val_adj, n = 100, with_ties = FALSE)
  de_split$down <- de_split$down |>
    dplyr::slice_min(order_by = p_val_adj, n = 100, with_ties = FALSE)

  return(de_split)
}

# Read DE results from cluster files with multiple sheets
read_de_cluster_top <- function(condition, sheets) {
  results <- list()

  for (sheet in sheets) {
    # Read and filter for significant genes
    de <- readxl::read_xlsx(
      file.path(
        "results",
        "de",
        paste0("de_", condition, "_cluster.xlsx")
      ),
      sheet = sheet
    ) |>
      dplyr::filter(
        p_val_adj < 0.05,
        abs(avg_log2FC) > 1
      )

    # Add entrez IDs and split by direction
    de_split <- de |>
      dplyr::mutate(entrez_id = map_to_entrez(gene)) |>
      dplyr::filter(!is.na(entrez_id))

    # Get top 100 for each direction by p-value
    nrows_up <- de_split |>
      dplyr::filter(avg_log2FC > 0) |>
      nrow()

    nrows_down <- de_split |>
      dplyr::filter(avg_log2FC < 0) |>
      nrow()

    if (nrows_down == 0) {
      results$down <- NULL
    } else {
      results$down <- de_split |>
        dplyr::filter(avg_log2FC < 0) |>
        dplyr::slice_min(
          order_by = p_val_adj,
          n = 100,
          with_ties = FALSE
        )
    }

    if (nrows_up == 0) {
      results$up <- NULL
    } else {
      results$up <- de_split |>
        dplyr::filter(avg_log2FC > 0) |>
        dplyr::slice_min(
          order_by = p_val_adj,
          n = 100,
          with_ties = FALSE
        )
    }
  }

  return(results)
}

# Helper function for GO enrichment analysis
run_go_enrichment <- function(gene_list, name, universe = background_genes) {
  # Function to process single direction
  process_direction <- function(genes, direction) {
    if (is.null(genes) || nrow(genes) == 0) {
      cat(sprintf("%s %s: No genes to analyze\n", name, direction))
      return(NULL)
    }

    # Check if we have enough genes for meaningful enrichment
    if (length(genes$entrez_id) < 10) {
      cat(sprintf(
        "%s %s: Only %d genes, skipping enrichment (minimum 10 required)\n",
        name,
        direction,
        length(genes$entrez_id)
      ))
      return(NULL)
    }

    prefix <- paste0(name, "_", direction, "_go_ora")

    # Run enrichment
    result <- enrichGO(
      gene = genes$entrez_id,
      universe = universe,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05,
      readable = TRUE
    )

    # If we have results, save and plot them
    if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
      n_terms <- nrow(as.data.frame(result))
      cat(sprintf("%s: %d enriched terms\n", prefix, n_terms))

      # Create plots
      plot_enrichment_results(
        result,
        fold_change = NULL,
        prefix = prefix,
        width_dp = 7,
        height_dp = 6
      )

      # Save results to Excel
      write_xlsx(
        data.frame(result),
        file.path("results", "enrich", paste0(prefix, ".xlsx"))
      )
    } else {
      cat(sprintf("%s: No enrichment results\n", prefix))
    }

    return(result)
  }

  # Process both directions if gene_list exists
  if (!is.null(gene_list)) {
    return(list(
      up = process_direction(gene_list$up, "up"),
      down = process_direction(gene_list$down, "down")
    ))
  } else {
    cat(sprintf("%s: No results to analyze\n", name))
    return(NULL)
  }
}
