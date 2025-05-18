###############################################################################
# perform enrichment analysis using clusterProfiler
# requires running annotate.R first
###############################################################################

# Required libraries ----
library(tidyverse)
library(ggplot2)
library(qs)
library(Seurat)
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(viridis)
library(writexl)

# Load helper functions ----
source(file.path("scripts", "enrichment_helper.R"))

# Load and prepare data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# Extract CD8/NK markers ----
topmarkers_cd8_nk <- read_xlsx(
    file.path("results", "de", "topmarkers.xlsx"),
    sheet = "CD8_NK"
)

# Get significant top markers with Entrez IDs
cd8_nk_markers_filtered <- topmarkers_cd8_nk |>
    mutate(entrez_id = map_to_entrez(gene)) |>
    filter(!is.na(entrez_id), p_val_adj < 0.05, avg_log2FC > 1) |>
    slice_min(order_by = p_val_adj, n = 200)

cd8_nk_logfc <- cd8_nk_markers_filtered$avg_log2FC
names(cd8_nk_logfc) <- cd8_nk_markers_filtered$entrez_id
cd8_nk_logfc <- sort(cd8_nk_logfc, decreasing = TRUE)

cd8_nk_logfc_gene <- cd8_nk_markers_filtered$avg_log2FC
names(cd8_nk_logfc_gene) <- cd8_nk_markers_filtered$gene
cd8_nk_logfc_gene <- sort(cd8_nk_logfc_gene, decreasing = TRUE)

# Read DE genes of all clusters combined
comparisons <- c(
    "cidp_ctrl_csf",
    "gbs_ctrl_csf",
    "cidp_ctrl_pbmc",
    "gbs_ctrl_pbmc"
)

de_top_combined_list <-
    lapply(comparisons, read_de_combined_top) |>
    setNames(comparisons)

# Read DE genes of specific clusters
# Define analysis configurations
de_cluster_parameters <- list(
    list(
        condition = "cidp_ctrl_csf",
        clusters = c("CD4TCM_2", "CD8_NK")
    ),
    list(
        condition = "gbs_ctrl_csf",
        clusters = c("pDC", "CD8_NK")
    ),
    list(
        condition = "cidp_ctrl_pbmc",
        clusters = c("CD4TCM_2", "NKCD56dim", "CD8_NK")
    ),
    list(
        condition = "gbs_ctrl_pbmc",
        clusters = c("Plasma", "CD16Mono")
    )
)

de_top_cluster_list <-
    lapply(
        de_cluster_parameters,
        function(config) {
            cluster_results <- lapply(
                config$clusters,
                function(cluster) {
                    read_de_cluster_top(
                        condition = config$condition,
                        sheets = cluster
                    )
                }
            )
            names(cluster_results) <- config$clusters
            cluster_results
        }
    )
names(de_top_cluster_list) <- sapply(de_cluster_parameters, `[[`, "condition")

# Prepare background gene set
background_genes <- map_to_entrez(rownames(sc_merge))
background_genes <- background_genes[!is.na(background_genes)]

# Prepare GSEA input ----
# Get all CD8/NK markers for ranking
all_cd8_nk_markers <- FindMarkers(
    object = sc_merge,
    ident.1 = "CD8_NK",
    only.pos = FALSE,
    min.pct = 0.1,
    logfc.threshold = 0,
    assay = "RNA"
) |>
    rownames_to_column("gene") |>
    mutate(entrez_id = map_to_entrez(gene)) |>
    filter(!is.na(entrez_id))

# Create ranked gene list for GSEA
# CD8_NK
ranked_genes <- all_cd8_nk_markers$avg_log2FC
names(ranked_genes) <- all_cd8_nk_markers$entrez_id
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Get all DE genes of all clusters combined
de_combined_all <-
    lapply(
        X = comparisons,
        FUN = function(comparison) {
            read_xlsx(
                file.path(
                    "results",
                    "de",
                    paste0("de_", comparison, "_combined.xlsx")
                )
            ) |>
                dplyr::mutate(entrez_id = map_to_entrez(gene)) |>
                dplyr::filter(!is.na(entrez_id))
        }
    ) |>
    setNames(comparisons)

# Create ranked gene list for de combined all
ranked_genes_combined_all <- lapply(
    de_combined_all,
    function(x) {
        ranked_genes <- x$avg_log2FC
        names(ranked_genes) <- x$entrez_id
        ranked_genes <- sort(ranked_genes, decreasing = TRUE)
        return(ranked_genes)
    }
) |>
    setNames(names(de_combined_all))

# Gene Ontology Analysis ----
# Over-representation analysis (ORA)
cd8_nk_go_ora <- enrichGO(
    gene = cd8_nk_markers_filtered$entrez_id,
    universe = background_genes,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
)

write_xlsx(
    data.frame(cd8_nk_go_ora),
    file.path("results", "enrich", "cd8_nk_go_ora_results.xlsx")
)

plot_enrichment_results(
    cd8_nk_go_ora,
    fold_change = cd8_nk_logfc,
    prefix = "cd8_nk_go_ora",
    width_dp = 6,
    height_dp = 6,
    width_hm = 7,
    height_hm = 4
)

# Create and plot GO term similarity tree
cd8_nk_go_ora_sim <- enrichplot::pairwise_termsim(cd8_nk_go_ora)
cd8_nk_tree_plot <- enrichplot::treeplot(cd8_nk_go_ora_sim, showCategory = 10)
save_plot(cd8_nk_tree_plot, "cd8_nk_go_ora_tree.pdf", width = 12, height = 8)

# Run GO enrichment analysis for all conditions
de_go_ora <- lapply(
    names(de_top_combined_list),
    function(condition) {
        run_go_enrichment(de_top_combined_list[[condition]], condition)
    }
) |>
    setNames(names(de_top_combined_list))

# Run GO enrichment analysis cluster specific DEG
de_go_ora_cluster <- lapply(
    names(de_top_cluster_list),
    function(condition) {
        cluster_results <- de_top_cluster_list[[condition]]
        # For each cluster in this condition
        cluster_enrichment <- lapply(
            names(cluster_results),
            function(cluster) {
                # Get DEG for this specific cluster
                deg_cluster <- cluster_results[[cluster]]
                # Run enrichment and name results with condition_cluster
                result <- run_go_enrichment(
                    deg_cluster,
                    paste0(condition, "_", cluster)
                )
                return(result)
            }
        ) |>
            setNames(names(cluster_results))
        return(cluster_enrichment)
    }
) |>
    setNames(names(de_top_cluster_list))

# Gene Set Enrichment Analysis (GSEA) GO
go_gsea <- gseGO(
    geneList = ranked_genes,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = TRUE,
    seed = TRUE
)

plot_enrichment_results(go_gsea, prefix = "go_gsea", fold_change = ranked_genes)

# GSEO GO all clusters combined
de_combined_gsea <- lapply(
    names(ranked_genes_combined_all),
    function(x) {
        gseGO(
            geneList = ranked_genes_combined_all[[x]],
            OrgDb = org.Hs.eg.db,
            ont = "BP",
            minGSSize = 10,
            maxGSSize = 500,
            pvalueCutoff = 0.05,
            verbose = TRUE,
            seed = TRUE
        )
    }
) |>
    setNames(names(ranked_genes_combined_all))

de_combined_gsea_readable <-
    lapply(
        de_combined_gsea,
        function(x) {
            setReadable(
                x,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID"
            )
        }
    )

lapply(
    names(de_combined_gsea_readable),
    function(x) {
        data.frame(de_combined_gsea_readable[[x]])
    }
) |>
    setNames(names(de_combined_gsea_readable)) |>
    write_xlsx(
        file.path("results", "enrich", "de_combined_gsea_results.xlsx")
    )


# Plot GSEA results for all clusters combined
lapply(
    names(de_combined_gsea),
    function(x) {
        if (nrow(de_combined_gsea[[x]]) == 0) {
            return(NULL)
        } else {
            plot_enrichment_results(
                de_combined_gsea[[x]],
                prefix = paste0(x, "_go_gsea"),
                fold_change = ranked_genes_combined_all[[x]]
            )
        }
    }
)

# Cell Marker Enrichment Analysis ----
cell_markers <- read_xlsx(file.path("lookup", "Cell_marker_Human.xlsx")) |>
    dplyr::select(cell_name, GeneID)

cd8_nk_cell_markers_enrichment <- enricher(
    cd8_nk_markers_filtered$entrez_id,
    TERM2GENE = cell_markers
)

cd8_nk_cell_markers_enrichment_readable <- setReadable(
    cd8_nk_cell_markers_enrichment,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
)

plot_enrichment_results(
    cd8_nk_cell_markers_enrichment_readable,
    prefix = "cd8_nk_cell_markers",
    fold_change = cd8_nk_logfc,
    width_dp = 6,
    height_dp = 3,
    width_hm = 7,
    height_hm = 3
)

write_xlsx(
    data.frame(cd8_nk_cell_markers_enrichment_readable),
    file.path("results", "enrich", "cd8_nk_cell_markers_results.xlsx")
)

# MSigDB Analysis ----
# Get MSigDB gene sets
msigdb_c8 <- msigdbr(species = "Homo sapiens", category = "C8") |>
    dplyr::select(gs_name, gene_symbol)

# Perform enrichment analysis for each MSigDB category
cd8_nk_msigdb_c8 <- enricher(
    cd8_nk_markers_filtered$gene,
    TERM2GENE = msigdb_c8,
    universe = rownames(sc_merge)
)

plot_enrichment_results(
    cd8_nk_msigdb_c8,
    prefix = "cd8_nk_msigdb_c8",
    fold_change = cd8_nk_logfc_gene,
    width_dp = 7,
    height_dp = 6,
    width_hm = 17,
    height_hm = 5
)

write_xlsx(
    data.frame(cd8_nk_msigdb_c8),
    file.path("results", "enrich", "cd8_nk_msigdb_c8_results.xlsx")
)
