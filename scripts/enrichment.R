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

# Helper functions ----
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
    if (prefix != "")
        save_plot(
            dp,
            paste0(prefix, "_dotplot.pdf"),
            width = width_dp,
            height = height_dp
        )

    # Create heatplot
    hp <- heatplot(
        enrichment_result,
        showCategory = n_categories,
        foldChange = fold_change
    ) +
        viridis::scale_fill_viridis()
    if (prefix != "")
        save_plot(
            hp,
            paste0(prefix, "_heatplot.pdf"),
            width = width_hm,
            height = height_hm
        )
}

map_to_entrez <- function(genes, from_type = "SYMBOL") {
    entrez_ids <- mapIds(
        org.Hs.eg.db,
        keys = genes,
        keytype = from_type,
        column = "ENTREZID"
    )
    return(entrez_ids)
}

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
    filter(!is.na(entrez_id), p_val_adj < 0.05, avg_log2FC > 2) |>
    slice_min(order_by = p_val_adj, n = 200)

cd8_nk_logfc <- cd8_nk_markers_filtered$avg_log2FC
names(cd8_nk_logfc) <- cd8_nk_markers_filtered$entrez_id
cd8_nk_logfc <- sort(cd8_nk_logfc, decreasing = TRUE)

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
ranked_genes <- all_cd8_nk_markers$avg_log2FC
names(ranked_genes) <- all_cd8_nk_markers$entrez_id
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Gene Ontology Analysis ----
# Over-representation analysis (ORA)
go_ora <- enrichGO(
    gene = cd8_nk_markers_filtered$entrez_id,
    universe = background_genes,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
)

writexl::write_xlsx(
    data.frame(go_ora),
    file.path("results", "enrich", "cd8_nk_go_ora_results.xlsx")
)

plot_enrichment_results(
    go_ora,
    fold_change = cd8_nk_logfc,
    prefix = "cd8_nk_go_ora",
    width_dp = 6,
    height_dp = 6,
    width_hm = 7,
    height_hm = 4
)

# Create and plot GO term similarity tree
go_ora_sim <- enrichplot::pairwise_termsim(go_ora)
tree_plot <- enrichplot::treeplot(go_ora_sim, showCategory = 10)
save_plot(tree_plot, "cd8_nk_go_ora_tree.pdf", width = 12, height = 8)

# Gene Set Enrichment Analysis (GSEA)
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

plot_enrichment_results(go_gsea, prefix = "go_gsea")

# Cell Marker Enrichment Analysis ----
cell_markers <- read_xlsx(file.path("lookup", "Cell_marker_Human.xlsx")) |>
    select(cell_name, GeneID)

cell_enrichment <- enricher(
    cd8_nk_markers_filtered$entrez_id,
    TERM2GENE = cell_markers
)

plot_enrichment_results(cell_enrichment, prefix = "cell_markers")

# MSigDB Analysis ----
# Get MSigDB gene sets
msigdb_sets <- list(
    C7 = msigdbr(species = "Homo sapiens", category = "C7") |>
        select(gs_name, gene_symbol),
    C8 = msigdbr(species = "Homo sapiens", category = "C8") |>
        select(gs_name, gene_symbol)
)

# Perform enrichment analysis for each MSigDB category
msigdb_results <- map(
    msigdb_sets,
    ~ enricher(
        cd8_nk_markers_filtered$gene,
        TERM2GENE = .x,
        universe = rownames(sc_merge)
    )
)

# Plot MSigDB results
walk2(
    msigdb_results,
    names(msigdb_results),
    ~ plot_enrichment_results(.x, prefix = paste0("msigdb_", .y))
)
