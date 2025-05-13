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

# get DE genes of all clusters combined
de_combined <- read_xlsx(
    file.path("results", "de", "de_cidp_ctrl_csf_combined.xlsx")
)

conditions <- c(
    "cidp_ctrl_csf",
    "gbs_ctrl_csf",
    "cidp_ctrl_pbmc",
    "gbs_ctrl_pbmc"
)

de_top_combined_list <-
    lapply(conditions, read_de_combined_top) |>
    setNames(conditions)

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

plot_enrichment_results(go_gsea, prefix = "go_gsea", fold_change = ranked_genes)

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
