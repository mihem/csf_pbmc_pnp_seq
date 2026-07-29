#===============================================================================
# Cerebro Export Script
#===============================================================================
# Purpose: Prepare and export Seurat object for CerebroApp visualization:
# - Convert Seurat object to Cerebro-compatible format
# - Clean up unnecessary data and optimize memory usage
# - Export selected metadata and marker genes
# - Generate interactive visualization file
#===============================================================================

# load libraries ----
library(Seurat)
library(tidyverse)
library(qs)
library(cerebroAppLite)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))

slots <- slotNames(sc_merge_cerebro)
for (slot in slots) {
    size <- object.size(slot(sc_merge_cerebro, slot))
    print(paste0("Slot ", slot, ": ", format(size, units = "auto")))
}

# prepare data for cerebro ----
sc_merge_cerebro <- sc_merge

# Remove unnecessary data to optimize memory usage
sc_merge_cerebro <- DietSeurat(
    sc_merge_cerebro,
    layers = c("data"),
    assays = "RNA",
    dimreducs = c("umap.stacas.ss.all")
)

# Additional cleanup
sc_merge_cerebro[["RNA"]]$scale.data <- NULL
sc_merge_cerebro[["RNA"]]$counts <- NULL

format(object.size(sc_merge_cerebro), units = "auto")

# Prepare metadata for export ----
# Select and rename relevant metadata columns
sc_merge_cerebro@meta.data <-
    sc_merge_cerebro@meta.data |>
    tibble::rownames_to_column("barcode") |>
    select(
        barcode,
        cluster,
        diagnosis,
        tissue,
        nCount_RNA,
        nFeature_RNA,
        percent_mt,
        sample,
        patient
    ) |>
    tibble::column_to_rownames(var = "barcode")

# Calculate marker genes ----
sc_merge_cerebro <- getMarkerGenes(
    sc_merge_cerebro,
    organism = "hg",
    groups = c("cluster"),
    only_pos = TRUE,
    assay = "RNA",
    min_pct = 0.1,
    thres_logFC = 0.25,
    thres_p_val = 0.05
)

# Sort marker genes by avg_log2FC descending ----
sc_merge_cerebro@misc$marker_genes$cerebro_seurat$cluster <-
    sc_merge_cerebro@misc$marker_genes$cerebro_seurat$cluster |>
    dplyr::arrange(dplyr::desc(avg_log2FC))

# Export for CerebroApp ----
# Save processed object and export to Cerebro format
qs::qsave(sc_merge_cerebro, file = file.path("objects", "sc_merge_cerebro.qs"))

cerebroAppLite::exportFromSeurat(
    object = sc_merge_cerebro,
    file = file.path("objects", "sc_merge_cerebro.crb"),
    experiment_name = "IN-seq",
    groups = c("cluster", "sample", "tissue", "diagnosis"),
    organism = "hg",
    nUMI = "nCount_RNA",
    nGene = "nFeature_RNA",
    use_delayed_array = FALSE
)

# Launch CerebroApp ----
cerebroAppLite::launchCerebro()
