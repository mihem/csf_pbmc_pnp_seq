########################################################
# differential gene expression analysis using pseudobulk
# requires running annotate.R first
########################################################

# libraries  ----
library(Seurat)
library(tidyverse)
library(qs)
library(limma)
library(DESeq2)
library(viridis)
library(Libra)
library(EnhancedVolcano)
library(readxl)
library(purrr)
library(writexl)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

lookup <-
    read_excel(file.path("lookup", "SEED_lookup_v6.xlsx")) |>
    janitor::clean_names() |>
    mutate(age = lubridate::time_length(difftime(date, birth_date), "years")) |>
    mutate(
        diagnosis = factor(
            diagnosis,
            levels = c(
                "CTRL",
                "CIAP",
                "CIDP",
                "GBS",
                "MAG",
                "MFS",
                "PNC",
                "CAN",
                "PPN"
            )
        )
    )

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
    my_contrasts <- do.call(makeContrasts, my_args)
    dge_voom <- limma::voomWithQualityWeights(dge, designMat, plot = FALSE)

    dge_voom <- dge_voom |>
        limma::lmFit(design = designMat, block = NULL) |>
        limma::contrasts.fit(my_contrasts) |>
        eBayes(robust = TRUE)

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
        path = file.path("results", "de", paste0(name, ".xlsx"))
    )
    cl_sig <-
        lapply(
            sheets,
            function(sheet) {
                read_xlsx(
                    path = file.path("results", "de", paste0(name, ".xlsx")),
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

    # Name the results with cluster names
    names(de_results) <- levels(sc_subset)

    # Remove NULL results
    de_results <- purrr::compact(de_results)

    # Save results to Excel
    writexl::write_xlsx(
        de_results,
        file.path("results", "de", paste0(comparison_name, ".xlsx"))
    )

    # Plot number of DEGs per cluster
    plotDE(
        comparison_name,
        title = paste(condition1, "vs", condition2, tissue_type)
    )

    return(de_results)
}

# Run DE analysis for CIDP vs CTRL in CSF
de_cidp_ctrl_csf <- performDEAnalysis(sc_merge, "CIDP", "CTRL", "CSF")

# Run DE analysis for GBS vs CTRL in CSF
de_gbs_ctrl_csf <- performDEAnalysis(sc_merge, "GBS", "CTRL", "CSF")

# Run DE analysis for CIDP vs CTRL in PBMC
de_cidp_ctrl_pbmc <- performDEAnalysis(sc_merge, "CIDP", "CTRL", "PBMC")

# Run DE analysis for GBS vs CTRL in PBMC
de_gbs_ctrl_pbmc <- performDEAnalysis(sc_merge, "GBS", "CTRL", "PBMC")

# volcano plot
volcanoPlot <- function(
    filename,
    sheet,
    FCcutoff = 2,
    selectLab = NULL,
    drawConnectors = TRUE,
    condition1,
    condition2
) {
    input <- readxl::read_excel(
        file.path("results", "de", paste0(filename, ".xlsx")),
        sheet = sheet
    )
    if (nrow(input) != 0) {
        volcano <- EnhancedVolcano::EnhancedVolcano(
            data.frame(input),
            lab = paste0("italic('", input[["gene"]], "')"),
            x = "avg_log2FC",
            y = "p_val_adj",
            xlim = c(min(input[["avg_log2FC"]], max(input[["avg_log2FC"]]))),
            ylim = c(0, max(-log10(input[["p_val_adj"]]))),
            pCutoff = 0.1,
            FCcutoff = FCcutoff,
            axisLabSize = 15,
            pointSize = 2,
            labSize = 5,
            subtitle = NULL,
            caption = NULL,
            border = "full",
            gridlines.major = FALSE,
            gridlines.minor = FALSE,
            drawConnectors = drawConnectors,
            lengthConnectors = grid::unit(0.0001, "npc"),
            title = paste(condition1, "vs", condition2, "in ", sheet),
            boxedLabels = TRUE,
            selectLab = selectLab,
            xlab = bquote(~ Log[2] ~ "fold change"),
            ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
            parseLabels = TRUE,
            legendLabels = c(
                "NS",
                "logFC",
                "p-val",
                "p-val + logFC"
            ),
            legendPosition = "right",
        )
        pdf(
            file.path("results", "de", paste0(filename, "_", sheet, ".pdf")),
            width = 8,
            height = 6
        )
        print(volcano)
        dev.off()
    }
}

# Define analysis configurations
volcano_parameters <- list(
    list(
        filename = "de_cidp_ctrl_csf",
        clusters = c("CD4TCM_2", "CD8_NK"),
        condition1 = "CIDP",
        condition2 = "CTRL"
    ),
    list(
        filename = "de_gbs_ctrl_csf",
        clusters = c("pDC", "CD8_NK"),
        condition1 = "GBS",
        condition2 = "CTRL"
    ),
    list(
        filename = "de_cidp_ctrl_pbmc",
        clusters = c("CD4TCM_2", "NKCD56dim", "CD8_NK"),
        condition1 = "CIDP",
        condition2 = "CTRL"
    ),
    list(
        filename = "de_gbs_ctrl_pbmc",
        clusters = c("Plasma", "CD16Mono"),
        condition1 = "GBS",
        condition2 = "CTRL"
    )
)

# Common labels for plots if needed
lab_pnp_ctrl <- list(
    "mySC" = paste0(
        "italic('",
        c("DCN", "TNXB", "COL1A1", "COL15A1", "CD53", "IL4R", "CD74"),
        "')"
    ),
    "nmSC" = paste0(
        "italic('",
        c("IL10RA", "IL13RA1", "CSF2RA", "TGFBI"),
        "')"
    ),
    "repairSC" = paste0("italic('", c("GALR1", "TMEM47"), "')"),
    "PC2" = paste0(
        "italic('",
        c("MFAP5", "NLGN4Y", "PCDH11Y", "IFIT3", "OASL", "MX1"),
        "')"
    )
)

# Generate volcano plots for all configurations
lapply(volcano_parameters, function(config) {
    lapply(
        config$clusters,
        function(cluster) {
            volcanoPlot(
                filename = config$filename,
                sheet = cluster,
                FCcutoff = 2,
                condition1 = config$condition1,
                condition2 = config$condition2,
                selectLab = NULL 
            )
        }
    )
})

# function to compare gene expression between conditions using VlnPlot -----
compareGeneExpression <- function(seu_obj, gene, seu_obj_name) {
    plot <- VlnPlot(
        seu_obj,
        features = gene,
        group.by = "level2",
        cols = seu_obj@misc$level2_cols,
        pt.size = 0
    ) +
        NoLegend() +
        xlab("") +
        ylab("") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(
        plot = plot,
        filename = file.path(
            "results",
            "de",
            paste0(seu_obj_name, "_", gene, ".pdf")
        ),
        width = 3,
        height = 3
    )
}

# compare CXCL14 expression in perneurial cells between conditions
perineurial <- subset(
    sc_merge,
    cluster %in%
        c("periC1", "periC2", "periC3") &
        level2 %in% c("CTRL", "CIDP", "VN", "CIAP")
)
perineurial$level2 <- factor(
    perineurial$level2,
    levels = c("CTRL", "CIDP", "VN", "CIAP")
)
compareGeneExpression(perineurial, "CXCL14", "perineurial")

# compare GRIK3 and PRIMA1 expression in nmSC between conditions
nmSC <- subset(
    sc_merge,
    cluster %in% c("nmSC") & level2 %in% c("CTRL", "CIDP", "VN", "CIAP")
)
nmSC$level2 <- factor(nmSC$level2, levels = c("CTRL", "CIDP", "VN", "CIAP"))

lapply(
    c("GRIK3", "PRIMA1"),
    function(gene) compareGeneExpression(nmSC, gene, "nmSC")
)

# compare MLIP expression in mySC between conditions
mySC <- subset(
    sc_merge,
    cluster %in% c("mySC") & level2 %in% c("CTRL", "CIDP", "VN", "CIAP")
)
mySC$level2 <- factor(mySC$level2, levels = c("CTRL", "CIDP", "VN", "CIAP"))
compareGeneExpression(mySC, "MLIP", "mySC")
