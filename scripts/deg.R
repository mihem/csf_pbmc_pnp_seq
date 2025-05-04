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

    topgenes_sig <- limma::topTable(dge_voom, n = Inf, adjust.method = "BH") |>
        dplyr::filter(adj.P.Val < 0.05) |>
        tibble::rownames_to_column("gene") |>
        tibble::tibble() |>
        dplyr::arrange(desc(logFC)) |>
        dplyr::rename(avg_log2FC = logFC, p_val_adj = adj.P.Val)

    return(topgenes_sig)
}

# CIDP vs CTRL pseudobulk ---
sc_cidp_ctrl_csf <- subset(
    sc_merge,
    subset = diagnosis %in% c("CIDP", "CTRL") & tissue %in% c("CSF")
)
sc_cidp_ctrl_csf$diagnosis <- droplevels(sc_cidp_ctrl_csf$diagnosis)

# sanity check
table(sc_cidp_ctrl_csf$diagnosis)
table(sc_cidp_ctrl_csf$tissue)
table(sc_cidp_ctrl_csf$cluster, sc_cidp_ctrl_csf$diagnosis)

# Create error-safe version of runLimma using purrr
safe_runLimma <- purrr::possibly(runLimma, otherwise = NULL)

de_cidp_ctrl_csf <-
    purrr::map(
        levels(sc_cidp_ctrl_csf),
        function(cluster) {
            result <- safe_runLimma(
                cluster = cluster,
                seurat_object = sc_cidp_ctrl_csf,
                condition1 = "CIDP",
                condition2 = "CTRL",
                lookup = lookup
            )
            if (is.null(result)) {
                message(paste("Failed to run DE for cluster:", cluster))
            }
            return(result)
        }
    )

# Name the results with cluster names
names(de_cidp_ctrl_csf) <- levels(sc_cidp_ctrl_csf)

# Remove NULL results
de_cidp_ctrl_csf <- purrr::compact(de_cidp_ctrl_csf)
writexl::write_xlsx(
    de_cidp_ctrl_csf,
    file.path("results", "de", "de_cidp_ctrl_csf.xlsx")
)


# AggregateExpression(
#     sc_cidp_ctrl_csf_cd8_nk,
#     assay = "RNA",
#     group.by = "diagnosis",
#     features = c("SYNE1", "ITGB1", "CD44")
# )

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
    ) |>
        dplyr::filter(n > 0) # Fixed: Filter applied to the data frame, not the vector

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

plotDE("de_cidp_ctrl_csf", title = "CIDP vs CTRL CSF")

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
            x = "avg_logFC",
            y = "p_val_adj",
            xlim = c(min(input[["avg_logFC"]], max(input[["avg_logFC"]]))),
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
            lengthConnectors = unit(0.0001, "npc"),
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

cluster_de <- c("repairSC", "mySC", "nmSC", "PC2")

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

# PNP vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "pnp_ctrl_pseudobulk",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "PNP",
            condition2 = "CTRL",
            selectLab = lab_pnp_ctrl[[cluster]]
        )
    }
)

# VN vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "vn_ctrl",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "VN",
            condition2 = "CTRL"
        )
    }
)

# ciap vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "ciap_ctrl",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "ciap",
            condition2 = "CTRL"
        )
    }
)
# CIAP vs CTRL
lapply(
    cluster_de,
    function(cluster) {
        volcanoPlot(
            filename = "ciap_ctrl",
            sheet = cluster,
            FCcutoff = 2,
            condition1 = "CIAP",
            condition2 = "CTRL"
        )
    }
)

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
