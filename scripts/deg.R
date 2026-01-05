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

# functions ----
source(file.path("scripts", "dge_helper.R"))

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

lookup <- qs::qread(file.path("objects", "lookup.qs"), nthread = 4)

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
    filename = "de_cidp_ctrl_csf_cluster",
    clusters = c("CD4TCM_2", "CD8_NK"),
    condition1 = "CIDP",
    condition2 = "CTRL"
  ),
  list(
    filename = "de_cidp_ctrl_csf_combined",
    clusters = c("combined"),
    condition1 = "CIDP",
    condition2 = "CTRL"
  ),
  list(
    filename = "de_gbs_ctrl_csf_cluster",
    clusters = c("pDC", "CD8_NK"),
    condition1 = "GBS",
    condition2 = "CTRL"
  ),
  list(
    filename = "de_gbs_ctrl_csf_combined",
    clusters = c("combined"),
    condition1 = "GBS",
    condition2 = "CTRL"
  ),
  list(
    filename = "de_cidp_ctrl_pbmc_cluster",
    clusters = c("CD4TCM_2", "NKCD56dim", "CD8_NK"),
    condition1 = "CIDP",
    condition2 = "CTRL"
  ),
  list(
    filename = "de_cidp_ctrl_pbmc_combined",
    clusters = c("combined"),
    condition1 = "CIDP",
    condition2 = "CTRL"
  ),
  list(
    filename = "de_gbs_ctrl_pbmc_cluster",
    clusters = c("Plasma", "CD16Mono"),
    condition1 = "GBS",
    condition2 = "CTRL"
  ),
  list(
    filename = "de_gbs_ctrl_pbmc_combined",
    clusters = c("combined"),
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
