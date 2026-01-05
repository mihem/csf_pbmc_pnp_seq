##################################################
# cellular interaction analysis using CellChat
# requires running annotate.R first
##################################################

#libraries ---
library(Seurat)
library(tidyverse)
library(qs)
library(CellChat)
library(tidyverse)
library(readxl)

# # set parallelization ---
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 8000 * 1024^2) # 8GB

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# cellchat requires sample names to be stored in "samples"
sc_merge$samples <- sc_merge$sample

# downsample sc_merge
sc_merge_small <- subset(sc_merge, downsample = 1000)

# read olink data ----
olink_quant_file <- file.path(
  "raw",
  "olink",
  "olink_quant_long_filtered.xlsx"
)

olink_quant <- read_xlsx(olink_quant_file) |>
  mutate(Quantified_value = as.numeric(Quantified_value))

olink_markers <- unique(olink_quant$Assay)

# create CellChat object ----
cellchat <- createCellChat(
  object = sc_merge_small,
  group.by = "cluster",
  assay = "RNA"
)

# # set the ligand receptor database ----
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# # use a subset of CellChatDB for cell-cell communication analysis
CellChatDB_use <- subsetDB(
  CellChatDB,
  search = "Secreted Signaling",
  key = "annotation"
)

cellchat@DB <- CellChatDB_use

# cellchatdb_olink_pre <-
#     CellChatDB.human$interaction |>
#     filter(ligand %in% olink_markers | receptor %in% olink_markers)

# cellchatdb_olink <-
#     updateCellChatDB(
#         db = cellchatdb_olink_pre,
#         species_target = "human",
#         merge = FALSE
#     )

# # # subset the expression data of signaling genes for saving computation cost
# cellchat@DB <- cellchatdb_olink

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

qsave(cellchat, file.path("objects", "cellchat.qs"))
cellchat <- qread(file.path("objects", "cellchat.qs"))

# compute communication probability ----
cellchat <- computeCommunProb(cellchat, type = "triMean")
# cellchat <- filterCommunication(cellchat, min.cells = 10)

# infer the cell-cell communication at signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


# Filter CellChat database for interactions involving Olink markers
# Get all ligand-receptor interactions from CellChat
cellchat_interactions <- cellchat@DB$interaction

unique(cellchat_interactions$interaction_name)

# Find interactions where ligand or receptor matches Olink markers
olink_interactions <- cellchat_interactions |>
  filter(ligand %in% olink_markers | receptor %in% olink_markers)

# Create pairLR.use dataframe with interaction names
pairLR_olink <- data.frame(
  interaction_name = olink_interactions$interaction_name
)

# Create bubble plot with Olink-related interactions
pdf(
  file.path("results", "interaction", "olink_bubble.pdf"),
  width = 10,
  height = 8
)
netVisual_bubble(
  cellchat,
  pairLR.use = pairLR_olink,
  remove.isolate = TRUE
)
dev.off()

# Create bubble plot with Olink-related interactions
pdf(
  file.path("results", "interaction", "olink_bubble.pdf"),
  width = 10,
  height = 8
)
netVisual_bubble(
  cellchat,
  remove.isolate = TRUE
)
dev.off()

levels(cellchat@idents)

pdf(
  file.path("results", "interaction", "cd8tem3_source_bubble.pdf"),
  width = 10,
  height = 8
)
netVisual_bubble(
  cellchat,
  sources.use = 13,
  remove.isolate = TRUE
)
dev.off()

pdf(
  file.path("results", "interaction", "cd8tem3_target_bubble.pdf"),
  width = 10,
  height = 8
)
netVisual_bubble(
  cellchat,
  targets.use = 13,
  remove.isolate = TRUE
)
dev.off()
