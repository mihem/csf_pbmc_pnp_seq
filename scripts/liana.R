##################################################
# cellular interaction analysis using LIANA
# requires running annotate.R first
##################################################
# Load required libraries ----
library(liana)
library(tidyverse)
library(Seurat)
library(qs)

# LIANA analysis ----
# Load and preprocess data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

# downsample sc_merge
sc_merge_small <- subset(sc_merge, downsample = 1000)

# Interaction analysis setup ----
# Get consensus interactions and filter for CXCL-related pairs
consensus_interactions <- select_resource("Consensus")$Consensus

olink_pairs <- consensus_interactions |>
  dplyr::filter(
    source_genesymbol %in%
      olink_markers |
      target_genesymbol %in% olink_markers
  )


# Perform cell-cell communication analysis using multiple methods
# using only CXCL interactions
liana_results <- liana_wrap(
  sc_merge_small,
  method = c("natmi", "connectome", "logfc", "sca"),
  resource = "custom",
  external_resource = olink_pairs,
  expression_pct = 0.1,
  parallel = TRUE,
  interactions = olink_pairs
)

qs::qsave(liana_results, file.path("objects", "liana_results.qs"))

# Results processing and visualization ----
liana_results$natmi |>
  print(width = Inf)

liana_results_aggregate <-
  liana_results |>
  liana_aggregate()


liana_results_aggregate |>
  #   dplyr::filter(ligand.complex == "CXCL14") |>
  print(width = Inf)


# Generate and save dotplot visualization
liana_plot_all <-
  liana_results_aggregate |>
  liana_dotplot() +
  theme(
    axis.text.x = element_text(size = 14, face = "plain"),
    strip.text = element_text(size = 14),
  )

ggsave(
  file.path("results", "interaction", "liana_all_dotplot.pdf"),
  liana_plot_all,
  width = 100,
  height = 10,
  limitsize = FALSE
)

liana_plot_CD8TEM3_target <-
  liana_results_aggregate |>
  liana_dotplot(target_groups = c("CD8TEM_3")) +
  theme(
    axis.text.x = element_text(size = 14, face = "plain"),
    strip.text = element_text(size = 14),
  )

ggsave(
  file.path("results", "interaction", "liana_cd8tem3_target_dotplot.pdf"),
  liana_plot_CD8TEM3_target,
  width = 40,
  height = 7,
  limitsize = FALSE
)
liana_plot_CD8TEM3_source <-
  liana_results_aggregate |>
  liana_dotplot(source_groups = c("CD8TEM_3")) +
  theme(
    axis.text.x = element_text(size = 14, face = "plain"),
    strip.text = element_text(size = 14),
  )

ggsave(
  file.path("results", "interaction", "liana_cd8tem3_source_dotplot.pdf"),
  liana_plot_CD8TEM3_source,
  width = 20,
  height = 5,
  limitsize = FALSE
)

selected_olink_markers <- c(
  "CC7",
  "CXCL8",
  "CCL3",
  "CCL2",
  "CXCL10",
  "TREM2",
  "IFNG",
  "TNFSF14"
)

liana_plot_selected_markers <-
  liana_results_aggregate |>
  dplyr::filter(
    ligand.complex %in%
      selected_olink_markers |
      receptor.complex %in% selected_olink_markers
  ) |>
  liana_dotplot() +
  theme(
    axis.text.x = element_text(
      size = 10,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      face = "plain"
    ),
    strip.text = element_text(size = 10),
  )

ggsave(
  file.path("results", "interaction", "liana_selected_markers_dotplot.pdf"),
  liana_plot_selected_markers,
  width = 50,
  height = 7,
  limitsize = FALSE
)

# selected ligands and receptors that are validated interactions
selected_ligands <- c("APOE", "CCL3", "TNFSF14")
selected_receptors <- c("TREM2", "CCR4", "CCR1", "TNFSFRSF14", "LTB")

liana_dotplot_selected <-
  liana_results_aggregate |>
  dplyr::filter(
    ligand.complex %in%
      selected_ligands |
      receptor.complex %in% selected_receptors
  ) |>
  liana_dotplot() +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(
    axis.text.x = element_text(
      size = 10,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      colour = "black",
      face = "plain"
    ),
    axis.text.y = element_text(colour = "black", size = 14),
    strip.text = element_text(size = 14, angle = 90),
  )

ggsave(
  file.path("results", "interaction", "liana_selected_dotplot.pdf"),
  liana_dotplot_selected,
  width = 14,
  height = 10
)

# Chord diagram filtered for selected interactions
pdf(
  file.path("results", "interaction", "liana_chord_freq_selected.pdf"),
  width = 20,
  height = 20
)
liana_results_aggregate |>
  dplyr::filter(
    ligand.complex %in%
      selected_ligands |
      receptor.complex %in% selected_receptors
  ) |>
  chord_freq()
dev.off()

# # Visualize CXCR4 expression ----
# dplot_cxcr4 <-
#   DotPlot(sc_merge, features = c("CXCR4"), dot.min = 0.01) +
#   viridis::scale_color_viridis(option = "viridis") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic")) +
#   xlab(NULL) +
#   ylab(NULL) +
#   scale_size(range = c(0, 5))

# ggsave(
#   file.path("results", "liana", "cxcr4_dotplot.pdf"),
#   dplot_cxcr4,
#   width = 3.5,
#   height = 5
# )
