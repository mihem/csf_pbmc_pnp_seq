# GLIPH2 manuscript figures

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(readxl)
library(patchwork)
library(ggrepel)
library(pheatmap)
library(qs)

set.seed(42)

# ============================================================================
# Color palettes and ordering (must match main analysis)
# ============================================================================
sc_tcr <- qs::qread("sc_tcr.qs", nthreads = 6)
diagnosis_col <- sc_tcr@misc$diagnosis_col
rm(sc_tcr); gc()

tissue_col <- c("CSF" = "#E41A1C", "PBMC" = "#377EB8")
neuropathy_dx <- c("CIAP", "CIDP", "GBS", "MAG", "MFS", "PNC", "CAN", "PPN")
dx_order <- c("CTRL", neuropathy_dx)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "grey95")
    )
}

dir.create("results/tcr_comparison/manuscript_figures", showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Figure 1: Disease-enriched GLIPH cluster overview
# Bubble plot: cluster size vs FDR, colored by diagnosis, with CDR3 labels
# ============================================================================
message("Figure 1: Disease-enriched cluster overview")

tissue_dx <- readxl::read_xlsx("results/tcr_comparison/gliph/tissue_specificity_analysis.xlsx",
                                sheet = "cluster_tissue_dx")
cdr3_motifs <- readxl::read_xlsx("results/tcr_comparison/gliph/disease_enriched_cdr3_motifs.xlsx")

# Merge CDR3 examples into tissue_dx
cluster_plot_data <- tissue_dx |>
  dplyr::mutate(
    odds_ratio = as.numeric(odds_ratio),
    p.adj = as.numeric(p.adj),
    cluster_size = as.numeric(cluster_size),
    csf_fraction = as.numeric(csf_fraction),
    n_ctrl_in_cluster = as.numeric(n_ctrl_in_cluster)
  ) |>
  dplyr::filter(p.adj < 0.1) |>
  dplyr::left_join(
    cdr3_motifs |> dplyr::select(cluster, example_CDR3s),
    by = "cluster"
  ) |>
  dplyr::mutate(
    # Extract first CDR3 as label
    lead_CDR3 = gsub(",.*", "", example_CDR3s),
    lead_CDR3 = trimws(lead_CDR3),
    diagnosis = factor(diagnosis, levels = neuropathy_dx),
    tissue_bias = factor(tissue_bias, levels = c("CSF-enriched", "No bias", "PBMC-enriched")),
    neg_log10_fdr = -log10(p.adj + 1e-300),
    # Cap infinite OR for plotting
    plot_log2OR = pmin(log2(ifelse(is.infinite(odds_ratio), 1e6, odds_ratio) + 0.01), 12)
  )

# Label top 2 per diagnosis by cluster_size
top_labels <- cluster_plot_data |>
  dplyr::group_by(diagnosis) |>
  dplyr::slice_max(order_by = cluster_size, n = 2) |>
  dplyr::ungroup()

p1 <- ggplot(cluster_plot_data, aes(x = diagnosis, y = cluster_size)) +
  geom_jitter(aes(fill = diagnosis, size = neg_log10_fdr, shape = tissue_bias),
              width = 0.2, alpha = 0.75, stroke = 0.4) +
  scale_fill_manual(values = diagnosis_col, guide = "none") +
  scale_shape_manual(values = c("CSF-enriched" = 24, "No bias" = 21, "PBMC-enriched" = 25),
                     name = "Tissue compartment") +
  scale_size_continuous(range = c(2, 8), name = "-log10(FDR)",
                        breaks = c(2, 5, 10, 50)) +
  geom_text_repel(data = top_labels, aes(label = lead_CDR3),
                  size = 2.3, max.overlaps = 20, fontface = "italic",
                  segment.color = "grey50", segment.size = 0.3,
                  nudge_y = 5, seed = 42) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Disease-enriched GLIPH2 specificity groups",
    x = "", y = "Cluster size (# clonotypes)"
  )

ggsave("results/tcr_comparison/manuscript_figures/fig_gliph_clusters_overview.pdf",
       p1, width = 9, height = 6.5)
message("  Saved fig_gliph_clusters_overview.pdf")


# ============================================================================
# Figure 2: Driving CDR3 motifs — one plot per diagnosis
# ============================================================================
message("Figure 2: Driving CDR3 motifs per diagnosis (individual plots)")

motif_data <- readxl::read_xlsx("results/tcr_comparison/gliph/motif_enrichment_by_diagnosis.xlsx")

# Top 5 3-mer motifs per diagnosis by OR (FDR < 0.05)
top_motifs <- motif_data |>
  dplyr::filter(k == 3, p.adj < 0.05, odds_ratio > 1) |>
  dplyr::group_by(diagnosis) |>
  dplyr::slice_max(order_by = odds_ratio, n = 5) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

for (dx in intersect(neuropathy_dx, unique(top_motifs$diagnosis))) {
  dx_motifs <- top_motifs |> dplyr::filter(diagnosis == dx)

  p_dx <- ggplot(dx_motifs, aes(x = reorder(motif, log2_fc), y = log2_fc)) +
    geom_bar(stat = "identity", fill = diagnosis_col[dx], width = 0.7, alpha = 0.85) +
    geom_text(aes(label = signif), hjust = -0.3, size = 3.5) +
    coord_flip() +
    theme_pub(base_size = 11) +
    labs(
      title = paste0(dx, " — CDR3 motifs enriched in GLIPH2 clusters"),
      x = "CDR3 motif", y = "log2(fold change vs background)"
    )

  ggsave(paste0("results/tcr_comparison/manuscript_figures/fig_gliph_motifs_", dx, ".pdf"),
         p_dx, width = 7, height = 4.5)
  message("  Saved fig_gliph_motifs_", dx, ".pdf")
}


# ============================================================================
# Figure 3: Tissue compartmentalization of disease-enriched clusters
# Stacked bar + dot overlay
# ============================================================================
message("Figure 3: Tissue compartmentalization")

tissue_summary <- readxl::read_xlsx("results/tcr_comparison/gliph/tissue_specificity_analysis.xlsx",
                                     sheet = "tissue_bias_summary")

tissue_long <- tissue_summary |>
  tidyr::pivot_longer(cols = -diagnosis, names_to = "tissue_bias", values_to = "n_clusters") |>
  dplyr::mutate(
    diagnosis = factor(diagnosis, levels = neuropathy_dx),
    tissue_bias = factor(tissue_bias, levels = c("CSF-enriched", "No bias", "PBMC-enriched"))
  )

# Calculate totals for label
totals <- tissue_long |>
  dplyr::group_by(diagnosis) |>
  dplyr::summarize(total = sum(n_clusters), .groups = "drop")

p3 <- ggplot(tissue_long, aes(x = diagnosis, y = n_clusters, fill = tissue_bias)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(data = totals, aes(x = diagnosis, y = total, label = total, fill = NULL),
            vjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(
    values = c("CSF-enriched" = "#E41A1C", "No bias" = "grey70", "PBMC-enriched" = "#377EB8"),
    name = "Tissue compartment"
  ) +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Tissue compartmentalization of disease-enriched GLIPH2 clusters",
    x = "", y = "Number of disease-enriched clusters"
  )

ggsave("results/tcr_comparison/manuscript_figures/fig_gliph_tissue_bias.pdf",
       p3, width = 8, height = 5.5)
message("  Saved fig_gliph_tissue_bias.pdf")


# ============================================================================
# Figure 4: Tissue-specific CDR3 motifs (CSF vs PBMC clusters)
# ============================================================================
message("Figure 4: Tissue-specific motifs")

tissue_motifs <- readxl::read_xlsx("results/tcr_comparison/gliph/tissue_specific_motif_comparison.xlsx")

# Top tissue-differential 3-mers
top_tissue_motifs <- tissue_motifs |>
  dplyr::filter(k == 3, p.adj < 0.1) |>
  dplyr::arrange(desc(abs(log2_fc))) |>
  dplyr::slice(1:min(20, n())) |>
  dplyr::mutate(
    tissue_direction = ifelse(log2_fc > 0, "CSF-enriched", "PBMC-enriched"),
    signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

if (nrow(top_tissue_motifs) > 0) {
  p4 <- ggplot(top_tissue_motifs, aes(x = reorder(motif, log2_fc), y = log2_fc,
                                       fill = tissue_direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = signif),
              hjust = ifelse(top_tissue_motifs$log2_fc > 0, -0.3, 1.3), size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = c("CSF-enriched" = "#E41A1C", "PBMC-enriched" = "#377EB8"),
                      name = "") +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.3) +
    theme_pub(base_size = 11) +
    labs(
      title = "CDR3 motifs differentially enriched in CSF vs PBMC GLIPH2 clusters",
      x = "CDR3 3-mer motif",
      y = "log2(CSF frequency / PBMC frequency)"
    )

  ggsave("results/tcr_comparison/manuscript_figures/fig_gliph_tissue_motifs.pdf",
         p4, width = 8, height = 5.5)
  message("  Saved fig_gliph_tissue_motifs.pdf")
}


# ============================================================================
# Figure 5a: MAG V-gene usage — one plot per gene
# ============================================================================
message("Figure 5a: MAG V-gene usage (individual plots)")

# Read gene usage stats
gene_stats <- readxl::read_xlsx("results/tcr_comparison/gene_usage/gene_usage_statistics.xlsx")

# Get all TRA V gene results across diagnoses (to show MAG in context)
trav_data <- gene_stats |>
  dplyr::filter(chain == "TRA", gene_type == "V") |>
  dplyr::mutate(
    diagnosis = factor(diagnosis, levels = neuropathy_dx),
    significant = !is.na(p.adj) & p.adj < 0.05
  )

# Focus on the two MAG-significant genes
mag_sig_genes <- c("TRAV29/DV5", "TRAV38-2/DV8")

for (gene_name in mag_sig_genes) {
  gene_data <- trav_data |>
    dplyr::filter(gene == gene_name, !is.na(estimate))

  # Clean gene name for filename
  gene_file <- gsub("/", "-", gene_name)

  p_gene <- ggplot(gene_data, aes(x = diagnosis, y = estimate, color = diagnosis)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = -log10(p.adj + 1e-10)), alpha = 0.8) +
    geom_errorbar(aes(ymin = estimate - 1.96 * SE, ymax = estimate + 1.96 * SE),
                  width = 0.2, alpha = 0.5) +
    scale_color_manual(values = diagnosis_col, guide = "none") +
    scale_size_continuous(range = c(2, 6), name = "-log10(FDR)") +
    theme_pub(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = paste0(gene_name, " usage across neuropathy diagnoses"),
      x = "", y = "Difference from CTRL (proportion)"
    )

  ggsave(paste0("results/tcr_comparison/manuscript_figures/fig_mag_vgene_", gene_file, ".pdf"),
         p_gene, width = 6.5, height = 5)
  message("  Saved fig_mag_vgene_", gene_file, ".pdf")
}

# ============================================================================
# Figure 5b: MAG-specific GLIPH motifs (standalone)
# ============================================================================
message("Figure 5b: MAG GLIPH motifs")

mag_motifs <- motif_data |>
  dplyr::filter(diagnosis == "MAG", k == 3, p.adj < 0.05, odds_ratio > 1) |>
  dplyr::slice_max(order_by = odds_ratio, n = 10) |>
  dplyr::mutate(
    signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

p_mag_motifs <- ggplot(mag_motifs, aes(x = reorder(motif, log2_fc), y = log2_fc)) +
  geom_bar(stat = "identity", fill = diagnosis_col["MAG"], width = 0.7, alpha = 0.85) +
  geom_text(aes(label = signif), hjust = -0.3, size = 3.5) +
  coord_flip() +
  theme_pub(base_size = 11) +
  labs(
    title = "MAG — CDR3 motifs enriched in GLIPH2 clusters",
    x = "CDR3 motif",
    y = "log2(fold change vs background)"
  )

ggsave("results/tcr_comparison/manuscript_figures/fig_mag_gliph_motifs.pdf",
       p_mag_motifs, width = 7, height = 5.5)
message("  Saved fig_mag_gliph_motifs.pdf")


# ============================================================================
# Figure 6: Disease-specific cluster summary table as plot
# Top clusters per diagnosis with CDR3, size, tissue bias
# ============================================================================
message("Figure 6: Summary table figure")

# Build summary: top 2 clusters per diagnosis
summary_data <- cluster_plot_data |>
  dplyr::group_by(diagnosis) |>
  dplyr::slice_max(order_by = cluster_size, n = 3) |>
  dplyr::ungroup() |>
  dplyr::select(diagnosis, cluster, cluster_size, n_ctrl_in_cluster, tissue_bias,
                csf_fraction, lead_CDR3) |>
  dplyr::arrange(diagnosis, desc(cluster_size)) |>
  dplyr::mutate(
    diagnosis = factor(diagnosis, levels = neuropathy_dx),
    ctrl_label = ifelse(n_ctrl_in_cluster == 0, "None", as.character(n_ctrl_in_cluster)),
    tissue_label = case_when(
      tissue_bias == "CSF-enriched" ~ paste0("CSF (", round(csf_fraction * 100), "%)"),
      tissue_bias == "PBMC-enriched" ~ paste0("PBMC (", round((1 - csf_fraction) * 100), "%)"),
      TRUE ~ paste0("Mixed (", round(csf_fraction * 100), "% CSF)")
    )
  )

# Tile plot: one row per cluster, columns = diagnosis properties
p6 <- ggplot(summary_data, aes(y = reorder(paste0(diagnosis, " #", cluster), cluster_size),
                                x = 1)) +
  geom_tile(aes(fill = diagnosis), width = 0.3, alpha = 0.3) +
  geom_text(aes(x = 0.5, label = lead_CDR3), hjust = 0, size = 2.5, fontface = "italic") +
  geom_text(aes(x = 2.5, label = cluster_size), size = 3, fontface = "bold") +
  geom_text(aes(x = 3.2, label = ctrl_label), size = 3) +
  geom_point(aes(x = 3.9, color = tissue_bias), size = 3) +
  scale_fill_manual(values = diagnosis_col, guide = "none") +
  scale_color_manual(values = c("CSF-enriched" = "#E41A1C", "No bias" = "grey50",
                                 "PBMC-enriched" = "#377EB8"), name = "Tissue") +
  scale_x_continuous(breaks = c(0.5, 2.5, 3.2, 3.9),
                     labels = c("Lead CDR3β", "Size", "CTRL\noverlap", "Tissue"),
                     limits = c(0, 4.5)) +
  facet_grid(diagnosis ~ ., scales = "free_y", space = "free_y") +
  theme_pub(base_size = 9) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Top disease-enriched GLIPH2 specificity groups",
    y = "", x = ""
  )

ggsave("results/tcr_comparison/manuscript_figures/fig_gliph_cluster_table.pdf",
       p6, width = 10, height = 12)
message("  Saved fig_gliph_cluster_table.pdf")


# ============================================================================
# Figure 7: GLIPH cluster sharing heatmap
# Which diagnoses share specificity groups?
# ============================================================================
message("Figure 7: Cluster sharing heatmap")

tissue_dx_all <- readxl::read_xlsx("results/tcr_comparison/gliph/tissue_specificity_analysis.xlsx",
                                    sheet = "cluster_tissue_dx")

# Also include CTRL-enriched clusters from main enrichment results
gliph_data <- readxl::read_xlsx("results/tcr_comparison/gliph/disease_enriched_cdr3_motifs.xlsx")

# Build a diagnosis x diagnosis sharing matrix from disease-enriched clusters
all_enriched_clusters <- tissue_dx_all$cluster
enriched_dx <- tissue_dx_all |>
  dplyr::select(cluster, diagnosis) |>
  dplyr::distinct()

# Create co-occurrence: clusters enriched for >1 diagnosis
dx_list <- neuropathy_dx
sharing_mat <- matrix(0, length(dx_list), length(dx_list),
                      dimnames = list(dx_list, dx_list))

for (i in seq_along(dx_list)) {
  cl_i <- enriched_dx$cluster[enriched_dx$diagnosis == dx_list[i]]
  sharing_mat[i, i] <- length(cl_i)  # private clusters
  for (j in seq_along(dx_list)) {
    if (i != j) {
      cl_j <- enriched_dx$cluster[enriched_dx$diagnosis == dx_list[j]]
      sharing_mat[i, j] <- length(intersect(cl_i, cl_j))
    }
  }
}

# Heatmap
pdf("results/tcr_comparison/manuscript_figures/fig_gliph_dx_sharing.pdf", width = 7, height = 6)
pheatmap::pheatmap(
  sharing_mat,
  color = colorRampPalette(c("white", "grey90", "firebrick"))(50),
  display_numbers = TRUE,
  number_format = "%d",
  fontsize_number = 10,
  main = "Shared disease-enriched GLIPH2 clusters between diagnoses",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "grey80"
)
dev.off()


