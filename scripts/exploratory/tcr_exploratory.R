##################################################
# exploratory analysis of TCR repertoire
# requires running annotate.R first
##################################################

sc_tcr <- qs::qread(file.path("objects", "sc_tcr.qs"))
str(sc_tcr@meta.data)

sc_tcr_cd8tem <- subset(
    sc_tcr,
    subset = cluster %in% c("CD8TEM_1", "CD8TEM_2", "CD8TEM_3")
)

clonal_bias <-
    clonalBias(
        sc_tcr_cd8tem,
        cloneCall = "aa",
        split.by = "sample",
        group.by = "cluster",
        n.boots = 10,
        min.expand = 2
    )

ggsave(file.path("results", "tcr", "clonal_bias.pdf"), clonal_bias)

percent_aa <-
    percentAA(
        sc_tcr_main_groups,
        chain = "TRB",
        aa.length = 20,
        group.by = "tissue_diagnosis"
    )

ggsave(
    file.path("results", "tcr", "percent_aa.pdf"),
    percent_aa,
    width = 10,
    height = 20
)

positional_entropy <-
    positionalEntropy(
        sc_tcr_main_groups,
        chain = "TRB",
        aa.length = 20,
        group.by = "tissue_diagnosis"
    )

ggsave(
    file.path("results", "tcr", "positional_entropy.pdf"),
    positional_entropy,
    width = 5,
    height = 5
)

viz_genes <-
    vizGenes(
        sc_tcr,
        x.axis = "TRBV",
        y.axis = NULL,
        plot = "heatmap",
        group.by = "tissue_diagnosis"
    )

ggsave(
    file.path("results", "tcr", "viz_genes.pdf"),
    viz_genes,
    width = 5,
    height = 5
)

df_genes <- percentGenes(
    sc_tcr_main_groups,
    chain = "TRB",
    gene = "Vgene",
    exportTable = TRUE,
    group.by = "tissue_diagnosis"
)

# Performing PCA on the gene usage matrix
pc <- prcomp(df_genes)

# Getting data frame to plot from
df_plot <- as.data.frame(cbind(pc$x[, 1:2], rownames(df_genes)))
colnames(df_plot) <- c("PC1", "PC2", "Sample")
df_plot$PC1 <- as.numeric(df_plot$PC1)
df_plot$PC2 <- as.numeric(df_plot$PC2)

paired_levels <- c(
    "PBMC_CTRL",
    "CSF_CTRL",
    "PBMC_CIDP",
    "CSF_CIDP",
    "PBMC_CIAP",
    "CSF_CIAP",
    "PBMC_GBS",
    "CSF_GBS"
)

paired_colors <- setNames(
    RColorBrewer::brewer.pal(8, "Paired"),
    paired_levels
)

pca_tcr <-
    ggplot(df_plot, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = Sample), shape = 21, size = 5) +
    guides(fill = guide_legend(title = "Samples")) +
    scale_fill_manual(values = paired_colors) +
    # scale_fill_manual(values = sc_tcr_main_groups@misc$tissue_diagnosis_col) +
    theme_bw() +
    labs(title = "PCA of TRBV Gene Usage")

ggsave(
    file.path("results", "tcr", "pca_tcr.pdf"),
    pca_tcr,
    width = 5,
    height = 5
)

# Run clustering, but group calculations by "Patient"
sc_tcr_clonal <- clonalCluster(
    sc_tcr,
    chain = "TRB",
    sequence = "aa",
    threshold = 0.85,
    group.by = "patient"
)

#Define color palette
num_clusters <- length(unique(na.omit(sc_tcr_clonal$TRA_cluster)))
cluster_colors <- hcl.colors(n = num_clusters, palette = "inferno")

umap_clonal <-
    DimPlot(sc_tcr_clonal, group.by = "TRA_cluster", raster = FALSE) +
    scale_color_manual(values = cluster_colors, na.value = "grey") +
    NoLegend()

ggsave(
    file.path("results", "tcr", "umap_trb.pdf"),
    umap_clonal,
    width = 5,
    height = 5
)

# Count the number of cloneSize per sample
pbmc_leftover <- read_xlsx(file.path("lookup", "SEED_lookup_leftover.xlsx"))

samples_expanded <-
    dplyr::count(sc_tcr@meta.data, cloneSize, sample, patient, tissue) |>
    dplyr::filter(tissue == "PBMC") |>
    dplyr::filter(grepl("Hyperexpanded|Large|Medium", cloneSize))

samples_expanded |>
    dplyr::left_join(pbmc_leftover, by = join_by(patient)) |>
    dplyr::select(
        patient,
        sample,
        cloneSize,
        n,
        diagnosis,
        leftover,
        same_date,
        date_sample,
        date_leftover
    ) |>
    dplyr::arrange(patient, desc(cloneSize)) |>
    dplyr::filter(leftover == "yes") |>
    write_xlsx(file.path("results", "tcr", "leftover_pbmc_expanded.xlsx"))

# analyze known self-reactive clones
self_reactive_clones <- c(
    "CVVNRGGGYQKVTF_CAWRGTGAEAFF",
    "CAAKREGSNYKLTF_CASSLAYEQYF"
)

sc_tcr@meta.data |>
    dplyr::filter(CTaa %in% self_reactive_clones) |>
    dplyr::count(sample, cloneSize, cluster, tissue_diagnosis, diagnosis)

str(tcr_contig_list[["PBMC_pool_7_P26"]])

# alpha chain of self reactive clone
p26_self_reactive <- strsplit(self_reactive_clones[2], "_")[[1]]

tcr_contig_list[["PBMC_pool_7_P26"]] |>
    dplyr::filter(cdr3 %in% p26_self_reactive) |>
    write_xlsx(file.path("results", "tcr", "p26_self_reactive_tcr.xlsx"))


##################################################
# Motif similarity analysis for disease clustering
# Using immApex functions to identify disease-specific patterns
##################################################

# Extract kmer data for different groupings ----
# Group by diagnosis to identify disease-specific motifs
# Get kmers for both TRA and TRB chains
kmer_table_trb <- percentKmer(
    sc_tcr,
    chain = "TRB",
    cloneCall = "aa",
    group.by = "sample",
    motif.length = 9,
    min.depth = 1,
    top.motifs = 30,
    exportTable = TRUE
)

# Filter columns (motifs) based on threshold:
# Keep only motifs where at least one sample has proportion > threshold
max_proportion_threshold <- 0.01
motif_colsums <- colSums(kmer_table_trb)
motifs_to_keep <- motif_colsums > max_proportion_threshold
kmer_table_trb <- kmer_table_trb[, motifs_to_keep]

# Hierarchical clustering of samples based on motif similarity ----
# Use tissue_diagnosis grouping for detailed analysis
library(pheatmap)

lookup_kmer <-
    lookup |>
    dplyr::select(patient, diagnosis)

# Create annotation for heatmap
annotation_df <- data.frame(
    sample = rownames(kmer_table_trb)
) |>
    mutate(
        tissue = factor(gsub("_.*", "", sample), levels = c("CSF", "PBMC")),
        patient = gsub(".*_", "", sample)
    ) |>
    left_join(lookup_kmer, by = join_by(patient)) |>
    mutate(diagnosis = factor(diagnosis)) |>
    select(tissue, diagnosis)

# Set rownames to match the matrix
rownames(annotation_df) <- rownames(kmer_table_trb)

# Define annotation colors - only include diagnoses present in the data
diagnosis_levels_present <- levels(annotation_df$diagnosis)
diagnosis_colors_subset <- sc_tcr@misc$diagnosis_col[diagnosis_levels_present]

annotation_colors <- list(
    tissue = c("CSF" = "#E41A1C", "PBMC" = "#377EB8"),
    diagnosis = diagnosis_colors_subset
)

# Create heatmap with clustering
pheatmap(
    t(kmer_table_trb),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    scale = "none",
    main = "TCR Motif Similarity (9-mers, TRB chain)",
    fontsize_row = 5,
    fontsize_col = 5,
    cellheight = 10,
    cellwidth = 10,
    annotation_col = annotation_df,
    annotation_colors = annotation_colors,
    color = viridis::inferno(100),
    filename = file.path("results", "tcr", "kmer_motif_trb_9mers.pdf"),
    border_color = NULL,
    cutree_cols = 6,
    width = 12,
    height = 5
)

# PCA analysis of motif patterns ----
# Perform PCA on the motif matrix (transposed so samples are rows)
# here it makes sense to use smaller motif length to capture more general patterns
kmer_table_tra <- percentKmer(
    sc_tcr,
    chain = "TRB",
    cloneCall = "aa",
    group.by = "sample",
    motif.length = 3,
    min.depth = 3,
    top.motifs = 100,
    exportTable = TRUE
)

pca_motifs <- prcomp(
    kmer_table_tra,
    scale. = TRUE,
    center = TRUE
)

# Create PCA data frame
pca_df <- as.data.frame(pca_motifs$x[, 1:4]) |>
    tibble::rownames_to_column("sample") |>
    dplyr::mutate(
        tissue = gsub("_.*", "", sample),
        patient = gsub(".*_", "", sample)
    ) |>
    dplyr::left_join(lookup_kmer, by = join_by(patient))

# PCA plot PC1 vs PC2
pca_plot_12 <- ggplot(
    pca_df,
    aes(x = PC1, y = PC2, fill = diagnosis, shape = tissue)
) +
    geom_point(size = 5, alpha = 0.8, color = "black") +
    scale_shape_manual(values = c("CSF" = 24, "PBMC" = 21)) +
    scale_fill_manual(values = sc_tcr@misc$diagnosis_col) +
    guides(
        fill = guide_legend(
            override.aes = list(shape = 21, color = "black", size = 5)
        )
    ) +
    labs(
        title = "PCA of TCR Motif Patterns",
        x = "PC1",
        y = "PC2"
    ) +
    theme_bw()

ggsave(
    file.path("results", "tcr", "kmer_pca_plot_pc12.pdf"),
    pca_plot_12,
    width = 8,
    height = 6
)

# Build the edge list: connect sequences with an edit distance of 3 or less

# Extract metadata with TCR sequences as a data frame
tcr_metadata <- sc_tcr@meta.data |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::filter(!is.na(CTaa)) # Keep only cells with TCR data

# Build the edge list using the metadata data frame
edge_list <- immApex::buildNetwork(
    tcr_metadata,
    seq_col = "CTaa",
    threshold = 3
    # dist_type = "hamming",
    # normalize = "maxlen",
    # threshold = 0.1
    # threshold = 2,
    # dist_type = "damerau"
)

# Replace numerical edge list with the sequences and add patient info
edge_list_patients <- data.frame(
    from = tcr_metadata$CTaa[as.numeric(edge_list$from)],
    to = tcr_metadata$CTaa[as.numeric(edge_list$to)],
    dist = edge_list$dist
) |>
    # Add patient info for each clone
    dplyr::left_join(
        tcr_metadata |>
            dplyr::distinct(CTaa, patient) |>
            dplyr::rename(patient_from = patient),
        by = c("from" = "CTaa")
    ) |>
    dplyr::left_join(
        tcr_metadata |>
            dplyr::distinct(CTaa, patient) |>
            dplyr::rename(patient_to = patient),
        by = c("to" = "CTaa")
    ) |>
    # Keep only edges between different patients (similar clones shared)
    dplyr::filter(patient_from != patient_to)

# Count similar clones shared between each pair of patients
patient_similarity <- edge_list_patients |>
    dplyr::count(patient_from, patient_to, name = "n_similar_clones") |>
    # Make symmetric matrix
    dplyr::bind_rows(
        edge_list_patients |>
            dplyr::count(patient_to, patient_from, name = "n_similar_clones") |>
            dplyr::rename(patient_from = patient_to, patient_to = patient_from)
    ) |>
    dplyr::group_by(patient_from, patient_to) |>
    dplyr::summarize(
        n_similar_clones = sum(n_similar_clones),
        .groups = "drop"
    ) |>
    tidyr::pivot_wider(
        names_from = patient_to,
        values_from = n_similar_clones,
        values_fill = 0
    ) |>
    tibble::column_to_rownames("patient_from") |>
    as.matrix()

# Make sure matrix is symmetric and square
all_patients <- sort(unique(c(
    rownames(patient_similarity),
    colnames(patient_similarity)
)))
full_matrix <- matrix(
    0,
    nrow = length(all_patients),
    ncol = length(all_patients),
    dimnames = list(all_patients, all_patients)
)
full_matrix[
    rownames(patient_similarity),
    colnames(patient_similarity)
] <- patient_similarity

# Convert similarity matrix to distance matrix for proper clustering
# Higher similarity should mean lower distance
dist_matrix <- max(full_matrix) - full_matrix
diag(dist_matrix) <- 0 # Distance to self should be 0

# Create annotation for patients with diagnosis
patient_annotation <- tcr_metadata |>
    dplyr::distinct(patient, diagnosis) |>
    dplyr::filter(patient %in% all_patients) |>
    dplyr::arrange(patient) |>
    tibble::column_to_rownames("patient")

# Define annotation colors using existing diagnosis colors
annotation_colors <- list(
    diagnosis = sc_tcr@misc$diagnosis_col
)

# Plot heatmap (display original similarity values, but cluster by distance)
pdf(
    file.path("results", "tcr", "tcr_similar_clones_heatmap.pdf"),
    width = 12,
    height = 10
)
pheatmap::pheatmap(
    full_matrix,
    color = colorRampPalette(c("white", "red"))(100),
    main = "Similar TCR clones shared between patients\n(edit distance ≤ 3)",
    display_numbers = TRUE,
    number_format = "%.0f",
    clustering_method = "ward.D2",
    clustering_distance_rows = as.dist(dist_matrix),
    clustering_distance_cols = as.dist(dist_matrix),
    cutree_cols = 7,
    cutree_rows = 7,
    annotation_row = patient_annotation,
    annotation_col = patient_annotation,
    annotation_colors = annotation_colors,
    border_color = "grey80"
)
dev.off()