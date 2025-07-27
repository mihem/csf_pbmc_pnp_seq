##################################################
# analyze T cell receptor (TCR) data
# requires running annotate.R first
##################################################

# load libraries ----
library(qs)
library(Seurat)
library(tidyverse)
library(scRepertoire)
library(scMisc)
library(writexl)

# load filtered contig annotations ---
tcr_files <- list.files(
    file.path("raw", "tcr"),
    pattern = "filtered_contig_annotations.csv",
    full.names = TRUE,
    recursive = TRUE
)

tcr_names <- gsub(
    "raw/tcr/(.*)/outs/filtered_contig_annotations.csv",
    "\\1",
    tcr_files
)

tcr_contig_list <- lapply(tcr_files, read.csv)
names(tcr_contig_list) <- tcr_names

# load donor list ----
donor_list <- qs::qread(file.path("objects", "donor_list.qs"))
lookup <- qs::qread(file.path("objects", "lookup.qs"))

#set donor id from donor list to contig list ----
for (sample in tcr_names[!tcr_names %in% c("PNP38")]) {
    tcr_contig_list[[sample]]$pseudonym <-
        tibble(cell = tcr_contig_list[[sample]]$barcode) |>
        dplyr::left_join(dplyr::select(
            donor_list[[sample]],
            cell,
            donor_id
        )) |>
        dplyr::pull(donor_id)
}

#sanity check
dplyr::count(tcr_contig_list[[tcr_names[1]]], pseudonym)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

#add doublet and unassigned
pseudonym_lookup <- lookup |>
    dplyr::select(pseudonym, patient) |>
    dplyr::add_row(pseudonym = "unassigned", patient = "unassigned") |>
    dplyr::add_row(pseudonym = "doublet", patient = "doublet")

#replace donor id with new patient names
for (sample in tcr_names[!tcr_names %in% c("PNP38")]) {
    tcr_contig_list[[sample]] <-
        tcr_contig_list[[sample]] |>
        dplyr::left_join(pseudonym_lookup) |>
        dplyr::select(-pseudonym)
}

#sanity check
dplyr::count(tcr_contig_list[[tcr_names[1]]], patient)
tcr_contig_list[["PNP38"]]

#assign patient names and remove doublets and unassigned
for (sample in tcr_names[!tcr_names %in% c("PNP38")]) {
    tcr_contig_list[[sample]] <- split(
        tcr_contig_list[[sample]],
        tcr_contig_list[[sample]][["patient"]]
    )
    tcr_contig_list[[sample]]$doublet <- NULL
    tcr_contig_list[[sample]]$unassigned <- NULL
}

#flatten list again
tcr_contig_list <- purrr::list_flatten(tcr_contig_list)

# sanity check
str(tcr_contig_list, max.level = 1)

#str_remove(names(tcr_contig_list), "[^_]+_[^_]_+[^_]_+")
#str_replace(names(tcr_contig_list), pattern = "([^_]+)_\\w+_([^_]+)", replacement = "\\1_\\2")
names(tcr_contig_list) <- gsub(
    x = names(tcr_contig_list),
    pattern = "([^_]+).*_([^_]+)",
    replacement = "\\1_\\2"
)

# manually rename P14/PNP38 (only CSF available, not multiplexing)
names(tcr_contig_list)[names(tcr_contig_list) == "PNP38"] <- "CSF_P14"

#sort by name alphabetically
tcr_contig_list <- tcr_contig_list[order(names(tcr_contig_list))]

#sanity check
str(tcr_contig_list, max.level = 1)

tcr_contig_sample <- gsub(
    x = names(tcr_contig_list),
    pattern = ".*_([^_]+)",
    replacement = "\\1"
)

tcr_contig_tissue <-
    gsub(
        x = names(tcr_contig_list),
        pattern = "([^_]+)_.*",
        replacement = "\\1"
    )

#combine tcr
combined_tcr <- scRepertoire::combineTCR(
    tcr_contig_list,
    samples = tcr_contig_tissue,
    ID = tcr_contig_sample,
    removeNA = TRUE,
    removeMulti = TRUE
)

qsave(
    combined_tcr,
    file = file.path("objects", "combined_tcr.qs")
)

#sanity check
str(combined_tcr, max.level = 1)

#################################################################################################################
# screpertoire basic analysis
#################################################################################################################
colors_samples <- Polychrome::createPalette(length(combined_tcr), colors_dutch)
names(colors_samples) <- names(combined_tcr)


## #abundance contig
scRepertoire::clonalQuant(combined_tcr, cloneCall = "aa", scale = FALSE) +
    scale_fill_manual(values = colors_samples) +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path("results", "tcr", "tcr_quant_abs.pdf"))

## #length contig CD3
scRepertoire::clonalLength(combined_tcr, cloneCall = "aa")
ggsave(file.path("results", "tcr", "tcr_length_aa.pdf"))

#################################################################################################################
# screpertoire analysis with seurat object
#################################################################################################################
sc_tcr <- subset(
    sc_merge,
    cluster %in%
        c(
            "CD4naive_1",
            "CD4TCM_1",
            "CD4TCM_2",
            "CD4TEM",
            "Treg",
            "MAIT",
            "gdT",
            "CD4CTL",
            "CD8naive",
            "CD8TCM",
            "CD8TEM_1",
            "CD8TEM_2",
            "CD8TEM_3",
            "NKCD56bright_1",
            "NKCD56bright_2",
            "NKCD56dim"
        )
)

#sanity check
DimPlot(sc_tcr, reduction = "umap.stacas.ss.all")
head(colnames(sc_tcr))
head(combined_tcr[[1]]$barcode)

#combine Seurat object with scRepertoire
sc_tcr <- scRepertoire::combineExpression(
    combined_tcr,
    sc_tcr,
    cloneCall = "aa",
    group.by = "sample",
    proportion = FALSE,
    cloneSize = c(
        Single = 1,
        Small = 5,
        Medium = 20,
        Large = 100,
        Hyperexpanded = 500
    )
)

# save Seurat object with TCR annotations
qs::qsave(
    sc_tcr,
    file = file.path("objects", "sc_tcr.qs")
)

# sanity checks
str(sc_tcr@meta.data)
table(sc_tcr$clonalFrequency)
table(sc_tcr$cloneSize)
# sum(!is.na(sc_tcr$CTaa)) / dim(sc_tcr)[2] T cells with TCR
# sum(sapply(combined_tcr, nrow)) - sum(!is.na(sc_tcr$CTaa)) # 6010 (8.6%) TCR not assigned to a cell

#find out most frequent aa
CTaa_sample <- dplyr::count(sc_tcr@meta.data, CTaa, sample) |>
    drop_na(CTaa) |>
    pivot_wider(names_from = "sample", values_from = "n") |>
    tibble() |>
    dplyr::mutate(across(everything(), function(x) replace(x, is.na(x), 0))) |>
    dplyr::mutate(sum = rowSums(across(where(is.numeric)))) |> #calculcate rowsum in all numeric
    dplyr::arrange(desc(sum)) |>
    dplyr::select(-sum) |>
    dplyr::select(sort(tidyselect::peek_vars())) |> # sort coolumns
    dplyr::relocate(CTaa) # move CTaa to start

write_xlsx(CTaa_sample, file.path("results", "tcr", "CTaa_sample.xlsx"))

#plot UMAP with frequency of clonotypes
clone_labels <- levels(sc_tcr$cloneSize)
clone_cols <- setNames(rev(viridis::turbo(length(clone_labels))), clone_labels)

umap_tcr_clone <- DimPlot(
    sc_tcr,
    group.by = "cloneSize",
    pt.size = .1,
    reduction = "umap.stacas.ss.all",
    raster = FALSE,
    cols = clone_cols
) +
    theme_rect() +
    xlab("UMAP1") +
    ylab("UMAP2") +
    ggtitle("")

umap_tcr_clone[[1]]$layers[[1]]$aes_params$alpha <-
    case_when(
        is.na(sc_tcr$cloneSize) ~ 0.01,
        sc_tcr$cloneSize == "Single (0 < X <= 1)" ~ 0.02,
        sc_tcr$cloneSize == "Small (1 < X <= 5)" ~ 0.2,
        sc_tcr$cloneSize == "Medium (5 < X <= 20)" ~ 0.3,
        sc_tcr$cloneSize == "Large (20 < X <= 100)" ~ 0.5,
        sc_tcr$cloneSize == "Hyperexpanded (X > 100)" ~ 1.0
    )

ggsave(
    file.path("results", "tcr", "umap_cloneSize.pdf"),
    plot = umap_tcr_clone,
    width = 10,
    height = 7
)

#abundance plot screpertoire
stackedPlot(
    object = sc_tcr,
    x_axis = "cluster",
    y_axis = "cloneSize",
    x_order = sc_tcr@misc$cluster_order,
    y_order = clone_labels,
    color = clone_cols,
    width = 5,
    height = 3
)

stackedPlot(
    object = sc_tcr,
    x_axis = "sample",
    y_axis = "cloneSize",
    x_order = unique(sc_tcr$sample),
    y_order = clone_labels,
    color = clone_cols,
    width = 10,
    height = 3
)

stackedPlot(
    object = sc_tcr,
    x_axis = "tissue_group",
    y_axis = "cloneSize",
    x_order = unique(sc_tcr$tissue_group),
    y_order = clone_labels,
    color = clone_cols,
    width = 5,
    height = 5
)

sc_tcr_main_groups <- subset(
    sc_tcr,
    subset = diagnosis %in% c("CTRL", "CIAP", "CIDP", "GBS")
)

sc_tcr_main_groups$tissue_diagnosis <- droplevels(
    sc_tcr_main_groups$tissue_diagnosis,
)

# sanity chec
table(sc_tcr_main_groups$diagnosis)

stackedPlot(
    object = sc_tcr_main_groups,
    x_axis = "tissue_diagnosis",
    y_axis = "cloneSize",
    x_order = sc_tcr_main_groups@misc$tissue_diagnosis_order,
    y_order = clone_labels,
    color = clone_cols,
    width = 5,
    height = 5
)

#clonal overlay
tcr_clonal_overlay <- clonalOverlay(
    sc_tcr,
    reduction = "umap.stacas.ss.all",
    cutpoint = 20,
    bins = 25,
    pt.size = 0.1,
    pt.alpha = 0.1
) +
    scale_color_manual(values = sc_tcr@misc$cluster_col) +
    theme_rect() +
    NoLegend() +
    xlab("UMAP1") +
    ylab("UMAP2")

ggsave(
    plot = tcr_clonal_overlay,
    file.path("results", "tcr", "tcr_clonal_overlay.pdf"),
    width = 10,
    height = 10
)

tcr_clonal_overlay_group <-
    clonalOverlay(
        sc_tcr,
        reduction = "umap.stacas.ss.all",
        cutpoint = 20,
        bins = 25,
        pt.size = 0.1,
        pt.alpha = 0.1,
        facet.by = "tissue_group"
    ) +
    scale_color_manual(values = sc_tcr@misc$cluster_col) +
    theme_rect() +
    NoLegend() +
    xlab("UMAP1") +
    ylab("UMAP2")

ggsave(
    plot = tcr_clonal_overlay_group,
    file.path("results", "tcr", "tcr_clonal_overlay_group.pdf"),
    width = 5,
    height = 7
)

tcr_clonal_overlay_diagnosis <-
    clonalOverlay(
        sc_tcr_main_groups,
        reduction = "umap.stacas.ss.all",
        cutpoint = 20,
        bins = 25,
        pt.size = 0.1,
        pt.alpha = 0.1,
    ) +
    scale_color_manual(values = sc_tcr_main_groups@misc$cluster_col) +
    theme_rect() +
    NoLegend() +
    xlab("UMAP1") +
    ylab("UMAP2") +
    facet_wrap(~tissue_diagnosis, nrow = 2)

ggsave(
    plot = tcr_clonal_overlay_diagnosis,
    file.path("results", "tcr", "tcr_clonal_overlay_diagnosis.pdf"),
    width = 10,
    height = 7
)



# alluvial plots main groups
tcr_top_clones <- dplyr::count(sc_tcr@meta.data, CTaa) |>
    slice_max(n, n = 10, with_ties = TRUE) |>
    drop_na() |>
    arrange(desc(n))

sc_tcr$CTaa_top <- ifelse(
    sc_tcr$CTaa %in% tcr_top_clones$CTaa,
    sc_tcr$CTaa,
    NA
)

# Check the results
table(sc_tcr_main_groups$CTaa_top)

# Create CTaa_top column - keep only clones in tcr_top_clones, set all others to NA
sc_tcr_main_groups$CTaa_top <- ifelse(
    sc_tcr_main_groups$CTaa %in% tcr_top_clones$CTaa,
    sc_tcr_main_groups$CTaa,
    NA
)

# alluvial plots all groups
tcr_alluvial_tissue_diagnosis <-
    alluvialClones(
        sc_tcr,
        cloneCall = "aa",
        y.axes = c("sample", "patient", "diagnosis", "cluster"),
        color = "CTaa_top"
    ) +
    scale_fill_manual(values = scales::hue_pal()(10))

ggsave(
    plot = tcr_alluvial_tissue_diagnosis,
    file.path(
        "results",
        "tcr",
        "tcr_alluvial_sample_tissue_diagnosis.pdf"
    ),
    width = 15,
    height = 10
)

# alluvial plots main groups
tcr_top_clones_main_groups <- dplyr::count(
    sc_tcr_main_groups@meta.data,
    CTaa
) |>
    slice_max(n, n = 6, with_ties = FALSE) |>
    drop_na() |>
    arrange(desc(n))

sc_tcr_main_groups$CTaa_top <- ifelse(
    sc_tcr_main_groups$CTaa %in% tcr_top_clones_main_groups$CTaa,
    sc_tcr_main_groups$CTaa,
    NA
)

tcr_alluvial_tissue_diagnosis_main_groups <-
    alluvialClones(
        sc_tcr_main_groups,
        cloneCall = "aa",
        y.axes = c("sample", "patient", "diagnosis", "cluster"),
        color = "CTaa_top"
    ) +
    scale_fill_manual(values = scales::hue_pal()(5))

ggsave(
    plot = tcr_alluvial_tissue_diagnosis_main_groups,
    file.path(
        "results",
        "tcr",
        "tcr_alluvial_sample_tissue_diagnosis_main_groups.pdf"
    ),
    width = 15,
    height = 10
)

tcr_overlap_tissue_diagnosis <-
    clonalOverlap(
        sc_tcr_main_groups,
        cloneCall = "aa",
        method = "overlap",
        group.by = "tissue_diagnosis",
    )

ggsave(
    plot = tcr_overlap_tissue_diagnosis,
    file.path(
        "results",
        "tcr",
        "tcr_clonal_overlap_tissue_diagnosis.pdf"
    ),
    width = 15,
    height = 10
)

tcr_overlap_tissue_sample <-
    clonalOverlap(
        sc_tcr_main_groups,
        cloneCall = "aa",
        method = "overlap",
        group.by = "sample",
    ) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

ggsave(
    plot = tcr_overlap_tissue_sample,
    file.path(
        "results",
        "tcr",
        "tcr_clonal_overlap_tissue_sample.pdf"
    ),
    width = 20,
    height = 17
)


#find clones that are the same between patients and create a table with samples
shared_clones <- sc_tcr@meta.data |>
    tibble() |>
    dplyr::select(sample, CTaa, tissue_diagnosis) |> # Add tissue_diagnosis
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::group_by(CTaa) |>
    dplyr::filter(n() > 1) |>
    dplyr::arrange(CTaa) |>
    dplyr::ungroup()

# Create a table showing which samples have each clone, sorted by number of occurrences
shared_clones_summary <- shared_clones |>
    dplyr::group_by(CTaa) |>
    dplyr::summarize(
        samples = paste(sample, collapse = ", "),
        tissue_diagnosis = paste(tissue_diagnosis, collapse = ", "), # Add tissue_diagnosis column
        sample_count = n() # Count how many samples have this clone
    ) |>
    dplyr::arrange(desc(sample_count), CTaa) |> # Sort by count (desc) then by CTaa
    dplyr::select(CTaa, samples, tissue_diagnosis, sample_count) # Keep the count column for reference

# Save the results to an Excel file
write_xlsx(
    shared_clones_summary,
    file.path("results", "tcr", "shared_clones_by_sample_summary.xlsx")
)

#
tcr_clonal_distribution <- clonalSizeDistribution(
    sc_tcr,
    cloneCall = "aa",
    method = "ward.D2",
    group.by = "sample"
)

ggsave(
    plot = tcr_clonal_distribution,
    file.path("results", "tcr", "tcr_clonal_distribution.pdf"),
    width = 10,
    height = 7
)

tcr_clonal_diversity <- clonalDiversity(
    sc_tcr_main_groups,
    cloneCall = "gene",
    group.by = "tissue_diagnosis"
) +
    scale_fill_manual(values = sc_tcr_main_groups@misc$tissue_diagnosis_col)

ggsave(
    plot = tcr_clonal_diversity,
    file.path("results", "tcr", "tcr_clonal_diversity.pdf"),
    width = 10,
    height = 7
)

Idents(sc_tcr_main_groups) <- sc_tcr_main_groups$cluster

tcr_startrac_diversity <- StartracDiversity(
    sc_tcr_main_groups,
    cloneCall = "strict",
    chain = "both",
    type = "tissue",
    group.by = "patient",
) +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

ggsave(
    plot = tcr_startrac_diversity,
    file.path("results", "tcr", "tcr_startrac_diversity.pdf"),
    width = 4,
    height = 5
)

# analyze MAIT cells
sc_tcr <- annotateInvariant(
    sc_tcr,
    type = "MAIT",
    species = "human"
)

# annotate iNKT cells
sc_tcr <- annotateInvariant(
    sc_tcr,
    type = "iNKT",
    species = "human"
)

# Function to create invariant cell plots with highlighting
create_invariant_plot <- function(
    seurat_obj,
    cell_type,
    width = 5,
    height = 5
) {
    score_column <- paste0(cell_type, ".score")
    title <- paste(cell_type, "Cells Highlighted")
    filename <- paste0("highlighted_", cell_type, "_cells.pdf")

    # Create the plot
    plot <- DimPlot(
        seurat_obj,
        reduction = "umap.stacas.ss.all",
        group.by = score_column,
        pt.size = 0.1,
        cols = c("0" = "grey", "1" = "red"),
        raster = FALSE
    ) +
        NoLegend() +
        theme_rect() +
        xlab("UMAP1") +
        ylab("UMAP2") +
        ggtitle(title)

    # Modify the plot to adjust transparency
    plot[[1]]$layers[[1]]$aes_params$alpha <- ifelse(
        seurat_obj@meta.data[[score_column]] == 1,
        1.0, # Full opacity for positive cells (score=1)
        0.01 # Low opacity for negative cells (score=0)
    )

    # Save the plot
    ggsave(
        plot = plot,
        file.path("results", "tcr", filename),
        width = width,
        height = height
    )

    return(plot)
}

# Create plots for both MAIT and iNKT cells
custom_MAIT_plot <- create_invariant_plot(sc_tcr, "MAIT")
custom_iNKT_plot <- create_invariant_plot(sc_tcr, "iNKT")

# Calculate mean clonal frequency by tissue type (CSF vs PBMC) for each patient (summarized)
correlation_data_mean <-
    sc_tcr@meta.data |>
    group_by(patient, tissue, diagnosis) |>
    summarize(
        clonalFrequency = mean(clonalFrequency, na.rm = TRUE),
        .groups = "drop"
    ) |>
    pivot_wider(
        names_from = tissue,
        values_from = clonalFrequency,
        names_prefix = "clonalFreq_"
    ) |>
    drop_na(clonalFreq_CSF, clonalFreq_PBMC)

# Plot CSF vs PBMC clonal frequency (mean values)
tissue_correlation_plot_mean <- correlation_data_mean |>
    ggplot(aes(x = clonalFreq_CSF, y = clonalFreq_PBMC, color = diagnosis)) +
    geom_point(size = 3) +
    scale_color_manual(values = sc_tcr@misc$diagnosis_col) +
    labs(
        x = "Mean Clonal Frequency (CSF)",
        y = "Mean Clonal Frequency (PBMC)",
        title = "Mean Clonal Frequency per Patient"
    ) +
    theme_classic()

# Save the mean plot
ggsave(
    plot = tissue_correlation_plot_mean,
    file.path("results", "tcr", "tcr_csf_pbmc_correlation_mean.pdf"),
    width = 8,
    height = 6
)

# Use all individual cell values - match cells from same clones across tissues
correlation_data_all <-
    sc_tcr@meta.data |>
    filter(!is.na(clonalFrequency) & !is.na(CTaa)) |>
    select(patient, tissue, diagnosis, clonalFrequency, CTaa) |>
    # Split into CSF and PBMC data
    group_split(tissue) |>
    # Name the list elements
    setNames(c("CSF", "PBMC")) |>
    # Join CSF and PBMC data by patient and clone (CTaa)
    reduce(inner_join, by = c("patient", "diagnosis", "CTaa")) |>
    # Rename columns to be clear
    rename(
        clonalFreq_CSF = clonalFrequency.x,
        clonalFreq_PBMC = clonalFrequency.y
    ) |>
    distinct()


# Plot CSF vs PBMC clonal frequency (matched clones)
tissue_correlation_plot_all <- correlation_data_all |>
    ggplot(aes(
        x = clonalFreq_CSF,
        y = clonalFreq_PBMC,
        color = diagnosis
    )) +
    geom_jitter(size = 1, alpha = 0.6) +
    scale_color_manual(values = sc_tcr@misc$diagnosis_col) +
    labs(
        x = "Clonal Frequency (CSF)",
        y = "Clonal Frequency (PBMC)",
        title = "Matched Clones Across Tissues"
    ) +
    theme_classic()

# Save the all values plot
ggsave(
    plot = tissue_correlation_plot_all,
    file.path("results", "tcr", "tcr_csf_pbmc_correlation_all.pdf"),
    width = 8,
    height = 6
)

# correlate TCR clonality with CSF protein
sc_tcr_csf <- subset(sc_tcr, tissue %in% c("CSF"))

csf_protein_lookup <-
    lookup |>
    select(pseudonym, csf_protein)

sc_tcr_csf@meta.data <-
    sc_tcr_csf@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(csf_protein_lookup) |>
    tibble::column_to_rownames("barcode")

# Calculate correlation between clonality and CSF protein
protein_correlation_data <-
    sc_tcr_csf@meta.data |>
    group_by(sample) |>
    summarize(
        clonalFrequency = mean(clonalFrequency, na.rm = TRUE),
        csf_protein = mean(csf_protein, na.rm = TRUE),
        diagnosis = unique(diagnosis)
    )

# Calculate correlation statistics
protein_correlation_test <- cor.test(
    protein_correlation_data$csf_protein,
    protein_correlation_data$clonalFrequency,
    method = "pearson"
)

# Extract statistics
r_value <- protein_correlation_test$estimate
p_value <- protein_correlation_test$p.value

# Format statistics for display
protein_stats_text <- paste0(
    "r = ",
    signif(r_value, 3),
    ", p = ",
    signif(p_value, 3)
)

protein_correlation_plot <-
    protein_correlation_data |>
    ggplot(aes(x = csf_protein, y = clonalFrequency, color = diagnosis)) +
    geom_point(size = 3) +
    scale_color_manual(values = sc_tcr_csf@misc$diagnosis_col) +
    geom_smooth(method = "lm", color = "blue") +
    labs(
        x = "CSF Protein (mg/L)",
        y = "Mean Clonal Frequency",
        subtitle = protein_stats_text
    ) +
    theme_classic()

# Save the correlation plot
ggsave(
    plot = protein_correlation_plot,
    file.path("results", "tcr", "tcr_clonality_csf_protein_correlation.pdf"),
    width = 8,
    height = 6
)

sc_tcr <- qread(file.path("objects", "sc_tcr.qs"))
str(sc_tcr@meta.data)

clonal_bias <-
    clonalBias(
        sc_tcr,
        cloneCall = "aa",
        split.by = "patient",
        group.by = "cluster",
        n.boots = 10,
        min.expand = 5
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
