##################################################
# analyze T cell receptor (TCR) data
# requires running annotate.R first
##################################################

# load libraries ----
library(qs)
library(Seurat)
library(tidyverse)
library(scRepertoire)

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

str(sc_merge@meta.data)

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

qsave(tcr_contig_list, file = file.path("objects", "tcr_contig_list.qs"))

#str_remove(names(tcr_contig_list), "[^_]+_[^_]_+[^_]_+")
#str_replace(names(tcr_contig_list), pattern = "([^_]+)_\\w+_([^_]+)", replacement = "\\1_\\2")
names(tcr_contig_list) <- gsub(
    x = names(tcr_contig_list),
    pattern = "([^_]+).*_([^_]+)",
    replacement = "\\1_\\2"
)

# manually rename P14 (only CSF available, not multiplexing)
names(tcr_contig_list)[61] <- "CSF_P14"

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
    removeNA = FALSE,
    removeMulti = FALSE
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
            "CD8_NK",
            "NKCD56dim",
            "NKCD56bright"
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
    group.by = "none",
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

# ## ## pay attention with frequency, probably based on all contigs from all clusters not only B cell clusters, therefore replace Frequency by counts of the CTaa
# # sanity check
# table(sc_tcr$clonalFrequency)
# sc_tcr$clonalFrequency_old <- sc_tcr$clonalFrequency
# sc_tcr$clonalFrequency <- NULL

# CTaa_freq_tcr <- dplyr::count(sc_tcr@meta.data, CTaa) |>
#     drop_na() |>
#     arrange(desc(n)) |>
#     tibble()

# sc_tcr@meta.data <-
#     sc_tcr@meta.data |>
#     tibble::rownames_to_column("barcode") |>
#     dplyr::left_join(CTaa_freq_tcr, by = "CTaa") |>
#     dplyr::rename(clonalFrequency = n) |>
#     tibble::column_to_rownames("barcode")

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
    x_axis = "sample",
    y_axis = "cloneSize",
    x_order = unique(sc_tcr$sample),
    y_order = clone_labels,
    color = clone_cols,
    width = 10,
    height = 3
)
str(sc_tcr@meta.data)

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

# sanity chec
table(sc_tcr_main_groups$diagnosis)
sc_tcr_main_groups@misc

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

# Custom clonalOverlay function with smaller point size and custom alpha
customClonalOverlay <- function(
    sc.data,
    reduction = NULL,
    cut.category = "clonalFrequency",
    cutpoint = 30,
    bins = 25,
    pt.size = 0.5,
    pt.alpha = 1,
    facet.by = NULL
) {
    if (!requireNamespace("scRepertoire", quietly = TRUE)) {
        stop("scRepertoire package is required")
    }

    #Forming the data frame to plot
    tmp <- data.frame(
        scRepertoire:::.grabMeta(sc.data),
        scRepertoire:::.get.coord(sc.data, reduction)
    )

    if (cut.category %!in% colnames(tmp)) {
        stop(
            "If filtering the data using a cutpoint, ensure the cut.category correspond to a variable in the meta data."
        )
    }
    #Add facet variable if present
    if (!is.null(facet.by)) {
        facet.by <- tmp[, facet.by]
        tmp <- data.frame(facet.by, tmp)
    }
    #If using cut.category for filtering
    if (!is.null(cut.category) & !is.null(cutpoint)) {
        tmp$include <- ifelse(tmp[, cut.category] >= cutpoint, "Yes", NA)
        tmp2 <- subset(tmp, include == "Yes")
    }

    #Plotting
    plot <- ggplot(
        tmp2,
        mapping = aes(
            x = tmp2[, (ncol(tmp2) - 2)],
            y = tmp2[, (ncol(tmp2) - 1)]
        )
    ) +
        geom_point(
            tmp,
            mapping = aes(
                x = as.numeric(tmp[, (ncol(tmp) - 2)]),
                y = as.numeric(tmp[, (ncol(tmp) - 1)]),
                color = tmp[, "ident"]
            ),
            size = pt.size,
            alpha = pt.alpha
        ) +
        geom_density_2d(color = "black", lwd = 0.25, bins = bins) +
        theme_classic() +
        labs(color = "Active Identity") +
        xlab("Dimension 1") +
        ylab("Dimension 2")
    if (!is.null(facet.by)) {
        plot <- plot + facet_wrap(~facet.by)
    }
    return(plot)
}

#clonal overlay
tcr_clonal_overlay <- customClonalOverlay(
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
    customClonalOverlay(
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
    customClonalOverlay(
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

sc_tcr$CTaa_top <- sc_tcr$CTaa[sc_tcr$Frequency > 20] #new column for top expanded clones
#sc_tcr$CTaa_top <- factor(sc_tcr$CTaa_top, levels = CTaa_freq_bcr$CTaa) # level based on size

alluvialClonotypes(
    sc_tcr,
    cloneCall = "aa",
    y.axes = c("sample", "tissue_level2"),
    color = "CTaa_top"
)

ggsave(
    file.path(
        "results",
        "repertoire",
        "screp_alluvial_sample_tissue_level2.pdf"
    ),
    width = 15,
    height = 10
)


alluvialClonotypes(
    sc_tcr,
    cloneCall = "aa",
    y.axes = c("sample"),
    color = "CTaa_top"
)

ggsave(
    file.path("results", "repertoire", "screp_alluvial_sample.pdf"),
    width = 15,
    height = 10
)


## clonalDiversity(sc_tcr,
##                 split.by = "tissue_level2",
##                 cloneCall = "aa") +
##      scale_color_manual(values = pals::cols25(25))

clonalOverlap(
    sc_tcr,
    cloneCall = "aa",
    split.by = "tissue_level2",
    method = "overlap"
)
ggsave(
    file.path(
        "results",
        "repertoire",
        "screp_clonal_overlap_tissue_level2.pdf"
    ),
    width = 15,
    height = 10
)


clonalOverlap(
    sc_tcr,
    cloneCall = "aa",
    split.by = "sample",
    method = "morisita"
) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(
    file.path("results", "repertoire", "screp_clonal_overlap_sample.pdf"),
    width = 15,
    height = 12
)


#find clones that are the same between patients
sc_tcr@meta.data |>
    tibble() |>
    dplyr::select(sample, CTaa) |>
    dplyr::filter(sample %in% c("CSF_PNP08", "CSF_PNP10")) |>
    dplyr::distinct() |>
    dplyr::mutate(overlap = duplicated(CTaa)) |>
    dplyr::filter(overlap)
