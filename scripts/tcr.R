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

sc_rep <- qs::qread("tcells.qs")

sc_rep <- tcells

#sanity check
head(colnames(sc_rep))
head(combined_tcr[[1]]$barcode)

#adjust colnames
#colnames_new <- paste0(colnames(sc_rep), "-1")
#sc_rep <- RenameCells(sc_rep, new.names = colnames_new)

#combine seurat and screpertoire
sc_rep <- scRepertoire::combineExpression(
    combined_tcr,
    sc_rep,
    cloneCall = "aa",
    group.by = "none",
    proportion = FALSE,
    cloneTypes = c(
        Single = 1,
        Small = 5,
        Medium = 20,
        Large = 100,
        Hyperexpanded = 500
    )
)


table(sc_rep$Frequency)
table(sc_rep$cloneType)

sum(!is.na(sc_rep$CTaa)) / dim(sc_rep)[2] #73.5% (21657) t cells with TCR
sum(sapply(combined_tcr, nrow)) - sum(!is.na(sc_rep$CTaa)) # around 2777 TCR lost

#find out most frequent aa
CTaa_sample <- dplyr::count(sc_rep@meta.data, CTaa, sample) |>
    drop_na(CTaa) |>
    pivot_wider(names_from = "sample", values_from = "n") |>
    tibble() |>
    dplyr::mutate(across(everything(), function(x) replace(x, is.na(x), 0))) |>
    dplyr::mutate(sum = rowSums(across(where(is.numeric)))) |> #calculcate rowsum in all numeric
    dplyr::arrange(desc(sum)) |>
    dplyr::select(-sum) |>
    dplyr::select(sort(tidyselect::peek_vars())) |> # sort coolumns
    dplyr::relocate(CTaa) # move CTaa to start

write_csv(CTaa_sample, file.path("results", "repertoire", "CTaa_sample.csv"))

## ## pay attention with frequency, probably based on all contigs from all clusters not only B cell clusters, therefore replace Frequency by counts of the CTaa
## table(sc_rep$Frequency)
## sc_rep$Frequency_old <- sc_rep$Frequency
## sc_rep$Frequency <- NULL
## sc_rep$barcode <- NULL

## sc_rep@meta.data <-
##     sc_rep@meta.data |>
##     tibble::rownames_to_column("barcode") |>
##     dplyr::left_join(CTaa_freq_tcr, by = "CTaa") |>
##     dplyr::rename(Frequency = n) |>
##     tibble::column_to_rownames("barcode")

## CTaa_sample |>
##     pivot_longer(-CTaa)

## sc_rep@meta.data |>
##     tibble::rownames_to_column("barcode") |>
##     dplyr::left_join(select(CTaa_sample, CTaa), by = "CTaa") |>
##     dplyr::rename(Frequency = n) |>
##     tibble::column_to_rownames("barcode")

## sc_rep$cloneType <-
##     factor(
##         cut(sc_rep$Frequency,
##             breaks = c(0,1,5,20,100, 2000),
##             labels = clone_labels),
##         levels = rev(clone_labels))

sc_rep$cloneType <-
    factor(
        sc_rep$cloneType,
        levels = rev(clone_labels)
    )

## table(sc_rep$Frequency)

#plot UMAP with frequency of clonotypes
umap_clone <- DimPlot(
    sc_rep,
    group.by = "cloneType",
    pt.size = .1,
    cols = clone_cols
) +
    theme_rect()
#umap_clone[[1]]$layers[[1]]$aes_params$alpha <- ifelse(is.na(sc_rep@meta.data$cloneType), .01, 0.5)
umap_clone[[1]]$layers[[1]]$aes_params$alpha <-
    case_when(
        is.na(sc_rep$cloneType) ~ 0.01,
        sc_rep$cloneType == "Single (0 < X <= 1)" ~ 0.2,
        sc_rep$cloneType == "Small (1 < X <= 5)" ~ 0.3,
        sc_rep$cloneType == "Medium (5 < X <= 20)" ~ 0.5,
        sc_rep$cloneType == "Large (20 < X <= 100)" ~ 0.7
    )
ggsave(
    file.path("results", "repertoire", "screp_cloneType.pdf"),
    width = 7,
    height = 5
)

#abundance plot screpertoire
stackedPlot(
    object = sc_rep,
    x_axis = "sample",
    y_axis = "cloneType",
    x_order = unique(sc_rep$sample),
    y_order = rev(clone_labels),
    color = clone_cols,
    width = 7,
    height = 3
)


stackedPlot(
    object = sc_rep,
    x_axis = "tissue_level1",
    y_axis = "cloneType",
    x_order = unique(sc_rep$tissue_level1),
    y_order = rev(clone_labels),
    color = clone_cols,
    width = 5,
    height = 5
)

stackedPlot(
    object = sc_rep,
    x_axis = "tissue_level2",
    y_axis = "cloneType",
    x_order = unique(sc_rep$tissue_level2),
    y_order = rev(clone_labels),
    color = clone_cols,
    width = 5,
    height = 5
)


sc_rep_csf <- subset(sc_rep, subset = tissue == "CSF")

stackedPlot(
    object = sc_rep_csf,
    x_axis = "tissue_level2",
    y_axis = "cloneType",
    x_order = unique(sc_rep_csf$tissue_level2),
    y_order = rev(clone_labels),
    color = clone_cols,
    width = 5,
    height = 5
)


sc_rep_pbmc <- subset(sc_rep, subset = tissue == "PBMC")

stackedPlot(
    object = sc_rep_pbmc,
    x_axis = "tissue_level2",
    y_axis = "cloneType",
    x_order = unique(sc_rep_pbmc$tissue_level2),
    y_order = rev(clone_labels),
    color = clone_cols,
    width = 5,
    height = 5
)


#clonal overlay
clonalOverlay(sc_rep, reduction = "umap", freq.cutpoint = 20, bins = 25) +
    scale_color_manual(values = cluster_col_tcells) +
    theme_rect()

ggsave(
    file.path("results", "repertoire", "screp_clonalOverlay.pdf"),
    width = 15,
    height = 10
)


clonalOverlay(
    sc_rep,
    reduction = "umap",
    freq.cutpoint = 20,
    bins = 50,
    facet = "tissue_level2"
) +
    theme_rect() +
    scale_color_manual(values = cluster_col_tcells)

ggsave(
    file.path("results", "repertoire", "screp_clonalOverlay_split.pdf"),
    width = 15,
    height = 10
)

sc_rep$CTaa_top <- sc_rep$CTaa[sc_rep$Frequency > 20] #new column for top expanded clones
#sc_rep$CTaa_top <- factor(sc_rep$CTaa_top, levels = CTaa_freq_bcr$CTaa) # level based on size

alluvialClonotypes(
    sc_rep,
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
    sc_rep,
    cloneCall = "aa",
    y.axes = c("sample"),
    color = "CTaa_top"
)

ggsave(
    file.path("results", "repertoire", "screp_alluvial_sample.pdf"),
    width = 15,
    height = 10
)


## clonalDiversity(sc_rep,
##                 split.by = "tissue_level2",
##                 cloneCall = "aa") +
##      scale_color_manual(values = pals::cols25(25))

clonalOverlap(
    sc_rep,
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
    sc_rep,
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
sc_rep@meta.data |>
    tibble() |>
    dplyr::select(sample, CTaa) |>
    dplyr::filter(sample %in% c("CSF_PNP08", "CSF_PNP10")) |>
    dplyr::distinct() |>
    dplyr::mutate(overlap = duplicated(CTaa)) |>
    dplyr::filter(overlap)
