##################################################
# preprocess scRNAseq data
##################################################

# load libraries ------
library(tidyverse)
library(readxl)
library(Seurat)
library(presto)
library(writexl)
library(scMisc)
library(qs)
library(Polychrome)

# optional libraries for better coding
library(languageserver)
library(httpgd)
library(codegrip)

# renv for better reproducibility
library(renv)

# general settings  ----
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)

# meta data ----
lookup <-
    read_excel(file.path("lookup", "SEED_lookup_v11.xlsx")) |>
    janitor::clean_names() |>
    dplyr::filter(cohort %in% c("scRNA", "scRNA_flow")) |>
    mutate(age = lubridate::time_length(difftime(date, birth_date), "years")) |>
    mutate(
        follow_up = lubridate::time_length(
            difftime(date_follow_up, date),
            "years"
        )
    ) |>
    mutate(
        incat_progress = (incat_follow_up - incat_at_lumbar_puncture) /
            follow_up
    ) |>
    mutate(
        onls_progress = (onls_follow_up - onls_at_lumbar_puncture) /
            follow_up
    ) |>
    mutate(
        mrc_sum_score_60_progress = (mrc_sum_score_60_follow_up -
            mrc_sum_score_60_at_lumbar_puncture) /
            follow_up
    ) |>
    mutate(group = factor(group, levels = c("CTRL", "PNP")), ) |>
    mutate(
        group2 = factor(
            group2,
            levels = c("CTRL", "NIN", "IN", "other")
        )
    ) |>
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
    ) |>
    mutate(across(dml_median_motoric:ncv_sural_sensory, as.numeric))

qsave(lookup, file.path("objects", "lookup.qs"))

# create overview table ---
overview_table <-
    lookup |>
    dplyr::mutate(sex_cat = if_else(sex == "male", 1, 0)) |>
    dplyr::group_by(diagnosis) |>
    dplyr::summarize(
        n = n(),
        age = mean(age, na.rm = TRUE),
        female = (1 - mean(sex_cat, na.rm = TRUE)) * 100
    )

writexl::write_xlsx(
    overview_table,
    file.path("results", "table", "overview_table.xlsx")
)

# find  raw files and match to names from lookup table ---
h5_path <- list.files(
    file.path("raw", "rna"),
    pattern = ".h5",
    recursive = TRUE
)

seq_names <-
    tibble(name = h5_path) |>
    dplyr::mutate(name = str_extract(name, pattern = "(CSF|PBMC|PNP)[^/]+")) |>
    pull(name)

# read in data and create Seurat object -----
sc_multiplex_matrix <-
    lapply(h5_path, scMisc::ReadCellBender_h5) |>
    setNames(seq_names)

sc_multiplex <- lapply(sc_multiplex_matrix, FUN = function(x) {
    CreateSeuratObject(
        counts = x,
        min.cells = 3,
        min.features = 200,
        project = "PNP"
    )
})

# match assignments from vireo ----
donor_path <-
    list.dirs("raw", recursive = TRUE) |>
    stringr::str_subset("vireo")

donor_list <-
    lapply(file.path(donor_path, "donor_ids.tsv"), read_tsv) |>
    setNames(seq_names[!seq_names %in% c("PNP38")])

#rename donor in pool1 (because they were name 1B, 2B, 3B, 4B and not PNP2, PNP3, PNP6 and PNP9)
donor_pool1_lookup <-
    tibble(
        donor_id = c("1B", "2B", "3B", "4B", "doublet", "unassigned"),
        donor_id_new = c(
            "PNP02",
            "PNP03",
            "PNP06",
            "PNP09",
            "doublet",
            "unassigned"
        )
    )

donor_list[["CSF_pool_1"]] <-
    donor_list[["CSF_pool_1"]] |>
    dplyr::left_join(donor_pool1_lookup) |>
    dplyr::select(-donor_id) |>
    dplyr::rename(donor_id = donor_id_new)

donor_list[["PBMC_pool_1"]] <-
    donor_list[["PBMC_pool_1"]] |>
    dplyr::left_join(donor_pool1_lookup) |>
    dplyr::select(-donor_id) |>
    dplyr::rename(donor_id = donor_id_new)

qsave(donor_list, file.path("objects", "donor_list.qs"))

#set donor id from donor list to seurat metadata
for (seq_name in seq_names[!seq_names %in% c("PNP38")]) {
    sc_multiplex[[seq_name]]$pseudonym <-
        tibble(cell = colnames(sc_multiplex[[seq_name]])) |>
        dplyr::left_join(dplyr::select(
            donor_list[[seq_name]],
            cell,
            donor_id
        )) |>
        dplyr::pull(donor_id)
}

# sanity check
dplyr::count(sc_multiplex[[seq_names[[1]]]]@meta.data, pseudonym)

#add doublet and unassigned
pseudonym_lookup <-
    lookup |>
    dplyr::select(pseudonym, patient) |>
    dplyr::add_row(pseudonym = "unassigned", patient = "unassigned") |>
    dplyr::add_row(pseudonym = "doublet", patient = "doublet")

# also add patient ids
# remove PNP52 and PNP62 and PNP56 because less than 100 cells, they do not have a patient id
for (seq_name in seq_names[!seq_names %in% c("PNP38")]) {
    sc_multiplex[[seq_name]]@meta.data <-
        sc_multiplex[[seq_name]]@meta.data |>
        tibble::rownames_to_column("barcode") |>
        dplyr::left_join(pseudonym_lookup) |>
        tibble::column_to_rownames(var = "barcode")
}

# sanity check
dplyr::count(sc_multiplex[[seq_names[[16]]]]@meta.data, patient)
dplyr::count(sc_multiplex[[seq_names[[8]]]]@meta.data, pseudonym)

#split list items based on vireo assignment ----
sc_list <- vector("list")

for (seq_name in seq_names[!seq_names %in% c("PNP38")]) {
    sc_list[[seq_name]] <- SplitObject(
        sc_multiplex[[seq_name]],
        split.by = "patient"
    )
}

#nested list to simple list
sc_list <- unlist(sc_list)

names(sc_list) <-
    names(sc_list) |>
    gsub(
        x = _,
        pattern = "(CSF|PBMC).*(P\\d|doublet|unassigned)",
        replacement = "\\1_\\2"
    )

# remove doublets and unassigned
sc_list <- sc_list[
    !names(sc_list) %in%
        c("CSF_doublet", "CSF_unassigned", "PBMC_doublet", "PBMC_unassigned")
]

# add PNP38 (no multiplexing, so no vireo donor assigment)
sc_list$CSF_P14 <- sc_multiplex$PNP38

# reorder based on new patient ids
sc_list <- sc_list[order(names(sc_list))]

# plot QC: MT and genes ----
for (i in seq_along(sc_list)) {
    sc_list[[i]][["percent_mt"]] <- PercentageFeatureSet(
        sc_list[[i]],
        pattern = "^MT"
    )
}

plot1 <- vector("list", length = length(sc_list))

for (i in seq_along(sc_list)) {
    plot1[[i]] <- FeatureScatter(
        object = sc_list[[i]],
        feature1 = "nCount_RNA",
        feature2 = "percent_mt",
        raster = FALSE
    ) +
        labs(title = names(sc_list)[[i]]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
        NoLegend()
}

plot1_patch <- patchwork::wrap_plots(plot1, ncol = 4)
ggsave(
    file.path("results", "qc", "mt.png"),
    width = 15,
    height = 60,
    limitsize = FALSE
)

plot2 <- vector("list", length = length(sc_list))

for (i in seq_along(sc_list)) {
    plot2[[i]] <- FeatureScatter(
        object = sc_list[[i]],
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA"
    ) +
        labs(title = names(sc_list)[[i]]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
        NoLegend()
}

plot2_patch <- patchwork::wrap_plots(plot2, ncol = 4)
ggsave(
    file.path("results", "qc", "genes.png"),
    width = 15,
    height = 60,
    limitsize = FALSE
)

qs::qsave(sc_list, file.path("objects", "sc_list.qs"))

# filter low quality cells and doublets ---
filter_df <- readr::read_csv(file.path("lookup", "filter_df.csv"))

sc_filter <- vector("list", length(sc_list))

for (i in seq_along(sc_list)) {
    sc_filter[[i]] <-
        subset(
            sc_list[[i]],
            subset = nFeature_RNA > 200 &
                nFeature_RNA < filter_df$rna[[i]] &
                percent_mt < filter_df$mt[[i]]
        )
}

names(sc_filter) <- names(sc_list)

# check filter settings ---

plot4 <- vector("list", length = length(sc_filter))

for (i in seq_along(sc_filter)) {
    plot4[[i]] <- FeatureScatter(
        object = sc_filter[[i]],
        feature1 = "nCount_RNA",
        feature2 = "percent_mt",
        raster = FALSE
    ) +
        labs(title = names(sc_filter)[[i]]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
        NoLegend()
}

plot4_patch <- patchwork::wrap_plots(plot4, ncol = 4)
ggsave(
    file.path("results", "qc", "mt_post.png"),
    width = 15,
    height = 60,
    limitsize = FALSE
)

plot5 <- vector("list", length = length(sc_filter))

for (i in seq_along(sc_filter)) {
    plot5[[i]] <- FeatureScatter(
        object = sc_filter[[i]],
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA"
    ) +
        labs(title = names(sc_filter)[[i]]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) +
        NoLegend()
}

plot5_patch <- patchwork::wrap_plots(plot5, ncol = 4)
ggsave(
    file.path("results", "qc", "genes_post.png"),
    width = 15,
    height = 60,
    limitsize = FALSE
)

# merge seurat objects  ------
sc_merge_pre <- merge(
    x = sc_filter[[1]],
    y = sc_filter[-1],
    merge.data = TRUE,
    add.cell.ids = names(sc_list)
)
sc_merge_pre$sample <- str_extract(colnames(sc_merge_pre), pattern = "[^_]+")

sc_merge_pre$tissue <- str_extract(
    colnames(sc_merge_pre),
    pattern = "(PBMC|CSF)"
)
sc_merge_pre$patient <- str_extract(colnames(sc_merge_pre), pattern = "P\\d+")
sc_merge_pre$sample <- str_extract(
    colnames(sc_merge_pre),
    pattern = "(PBMC|CSF)_P\\d+"
) # extract sample

# this needs to be rejoined and then split again to get the correct names (better naming and required for integration)
sc_merge_pre <- JoinLayers(sc_merge_pre)
sc_merge_pre <- split(x = sc_merge_pre, f = sc_merge_pre$sample)

# add metadata
seq_metadata <-
    lookup |>
    dplyr::select(
        patient,
        sex,
        age,
        group,
        diagnosis,
        incat_at_lumbar_puncture,
        incat_follow_up,
        onls_at_lumbar_puncture,
        onls_follow_up,
        mrc_sum_score_60_at_lumbar_puncture,
        mrc_sum_score_60_follow_up,
        icu
    )

seq_metadata <-
    seq_metadata |>
    dplyr::select(patient, batch)

sc_merge_pre@meta.data <-
    sc_merge_pre@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::left_join(seq_metadata, by = "patient") |>
    tibble::column_to_rownames(var = "barcode")

sc_merge_pre$batch <- paste0(sc_merge_pre$batch, "_", sc_merge_pre$tissue)

str(sc_merge@meta.data)
str(sc_merge_pre@meta.data)

# qc metrics  -----
metrics_files <- list.files(
    file.path("raw", "rna"),
    pattern = "metrics_summary.csv",
    recursive = TRUE,
    full.names = TRUE
)

metrics_data <-
    purrr::map_df(metrics_files, read_csv) |>
    mutate(sample = seq_names, .before = `Estimated Number of Cells`) |>
    arrange(sample)

write_csv(metrics_data, file.path("results", "qc", "cellranger_metrics.csv"))

count_cells <-
    purrr::map_df(sc_list, dim) |>
    dplyr::slice(2) |>
    tidyr::pivot_longer(everything(), names_to = "sample") |>
    dplyr::left_join(dplyr::count(sc_merge_pre@meta.data, sample)) |>
    dplyr::rename(before = value, after = n)

write_csv(count_cells, file.path("results", "qc", "count_cells.csv"))

count_genes <-
    dplyr::bind_cols(
        feature = sc_merge_pre@meta.data$nFeature_RNA,
        sample = sc_merge_pre@meta.data$sample
    ) |>
    dplyr::group_by(sample) |>
    dplyr::summarize(median_genes_after = median(feature))

write_csv(count_genes, file.path("results", "qc", "count_genes.csv"))

qs::qsave(sc_merge_pre, file.path("objects", "sc_merge_pre.qs"))

# normalize ----
sc_merge <- NormalizeData(
    sc_merge_pre,
    verbose = TRUE,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)

sc_merge <-
    sc_merge |>
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
    ScaleData() |>
    RunPCA()


# NMF instead of PCA ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)
library(singlet)

# convert seuratv5 to seuratv3 data, because singlet does not support seuratv5
sc_merge_nmf <- sc_merge
sc_merge_nmf[["RNA3"]] <- as(object = sc_merge_nmf[["RNA"]], Class = "Assay")
DefaultAssay(sc_merge_nmf) <- "RNA3"
sc_merge_nmf[["RNA"]] <- NULL
sc_merge_nmf <- RenameAssays(object = sc_merge_nmf, RNA3 = 'RNA')

# for reproducible NMF models
set.seed(123)
sc_merge_nmf <- RunNMF(sc_merge_nmf, k = 50)

# add NMF embeddings to sc_merge
sc_merge$nmf <- sc_merge_nmf@reductions$nmf
qsave(sc_merge, file.path("objects", "sc_merge.qs"))
