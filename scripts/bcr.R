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
bcr_files <- list.files(
    file.path("raw", "bcr"),
    pattern = "filtered_contig_annotations.csv",
    full.names = TRUE,
    recursive = TRUE
)

bcr_names <- gsub(
    "raw/bcr/(.*)/outs/filtered_contig_annotations.csv",
    "\\1",
    bcr_files
)

bcr_contig_list <- lapply(bcr_files, read.csv)
names(bcr_contig_list) <- bcr_names

# load donor list ----
donor_list <- qs::qread(file.path("objects", "donor_list.qs"))
lookup <- qs::qread(file.path("objects", "lookup.qs"))

#set donor id from donor list to contig list ----
for (sample in bcr_names[!bcr_names %in% c("PNP38")]) {
    bcr_contig_list[[sample]]$pseudonym <-
        tibble(cell = bcr_contig_list[[sample]]$barcode) |>
        dplyr::left_join(dplyr::select(
            donor_list[[sample]],
            cell,
            donor_id
        )) |>
        dplyr::pull(donor_id)
}

#sanity check
dplyr::count(bcr_contig_list[[bcr_names[1]]], pseudonym)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthreads = 6)

#add doublet and unassigned
pseudonym_lookup <- lookup |>
    dplyr::select(pseudonym, patient) |>
    dplyr::add_row(pseudonym = "unassigned", patient = "unassigned") |>
    dplyr::add_row(pseudonym = "doublet", patient = "doublet")

#replace donor id with new patient names
for (sample in bcr_names[!bcr_names %in% c("PNP38")]) {
    bcr_contig_list[[sample]] <-
        bcr_contig_list[[sample]] |>
        dplyr::left_join(pseudonym_lookup) |>
        dplyr::select(-pseudonym)
}

#sanity check
dplyr::count(bcr_contig_list[[bcr_names[1]]], patient)
bcr_contig_list[["PNP38"]]

#assign patient names and remove doublets and unassigned
for (sample in bcr_names[!bcr_names %in% c("PNP38")]) {
    bcr_contig_list[[sample]] <- split(
        bcr_contig_list[[sample]],
        bcr_contig_list[[sample]][["patient"]]
    )
    bcr_contig_list[[sample]]$doublet <- NULL
    bcr_contig_list[[sample]]$unassigned <- NULL
}

#flatten list again
bcr_contig_list <- purrr::list_flatten(bcr_contig_list)

# sanity check
str(bcr_contig_list, max.level = 1)

names(bcr_contig_list) <- gsub(
    x = names(bcr_contig_list),
    pattern = "([^_]+).*_([^_]+)",
    replacement = "\\1_\\2"
)

# # manually rename P14/PNP38 (only CSF available, not multiplexing)
names(bcr_contig_list)[names(bcr_contig_list) == "PNP38"] <- "CSF_P14"

#sort by name alphabetically
bcr_contig_list <- bcr_contig_list[order(names(bcr_contig_list))]

#sanity check
str(bcr_contig_list, max.level = 1)

bcr_contig_sample <- gsub(
    x = names(bcr_contig_list),
    pattern = ".*_([^_]+)",
    replacement = "\\1"
)

bcr_contig_tissue <-
    gsub(
        x = names(bcr_contig_list),
        pattern = "([^_]+)_.*",
        replacement = "\\1"
    )

#combine bcr
combined_bcr <- scRepertoire::combineBCR(
    bcr_contig_list,
    samples = bcr_contig_tissue,
    ID = bcr_contig_sample,
    removeNA = TRUE,
    removeMulti = TRUE
)

qsave(
    combined_bcr,
    file = file.path("objects", "combined_bcr.qs")
)

#sanity check
str(combined_bcr, max.level = 1)


