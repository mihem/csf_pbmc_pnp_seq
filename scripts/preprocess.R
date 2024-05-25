#####################################
# preprocess scRNAseq data
#####################################

# load libraries ------
library(tidyverse)
library(readxl)
library(Seurat)
library(writexl)
library(scMisc)

# optional libraries for better coding
library(languageserver)
library(httpgd)
library(codegrip)

# general settings  ----
future::plan("multicore", workers = 6)
options(future.globals.maxSize = 16000 * 1024^2)

# meta data ----
lookup <-
  read_excel(file.path("lookup", "SEED_lookup_v3.xlsx")) |>
  janitor::clean_names() |>
  mutate(age = lubridate::time_length(difftime(date, birth_date), "years")) |>
  mutate(diagnosis = factor(diagnosis, levels = c("CTRL", "CIAP", "CIDP", "GBS", "MAG", "MFS", "PNC", "CAN", "PPN")))

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

writexl::write_xlsx(overview_table, file.path("results", "table", "overview_table.xlsx"))

# find  raw files and match to names from lookup table ---
h5_path <- list.files(pattern = ".h5", recursive = TRUE)

seq_names <-
  tibble(name = h5_path) |>
  dplyr::mutate(name = str_extract(name, pattern = "(CSF|PBMC|PNP)[^/]+")) |>
  pull(name)

    
# read in data and create Seurat object -----
sc_multiplex_matrix <- lapply(h5_path, scMisc::ReadCellBender_h5) |>
    setNames(seq_names)

sc_multiplex <- lapply(sc_multiplex_matrix,
    FUN = function(x) {
        CreateSeuratObject(counts = x, min.cells = 3, min.features = 200, project = "PNP")
    }
)

scMisc::lss()

# match assignments from vireo ----
donor_path <-
    list.dirs("raw", recursive = TRUE) |>
    stringr::str_subset("vireo")

donor_list <-
    lapply(file.path(donor_path, "donor_ids.tsv"), read_tsv) |>
    setNames(seq_names[!seq_names %in% c("PNP38")])

donor_list$PBMC_pool_1

#rename donor in pool1 (because they were name 1B, 2B, 3B, 4B and not PNP2, PNP3, PNP6 and PNP9)
donor_pool1_lookup <- 
    tibble(donor_id = c("1B", "2B", "3B", "4B", "doublet", "unassigned"),
           donor_id_new = c("PNP02", "PNP03", "PNP06", "PNP09", "doublet", "unassigned"))

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

#set donor id from donor list to seurat metadata
for (seq_name in seq_names[!seq_names %in% c("PNP38")]) {
    sc_multiplex[[seq_name]]$pseudonym <-
        tibble(cell = colnames(sc_multiplex[[seq_name]])) |>
        dplyr::left_join(dplyr::select(donor_list[[seq_name]], cell, donor_id)) |>
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

#also add patient ids
for (seq_name in seq_names[!seq_names %in% c("PNP38")]) {
    sc_multiplex[[seq_name]]@meta.data <- 
        sc_multiplex[[seq_name]]@meta.data |>
        tibble::rownames_to_column("barcode") |>
        dplyr::left_join(pseudonym_lookup) |>
        tibble::column_to_rownames(var = "barcode")
}

dplyr::count(sc_multiplex[[seq_names[[16]]]]@meta.data, patient)
dplyr::count(sc_multiplex[[seq_names[[16]]]]@meta.data, pseudonym)
dplyr::count(sc_multiplex[[seq_names[[16]]]]@meta.data, pseudonym)

#split list items based on vireo assignment
sc_list <- vector("list")

for (seq_name in seq_names[!seq_names %in% c("PNP38")]) {
    sc_list[[seq_name]] <- SplitObject(sc_multiplex[[seq_name]], split.by = "patient")
}

for (seq_name in seq_names[16]) {
    sc_list[[seq_name]] <- SplitObject(sc_multiplex[[seq_name]], split.by = "patient")
}

seq_names

sc_list[[19]]

#nested list to simple list
sc_list <- unlist(sc_list)
sc_list <- sc_list[order(names(sc_list))]

names(sc_list) <- 
    names(sc_list) |>
    str_remove("(_[^_]+){1,}_[^_\\.]+") |>
    str_replace("\\.", "_")


#remove doublets and unassigned
sc_list <- sc_list[!names(sc_list) %in% c("CSF_doublet", "CSF_unassigned", "PBMC_doublet", "PBMC_unassigned")]

sc_list$CSF_PNP11 <- sc_multiplex$PNP38
sc_list <- sc_list[order(names(sc_list))]
