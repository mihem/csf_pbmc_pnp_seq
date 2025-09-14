##################################################
# compare data with Latorre Nature 2024
# https://www.nature.com/articles/s41586-023-06916-6
##################################################
library(tidyverse)
library(qs)


# latorre supp table 2
# PNS-myelin-reactivi clones
latorre_reactive_clones <-
    readr::read_csv(
        file.path("lookup", "latorre_nature_gbs_supp_table_2.csv")
    ) |>
    dplyr::rename(CTAA = cdr3b_aa) |>
    dplyr::rename_with(function(x) paste0("latorre_", x), .cols = -CTAA)

# latorre public files
latorre_files <- list.files(
    file.path("raw", "latorre"),
    pattern = ".tsv",
    # pattern = "AC",
    full.names = TRUE
)

# read in the files and combine into single data frame
latorre_data <- lapply(latorre_files, readr::read_tsv)
latorre_combined <- bind_rows(latorre_data) |>
    dplyr::select(sample_name, amino_acid) |>
    dplyr::rename(
        CTaa = amino_acid
    ) |>
    drop_na() |>
    dplyr::distinct()

# sanity check
latorre_combined |>
    group_by(sample_name, CTaa) |>
    filter(n() > 1)

qsave(
    latorre_combined,
    file = file.path("objects", "latorre_combined.qs")
)

# read this dataset
sc_tcr <- qread(file.path("objects", "sc_tcr.qs"))

# extract CTaa from TRB data from sc_tcr object
heming_tcr <- sc_tcr@meta.data |>
    dplyr::select(sample, CTaa, tissue_diagnosis) |>
    tidyr::separate_wider_delim(CTaa, names = c("TRA", "TRB"), delim = "_") |>
    dplyr::select(sample, TRB, tissue_diagnosis) |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::rename(CTaa = TRB) |>
    dplyr::rename_with(function(x) paste0("heming_", x), .cols = -CTaa)

# compare with latore supp table 2
latorre_reactive_shared <- heming_tcr |>
    dplyr::inner_join(
        latorre_reactive_clones,
        by = c("CTaa" = "CTAA")
    ) |>
    relocate(CTaa, .before = everything())

qsave(
    latorre_reactive_shared,
    file = file.path("objects", "latorre_reactive_shared.qs")
)

writexl::write_xlsx(
    latorre_reactive_shared,
    path = file.path("results", "tcr", "latorre_reactive_shared.xlsx")
)

# compare with Latorre public data
shared_clones_latorre <- heming_tcr |>
    dplyr::inner_join(latorre_combined, by = c("CTaa"))

# sanity check
heming_tcr |>
    dplyr::filter(CTaa == "CASSRTGPYGYTF")

latorre_combined |>
    dplyr::filter(CTaa == "CASSRTGPYGYTF")

shared_clones_latorre |>
    dplyr::group_by(tissue_diagnosis, sample_name) |>
    dplyr::summarise(n = n()) |>
    dplyr::filter(tissue_diagnosis == "CSF_CTRL") |>
    print(n = Inf)

# frequency of tissue_diagnosis in shared clones
table(shared_clones_latorre$tissue_diagnosis)

plot_shared_clones_latorre <-
    as.data.frame(
        table(shared_clones_latorre$tissue_diagnosis) /
            table(heming_tcr$tissue_diagnosis)
    ) |>
    ggplot(aes(
        x = Var1,
        y = Freq,
        fill = Var1
    )) +
    geom_col() +
    scale_fill_manual(values = sc_tcr@misc$tissue_diagnosis_col) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none"
    ) +
    xlab("")

ggsave(
    plot = plot_shared_clones_latorre,
    filename = file.path("results", "tcr", "shared_clones_latorre.pdf"),
    width = 6,
    height = 4
)
