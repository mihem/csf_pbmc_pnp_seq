##################################################
# compare data with Sukenikova Nature 2024
# https://www.nature.com/articles/s41586-023-06916-6
##################################################
library(tidyverse)
library(qs)
library(ggplot2)
library(ggalluvial)

# Sukenikova supp table 2
# PNS-myelin-reactivi clones
sukenikova_reactive_clones <-
    readr::read_csv(
        file.path("lookup", "sukenikova_nature_gbs_supp_table_2.csv")
    ) |>
    dplyr::rename(CTAA = cdr3b_aa) |>
    dplyr::rename_with(function(x) paste0("sukenikova_", x), .cols = -CTAA)

# read this dataset
sc_tcr <- qread(file.path("objects", "sc_tcr.qs"))

# extract CTaa from TRB data from sc_tcr object
heming_tcr <- sc_tcr@meta.data |>
    dplyr::select(sample, CTaa, tissue, diagnosis) |>
    tidyr::separate_wider_delim(CTaa, names = c("TRA", "TRB"), delim = "_") |>
    dplyr::select(sample, TRB, tissue, diagnosis) |>
    tidyr::drop_na() |>
    dplyr::distinct() |>
    dplyr::rename(CTaa = TRB) |>
    dplyr::rename_with(function(x) paste0("heming_", x), .cols = -CTaa)

# compare with latore supp table 2
sukenikova_reactive_shared <- heming_tcr |>
    dplyr::inner_join(
        sukenikova_reactive_clones,
        by = c("CTaa" = "CTAA")
    ) |>
    relocate(CTaa, .before = everything())

qsave(
    sukenikova_reactive_shared,
    file = file.path("objects", "sukenikova_reactive_shared.qs")
)

writexl::write_xlsx(
    sukenikova_reactive_shared,
    path = file.path("results", "tcr", "sukenikova_reactive_shared.xlsx")
)

# alluvial plot
sukenikova_reactive_plot <-
    sukenikova_reactive_shared |>
    ggplot(aes(
        axis1 = CTaa,
        axis2 = heming_tissue,
        axis3 = heming_diagnosis,
        axis4 = sukenikova_specificity
    )) +
    scale_x_discrete(
        limits = c("CTaa", "Tissue",  "Diagnosis", "Specificity"),
        expand = c(.1, .05)
    ) +
    geom_alluvium(aes(fill = sukenikova_specificity)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 90, size = 3) +
    theme_minimal() +
    ylab("Frequency") + 
    theme(legend.position = "none")

ggsave(
    plot = sukenikova_reactive_plot,
    filename = file.path("results", "tcr", "sukenikova_reactive_alluvial.pdf"),
    width = 6,
    height = 6
)

# # Sukenikova public files
# sukenikova_files <- list.files(
#     file.path("raw", "sukenikova"),
#     pattern = ".tsv",
#     # pattern = "AC",
#     full.names = TRUE
# )

# # read in the files and combine into single data frame
# sukenikova_data <- lapply(sukenikova_files, readr::read_tsv)
# sukenikova_combined <- bind_rows(sukenikova_data) |>
#     dplyr::select(sample_name, amino_acid) |>
#     dplyr::rename(
#         CTaa = amino_acid
#     ) |>
#     drop_na() |>
#     dplyr::distinct()

# # sanity check
# sukenikova_combined |>
#     group_by(sample_name, CTaa) |>
#     filter(n() > 1)

# qsave(
#     sukenikova_combined,
#     file = file.path("objects", "sukenikova_combined.qs")
# )

# # compare with sukenikova public data
# shared_clones_sukenikova <- heming_tcr |>
#     dplyr::inner_join(sukenikova_combined, by = c("CTaa"))

# # sanity check
# heming_tcr |>
#     dplyr::filter(CTaa == "CASSRTGPYGYTF")

# sukenikova_combined |>
#     dplyr::filter(CTaa == "CASSRTGPYGYTF")

# shared_clones_sukenikova |>
#     dplyr::group_by(tissue_diagnosis, sample_name) |>
#     dplyr::summarise(n = n()) |>
#     dplyr::filter(tissue_diagnosis == "CSF_CTRL") |>
#     print(n = Inf)

# # frequency of tissue_diagnosis in shared clones
# table(shared_clones_sukenikova$tissue_diagnosis)

# plot_shared_clones_sukenikova <-
#     as.data.frame(
#         table(shared_clones_sukenikova$tissue_diagnosis) /
#             table(heming_tcr$tissue_diagnosis)
#     ) |>
#     ggplot(aes(
#         x = Var1,
#         y = Freq,
#         fill = Var1
#     )) +
#     geom_col() +
#     scale_fill_manual(values = sc_tcr@misc$tissue_diagnosis_col) +
#     theme_classic() +
#     theme(
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#         legend.position = "none"
#     ) +
#     xlab("")

# ggsave(
#     plot = plot_shared_clones_sukenikova,
#     filename = file.path("results", "tcr", "shared_clones_sukenikova.pdf"),
#     width = 6,
#     height = 4
# )
