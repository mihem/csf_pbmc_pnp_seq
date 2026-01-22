##################################################
# compare data with Sukenikova Nature 2024
# https://www.nature.com/articles/s41586-023-06916-6
##################################################
library(tidyverse)
library(qs)
library(ggplot2)
library(gridExtra)
library(grid)

# Sukenikova supp table 2
# PNS-myelin-reactivi clones
sukenikova_reactive_clones <-
  readr::read_csv(
    file.path("lookup", "sukenikova_nature_gbs_supp_table_2.csv")
  ) |>
  dplyr::rename(CTaa_TRB = cdr3b_aa) |>
  dplyr::mutate(diagnosis = "GBS") |>
  dplyr::rename_with(function(x) paste0("sukenikova_", x), .cols = -CTaa_TRB)

# read this dataset
sc_tcr <- qread(file.path("objects", "sc_tcr.qs"))

# extract CTaa from TRB data from sc_tcr object
heming_tcr <- sc_tcr@meta.data |>
  dplyr::select(patient, CTaa, CTgene, CTnt, tissue, diagnosis) |>
  tidyr::separate_wider_delim(CTaa, names = c("TRA", "TRB"), delim = "_") |>
  tidyr::drop_na() |>
  dplyr::distinct() |>
  dplyr::rename(CTaa_TRB = TRB) |>
  dplyr::rename(CTaa_TRA = TRA) |>
  dplyr::rename_with(function(x) paste0("heming_", x), .cols = -CTaa_TRB)

# compare with latore supp table 2
sukenikova_reactive_shared <- heming_tcr |>
  dplyr::inner_join(
    sukenikova_reactive_clones,
    by = c("CTaa_TRB" = "CTaa_TRB")
  ) |>
  relocate(CTaa_TRB, .before = everything())

data.frame(sukenikova_reactive_shared)

qsave(
  sukenikova_reactive_shared,
  file = file.path("objects", "sukenikova_reactive_shared.qs")
)

writexl::write_xlsx(
  sukenikova_reactive_shared,
  path = file.path("results", "tcr", "sukenikova_reactive_shared.xlsx")
)


# Prepare clean data for the table
sukenikova_table_clean <- sukenikova_reactive_shared |>
  dplyr::arrange(heming_tissue, heming_diagnosis, CTaa_TRB) |>
  dplyr::select(
    CTaa_TRB,
    heming_patient,
    heming_tissue,
    heming_diagnosis,
    sukenikova_pt,
    sukenikova_diagnosis,
    sukenikova_source,
    sukenikova_specificity,
    sukenikova_clone_id
  ) |>
  dplyr::rename(
    "CDR3 beta" = CTaa_TRB,
    "Patient" = heming_patient,
    "Tissue" = heming_tissue,
    "Diagnosis" = heming_diagnosis,
    "Sukenikova Patient" = sukenikova_pt,
    "Sukenikova Diagnosis" = sukenikova_diagnosis,
    "Sukenikova Tissue" = sukenikova_source,
    "Sukenikova Specificity" = sukenikova_specificity
  )

table_plot <- gridExtra::tableGrob(
  sukenikova_table_clean,
  rows = NULL,
  theme = gridExtra::ttheme_default(
    core = list(
      fg_params = list(cex = 0.8),
      bg_params = list(fill = c("white", "lightgray"), alpha = 0.5)
    ),
    colhead = list(
      fg_params = list(cex = 0.9, fontface = "bold"),
      bg_params = list(fill = "darkgray", alpha = 0.8)
    )
  )
)

# Save as PDF using base R graphics (very reliable)
pdf(
  file.path("results", "table", "sukenikova_reactive_table_grid.pdf"),
  width = 12,
  height = 8
)
grid::grid.draw(table_plot)
dev.off()

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
