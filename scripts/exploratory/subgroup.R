# investigate subgroup of CIDP/GBS which has NKCD56bright_1 cells
sc_merge_csf_subgroup <- sc_merge_csf

nkbright_high_group <- c("P16", "P17", "P18", "P21", "P23")
sc_merge_csf_subgroup$subgroup <- ifelse(
  sc_merge_csf_subgroup$patient %in% nkbright_high_group,
  "NKhi",
  as.character(sc_merge_csf_subgroup$diagnosis)
)

sc_merge_csf_subgroup@misc$diagnosis_order <- c(
  "CTRL",
  "NKhi",
  "CIAP",
  "GBS",
  "CIDP",
  "MAG",
  "MFS",
  "PNC",
  "CAN",
  "PPN"
)

scMisc::abBoxPlot(
  object = sc_merge_csf_subgroup,
  cluster_idents = "cluster",
  sample = "sample",
  cluster_order = sc_merge_csf_subgroup@misc$cluster_order,
  group_by = "subgroup",
  group_order = sc_merge_csf_subgroup@misc$diagnosis_order,
  color = sc_merge_csf_subgroup@misc$diagnosis_col,
  number_of_tests = choose(9, 2),
  width = 16
)

lookup_subgroup <- lookup |>
  dplyr::mutate(
    subgroup = dplyr::case_when(
      patient %in% nkbright_high_group ~ "NKhi",
      diagnosis %in% c("CIDP", "GBS") ~ "IN",
      TRUE ~ as.character(diagnosis)
    )
  )
