##################################################
# analyze B cell receptor (BCR) data
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


#################################################################################################################
# screpertoire basic analysis
#################################################################################################################
colors_samples <- Polychrome::createPalette(length(combined_bcr), colors_dutch)
names(colors_samples) <- names(combined_bcr)


## #abundance contig
scRepertoire::clonalQuant(combined_bcr, cloneCall = "aa", scale = FALSE) +
  scale_fill_manual(values = colors_samples) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(file.path("results", "bcr", "bcr_quant_abs.pdf"))

## #length contig CD3
scRepertoire::clonalLength(combined_bcr, cloneCall = "aa")
ggsave(file.path("results", "bcr", "bcr_length_aa.pdf"))

#################################################################################################################
# screpertoire analysis with seurat object
#################################################################################################################
sc_bcr <- subset(
  sc_merge,
  cluster %in%
    c(
      "Bmemory",
      "Bnaive",
      "Plasma"
    )
)

#sanity check
sc_bcr_plot <- DimPlot(sc_bcr, reduction = "umap.stacas.ss.all", label = TRUE)
ggsave("sc_bcr_plot.pdf")
head(colnames(sc_bcr))
head(combined_bcr[[1]]$barcode)

# remove empty samples from combined_bcr (CSF_P25)
combined_bcr <- combined_bcr[sapply(combined_bcr, nrow) > 0]

#combine Seurat object with scRepertoire
sc_bcr <- scRepertoire::combineExpression(
  combined_bcr,
  sc_bcr,
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

# save Seurat object with bcr annotations
qs::qsave(
  sc_bcr,
  file = file.path("objects", "sc_bcr.qs")
)

# sanity checks
str(sc_bcr@meta.data)
table(sc_bcr$clonalFrequency)
table(sc_bcr$cloneSize)
# sum(!is.na(sc_bcr$CTaa)) / dim(sc_bcr)[2] T cells with bcr
# sum(sapply(combined_bcr, nrow)) - sum(!is.na(sc_bcr$CTaa)) # 6010 (8.6%) bcr not assigned to a cell

#find out most frequent aa
CTaa_sample <- dplyr::count(sc_bcr@meta.data, CTaa, sample) |>
  drop_na(CTaa) |>
  pivot_wider(names_from = "sample", values_from = "n") |>
  tibble() |>
  dplyr::mutate(across(everything(), function(x) replace(x, is.na(x), 0))) |>
  dplyr::mutate(sum = rowSums(across(where(is.numeric)))) |> #calculcate rowsum in all numeric
  dplyr::arrange(desc(sum)) |>
  dplyr::select(-sum) |>
  dplyr::select(sort(tidyselect::peek_vars())) |> # sort coolumns
  dplyr::relocate(CTaa) # move CTaa to start

write_xlsx(CTaa_sample, file.path("results", "bcr", "CTaa_sample.xlsx"))

# Count number of unique clones (clonotypes) per clone size category
cloneSize_count <- sc_bcr@meta.data |>
  dplyr::filter(!is.na(CTaa)) |> # Only include cells with BCR data
  dplyr::group_by(cloneSize) |>
  dplyr::summarize(
    n_clones = n_distinct(CTaa), # Number of unique clones
    n_cells = n(), # Number of cells
    .groups = "drop"
  ) |>
  dplyr::arrange(cloneSize)

# Optionally save to Excel
write_xlsx(
  cloneSize_count,
  file.path("results", "bcr", "cloneSize_count.xlsx")
)

#plot UMAP with frequency of clonotypes
clone_labels <- levels(sc_bcr$cloneSize)
clone_cols <- setNames(rev(viridis::turbo(length(clone_labels))), clone_labels)

umap_bcr_clone <- DimPlot(
  sc_bcr,
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

umap_bcr_clone[[1]]$layers[[1]]$aes_params$alpha <-
  case_when(
    is.na(sc_bcr$cloneSize) ~ 0.01,
    sc_bcr$cloneSize == "Single (0 < X <= 1)" ~ 0.02,
    sc_bcr$cloneSize == "Small (1 < X <= 5)" ~ 0.2,
    sc_bcr$cloneSize == "Medium (5 < X <= 20)" ~ 0.5,
    sc_bcr$cloneSize == "Large (20 < X <= 100)" ~ 1.0,
    sc_bcr$cloneSize == "Hyperexpanded (X > 100)" ~ 1.0
  )

ggsave(
  file.path("results", "bcr", "umap_cloneSize.pdf"),
  plot = umap_bcr_clone,
  width = 10,
  height = 7
)

#abundance plot screpertoire
stackedPlot(
  object = sc_bcr,
  x_axis = "cluster",
  y_axis = "cloneSize",
  x_order = sc_bcr@misc$cluster_order,
  y_order = clone_labels,
  color = clone_cols,
  width = 3,
  height = 3
)

stackedPlot(
  object = sc_bcr,
  x_axis = "sample",
  y_axis = "cloneSize",
  x_order = unique(sc_bcr$sample),
  y_order = clone_labels,
  color = clone_cols,
  width = 10,
  height = 3
)

stackedPlot(
  object = sc_bcr,
  x_axis = "tissue_group",
  y_axis = "cloneSize",
  x_order = unique(sc_bcr$tissue_group),
  y_order = clone_labels,
  color = clone_cols,
  width = 5,
  height = 5
)

sc_bcr_main_groups <- subset(
  sc_bcr,
  subset = diagnosis %in% c("CTRL", "CIAP", "CIDP", "GBS")
)

sc_bcr_main_groups$tissue_diagnosis <- droplevels(
  sc_bcr_main_groups$tissue_diagnosis,
)

stackedPlot(
  object = sc_bcr_main_groups,
  x_axis = "tissue_diagnosis",
  y_axis = "cloneSize",
  x_order = sc_bcr_main_groups@misc$tissue_diagnosis_order,
  y_order = clone_labels,
  color = clone_cols,
  width = 5,
  height = 5
)

#clonal overlay
bcr_clonal_overlay <- clonalOverlay(
  sc_bcr,
  reduction = "umap.stacas.ss.all",
  cutpoint = 20,
  bins = 25,
  pt.size = 0.1,
  pt.alpha = 0.1
) +
  scale_color_manual(values = sc_bcr@misc$cluster_col) +
  theme_rect() +
  NoLegend() +
  xlab("UMAP1") +
  ylab("UMAP2")

ggsave(
  plot = bcr_clonal_overlay,
  file.path("results", "bcr", "bcr_clonal_overlay.pdf"),
  width = 10,
  height = 10
)

bcr_clonal_overlay_group <-
  clonalOverlay(
    sc_bcr,
    reduction = "umap.stacas.ss.all",
    cutpoint = 20,
    bins = 25,
    pt.size = 0.1,
    pt.alpha = 0.1,
    facet.by = "tissue_group"
  ) +
  scale_color_manual(values = sc_bcr@misc$cluster_col) +
  theme_rect() +
  NoLegend() +
  xlab("UMAP1") +
  ylab("UMAP2")

ggsave(
  plot = bcr_clonal_overlay_group,
  file.path("results", "bcr", "bcr_clonal_overlay_group.pdf"),
  width = 5,
  height = 7
)

bcr_clonal_overlay_diagnosis <-
  clonalOverlay(
    sc_bcr_main_groups,
    reduction = "umap.stacas.ss.all",
    cutpoint = 20,
    bins = 25,
    pt.size = 0.1,
    pt.alpha = 0.1,
  ) +
  scale_color_manual(values = sc_bcr_main_groups@misc$cluster_col) +
  theme_rect() +
  NoLegend() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  facet_wrap(~tissue_diagnosis, nrow = 2)

ggsave(
  plot = bcr_clonal_overlay_diagnosis,
  file.path("results", "bcr", "bcr_clonal_overlay_diagnosis.pdf"),
  width = 10,
  height = 7
)

# alluvial plots main groups
# use 6 here because all other clones are less than 3
bcr_top_clones <- dplyr::count(sc_bcr@meta.data, CTaa) |>
  slice_max(n, n = 6, with_ties = TRUE) |>
  drop_na() |>
  arrange(desc(n))

sc_bcr$CTaa_top <- ifelse(
  sc_bcr$CTaa %in% bcr_top_clones$CTaa,
  sc_bcr$CTaa,
  NA
)

# Create CTaa_top column - keep only clones in bcr_top_clones, set all others to NA
sc_bcr_main_groups$CTaa_top <- ifelse(
  sc_bcr_main_groups$CTaa %in% bcr_top_clones$CTaa,
  sc_bcr_main_groups$CTaa,
  NA
)

# Check the results
table(sc_bcr_main_groups$CTaa_top)

# alluvial plots all groups
bcr_alluvial_tissue_diagnosis <-
  alluvialClones(
    sc_bcr,
    cloneCall = "aa",
    y.axes = c("sample", "patient", "diagnosis", "cluster"),
    color = "CTaa_top"
  ) +
  scale_fill_manual(values = scales::hue_pal()(6))

ggsave(
  plot = bcr_alluvial_tissue_diagnosis,
  file.path(
    "results",
    "bcr",
    "bcr_alluvial_sample_tissue_diagnosis.pdf"
  ),
  width = 15,
  height = 10
)


bcr_overlap_tissue_diagnosis <-
  clonalOverlap(
    sc_bcr_main_groups,
    cloneCall = "aa",
    method = "overlap",
    group.by = "tissue_diagnosis",
  )

ggsave(
  plot = bcr_overlap_tissue_diagnosis,
  file.path(
    "results",
    "bcr",
    "bcr_clonal_overlap_tissue_diagnosis.pdf"
  ),
  width = 15,
  height = 10
)

bcr_overlap_tissue_sample <-
  clonalOverlap(
    sc_bcr_main_groups,
    cloneCall = "aa",
    method = "overlap",
    group.by = "sample",
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave(
  plot = bcr_overlap_tissue_sample,
  file.path(
    "results",
    "bcr",
    "bcr_clonal_overlap_tissue_sample.pdf"
  ),
  width = 20,
  height = 17
)

#find clones that are the same between patients and create a table with samples
shared_clones <- sc_bcr@meta.data |>
  tibble() |>
  dplyr::select(sample, CTaa, tissue_diagnosis, cloneSize, clonalFrequency) |>
  tidyr::drop_na() |>
  dplyr::distinct() |>
  dplyr::group_by(CTaa) |>
  dplyr::filter(n() > 1) |>
  dplyr::arrange(CTaa) |>
  dplyr::ungroup()

# Create a table showing which samples have each clone, sorted by number of occurrences
shared_clones_summary <- shared_clones |>
  dplyr::group_by(CTaa) |>
  dplyr::reframe(
    samples = paste(sample, collapse = ", "),
    tissue_diagnosis = paste(tissue_diagnosis, collapse = ", "),
    sample_count = n(),
    clone_size = unique(cloneSize),
    clonal_frequency = unique(clonalFrequency)
  ) |>
  dplyr::arrange(desc(sample_count), desc(clonal_frequency), CTaa)


# Save the results to an Excel file
write_xlsx(
  shared_clones_summary,
  file.path("results", "bcr", "shared_clones_by_sample_summary.xlsx")
)

# clonal diversity
bcr_clonal_diversity <- clonalDiversity(
  sc_bcr_main_groups,
  cloneCall = "gene",
  group.by = "tissue_diagnosis"
) +
  scale_fill_manual(values = sc_bcr_main_groups@misc$tissue_diagnosis_col)

ggsave(
  plot = bcr_clonal_diversity,
  file.path("results", "bcr", "bcr_clonal_diversity.pdf"),
  width = 10,
  height = 7
)

str(sc_bcr@meta.data$CTgene)

scRepertoire::clonalAbundance(
  combined_bcr,
  cloneCall = "gene",
  scale = FALSE,
  exportTable = TRUE
)

str(sc_bcr@meta.data)

sc_bcr$CT

# abundance of antibody type
ab_type_abundance <- sc_bcr@meta.data |>
  dplyr::select(
    sample,
    patient,
    diagnosis,
    CTgene,
    CTnt,
    CTaa,
    clonalFrequency,
    cloneSize
  ) |>
  filter(str_detect(CTgene, "IGH")) |>
  mutate(IgType = str_extract(CTgene, "IGH[AGMD]")) |>
  filter(!is.na(IgType)) |>
  group_by(IgType, cloneSize) |>
  dplyr::summarize(n = n())

write_xlsx(
  ab_type_abundance,
  file.path("results", "bcr", "ab_type_abundance.xlsx")
)

ab_type_abundance_plot <-
  ggplot(ab_type_abundance, aes(x = cloneSize, y = n, fill = IgType)) +
  geom_col(color = "black", position = "fill") +
  theme_classic() +
  xlab(NULL) +
  ylab("Number of clonotypes") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

ggsave(
  file.path("results", "bcr", "ab_type_abundance.pdf"),
  plot = ab_type_abundance_plot,
  width = 3,
  height = 5,
)
