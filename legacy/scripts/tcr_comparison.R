# ==============================================================================
# TCR Repertoire Comparison with Sukenikova Integration
# ==============================================================================
#
# Consolidated GLIPH2 analysis and figure generation for CSF/PBMC neuropathy
# single-cell TCR data integrated with Sukenikova et al. bulk TCR-beta.
#
# Outputs (4 PDFs):
#   1. fig_gliph_enriched_clusters_overview_bubble.pdf
#   2. fig_gliph_tissue_bias_fraction_by_diagnosis.pdf
#   3. fig_gliph_myelin_reactive_fraction_by_diagnosis.pdf
#   4. sukenikova_reactive_table.pdf
#
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggrepel)
  library(gridExtra)
  library(grid)
  library(turboGliph)
  library(janitor)
})

set.seed(42)

# ==============================================================================
# Constants & helpers
# ==============================================================================

neuropathy_dx <- c("CIAP", "CIDP", "GBS", "MAG", "MFS", "PNC", "CAN", "PPN")
dx_order <- c("CTRL", neuropathy_dx)

# Output directory
fig_dir <- "results/tcr_comparison/manuscript_figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1),
      strip.text       = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "grey95")
    )
}

format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("")
}

safe_save <- function(filename, plot, width, height) {
  out <- file.path(fig_dir, filename)
  ggsave(out, plot = plot, width = width, height = height, device = "pdf")
  message(sprintf("  Saved: %s (%g x %g)", filename, width, height))
}

# ==============================================================================
# Section 1: Load Seurat & extract TCR metadata
# ==============================================================================

message("Loading single-cell TCR data...")
sc_tcr <- qs::qread("data/sc_tcr.qs", nthreads = 6)

# Color palettes
diagnosis_col <- sc_tcr@misc$diagnosis_col
diagnosis_col["GBS_Sukenikova"] <- "#FF7F00"

# Extract TCR metadata from Seurat object
tcr_metadata <- sc_tcr@meta.data |>
  tibble::rownames_to_column("cell_id") |>
  dplyr::filter(!is.na(CTaa)) |>
  dplyr::mutate(
    TRA_CDR3 = gsub("_.*", "", CTaa),
    TRB_CDR3 = gsub(".*_", "", CTaa),
    TRA_length = nchar(TRA_CDR3),
    TRB_length = nchar(TRB_CDR3),
    source = "Heming"
  )

message(sprintf("  Loaded %d cells with TCR, %d patients",
                nrow(tcr_metadata), length(unique(tcr_metadata$patient))))

rm(sc_tcr); invisible(gc())

# ==============================================================================
# Section 2: Sukenikova data import
# ==============================================================================

# --------------------------------------------------------------------------
# 2a. Import bulk immunoseq TSVs
# --------------------------------------------------------------------------
message("  [2a] Importing Adaptive immunoseq TSV files...")

sukenikova_dir <- "data/Sukenikova/TCR"
suk_files <- list.files(sukenikova_dir, pattern = "\\.tsv$", full.names = TRUE)
message(sprintf("    Found %d TSV files", length(suk_files)))

sukenikova_bulk <- lapply(suk_files, function(f) {
  fname <- basename(f)
  message(sprintf("    Reading: %s", fname))

  # Parse metadata from filename: PT##_timepoint_tissue_celltype.tsv
  parts <- str_match(fname, "^(PT\\d+)_(AC|REC)_(blood|CSF|biopsy)_(.+)\\.tsv$")

  df <- readr::read_tsv(f, show_col_types = FALSE) |>
    janitor::clean_names()

  # Filter to productive rearrangements
  if ("frame_type" %in% colnames(df)) {
    df <- df |> dplyr::filter(frame_type == "In")
  } else if ("productive_frequency" %in% colnames(df)) {
    df <- df |> dplyr::filter(!is.na(productive_frequency), productive_frequency > 0)
  }

  # Standardize V gene names: Adaptive TCRBV -> IMGT TRBV
  df <- df |>
    dplyr::mutate(
      v_raw = dplyr::if_else(is.na(v_gene), v_gene_ties, v_gene),
      j_raw = dplyr::if_else(is.na(j_gene), j_gene_ties, j_gene),
      v = v_raw |>
        str_replace_all("TCRBV", "TRBV") |>
        str_extract("TRBV[0-9\\-]+"),
      j = j_raw |>
        str_replace_all("TCRBJ", "TRBJ") |>
        str_extract("TRBJ[0-9\\-]+")
    )

  # CDR3 quality filters
  df <- df |>
    dplyr::filter(
      !is.na(amino_acid),
      amino_acid != "",
      nchar(amino_acid) >= 5,
      nchar(amino_acid) <= 25,
      !str_detect(amino_acid, "\\*|_|X")
    )

  # Deduplicate: group by CDR3 + V gene, sum templates
  count_col <- if ("templates" %in% colnames(df)) "templates" else "seq_reads"

  df_dedup <- df |>
    dplyr::group_by(amino_acid, v) |>
    dplyr::summarize(
      count = sum(.data[[count_col]], na.rm = TRUE),
      j = dplyr::first(na.omit(j)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      frequency = count / sum(count),
      cdr3_length = nchar(amino_acid)
    ) |>
    dplyr::arrange(desc(count))

  # Add metadata
  df_dedup$file_name <- fname
  df_dedup$suk_patient <- if (!is.na(parts[1, 2])) parts[1, 2] else str_extract(fname, "PT\\d+")
  df_dedup$suk_timepoint <- if (!is.na(parts[1, 3])) parts[1, 3] else NA_character_
  df_dedup$suk_tissue <- if (!is.na(parts[1, 4])) parts[1, 4] else NA_character_
  df_dedup$suk_celltype <- if (!is.na(parts[1, 5])) parts[1, 5] else NA_character_
  df_dedup$source <- "Sukenikova"
  df_dedup$diagnosis <- "GBS_Sukenikova"

  df_dedup
}) |> dplyr::bind_rows()

message(sprintf("  Total Sukenikova clonotypes: %d (across %d files)",
                nrow(sukenikova_bulk), length(suk_files)))
message(sprintf("  Unique CDR3b sequences: %d",
                n_distinct(sukenikova_bulk$amino_acid)))

# --------------------------------------------------------------------------
# 2b. Import supplementary table 2 (myelin-reactive clones)
# --------------------------------------------------------------------------
message("  [2b] Importing Sukenikova supp table 2 (myelin-reactive clones)...")

sukenikova_reactive <- readr::read_csv(
  "data/Sukenikova/sukenikova_nature_gbs_supp_table_2.csv",
  show_col_types = FALSE
)

message(sprintf("    %d myelin-reactive clones loaded", nrow(sukenikova_reactive)))
message(sprintf("    Antigens: %s",
                paste(unique(sukenikova_reactive$specificity), collapse = ", ")))

# Create lookup table: CDR3b -> reactivity info
myelin_lookup <- sukenikova_reactive |>
  dplyr::select(cdr3b_aa, trbv, specificity, source, pns_myelin_antigen,
                gliph_cluster_id, public_clonotype, pt) |>
  dplyr::rename(CDR3b = cdr3b_aa, suk_TRBV = trbv)

# --------------------------------------------------------------------------
# 2c. Size normalization
# --------------------------------------------------------------------------
message("  [2c] Size normalization analysis...")

n_heming_trb <- n_distinct(tcr_metadata$TRB_CDR3)
n_suk_total <- n_distinct(sukenikova_bulk$amino_acid)

message(sprintf("    Heming single-cell unique TRB CDR3: %d", n_heming_trb))
message(sprintf("    Sukenikova bulk unique CDR3b: %d", n_suk_total))
message(sprintf("    Ratio: %.1fx", n_suk_total / n_heming_trb))

sukenikova_for_gliph <- sukenikova_bulk |>
  dplyr::group_by(suk_patient, amino_acid, v) |>
  dplyr::summarize(
    total_count = sum(count, na.rm = TRUE),
    j = dplyr::first(na.omit(j)),
    tissues = paste(unique(na.omit(suk_tissue)), collapse = ","),
    timepoints = paste(unique(na.omit(suk_timepoint)), collapse = ","),
    .groups = "drop"
  )

n_suk_dedup <- n_distinct(sukenikova_for_gliph$amino_acid)
message(sprintf("    After per-patient dedup: %d unique CDR3b", n_suk_dedup))

max_per_patient <- 1000
sukenikova_for_gliph <- sukenikova_for_gliph |>
  dplyr::group_by(suk_patient) |>
  dplyr::slice_max(order_by = total_count, n = max_per_patient, with_ties = FALSE) |>
  dplyr::ungroup()

n_suk_downsampled <- n_distinct(sukenikova_for_gliph$amino_acid)
message(sprintf("    After top-%d/patient cap: %d unique CDR3b", max_per_patient,
                n_suk_downsampled))


# ==============================================================================
# Section 3: GLIPH2 analysis
# ==============================================================================

# --------------------------------------------------------------------------
# 3a. Prepare combined GLIPH2 input
# --------------------------------------------------------------------------
message("  [3a] Preparing combined GLIPH2 input...")

heming_trb <- tcr_metadata |>
  dplyr::select(cell_id, TRB_CDR3, sample, patient, tissue, diagnosis) |>
  dplyr::filter(!is.na(TRB_CDR3), TRB_CDR3 != "") |>
  dplyr::rename(CDR3b = TRB_CDR3) |>
  dplyr::mutate(source = "Heming")

heming_vgene <- tcr_metadata |>
  dplyr::mutate(TRBV = gsub("\\..*", "", gsub(".*_", "", CTgene))) |>
  dplyr::select(cell_id, TRBV)

heming_trb <- heming_trb |>
  dplyr::left_join(heming_vgene, by = "cell_id")

suk_trb <- sukenikova_for_gliph |>
  dplyr::transmute(
    cell_id = paste0("suk_", suk_patient, "_", row_number()),
    CDR3b = amino_acid,
    TRBV = v,
    sample = paste0("Suk_", suk_patient),
    patient = suk_patient,
    tissue = tissues,
    diagnosis = "GBS_Sukenikova",
    source = "Sukenikova"
  )

combined_trb <- dplyr::bind_rows(heming_trb, suk_trb)
message(sprintf("    Combined input: %d sequences (%d Heming, %d Sukenikova)",
                nrow(combined_trb), nrow(heming_trb), nrow(suk_trb)))

gliph_input <- data.frame(
  CDR3b = combined_trb$CDR3b,
  stringsAsFactors = FALSE
)

# --------------------------------------------------------------------------
# 3b. Run turboGliph
# --------------------------------------------------------------------------
message("  [3b] Running GLIPH2 via turboGliph...")

gliph_result <- turboGliph::gliph_combined(
  cdr3_sequences = gliph_input,
  local_similarities = TRUE,
  local_method = "fisher",
  global_similarities = TRUE,
  global_method = "fisher",
  motif_length = 3,
  gccutoff = 1,
  vgene_match = "none",
  n_cores = 1,
  result_folder = "",
  lcminove = 10
)

clusters <- gliph_result$cluster_properties
n_clusters <- if (!is.null(clusters)) nrow(clusters) else 0
message(sprintf("    GLIPH2 identified %d specificity groups", n_clusters))

# Parse cluster membership
cluster_membership <- lapply(seq_len(nrow(clusters)), function(i) {
  row <- clusters[i, ]
  members_str <- row$members
  if (is.null(members_str) || is.na(members_str) || members_str == "") return(NULL)

  cdr3_seqs <- strsplit(as.character(members_str), "\\s+")[[1]]
  cdr3_seqs <- cdr3_seqs[cdr3_seqs != ""]
  if (length(cdr3_seqs) == 0) return(NULL)

  data.frame(
    cluster = i,
    CDR3b = cdr3_seqs,
    cluster_tag = if ("tag" %in% names(row)) row$tag else NA_character_,
    cluster_size = if ("cluster_size" %in% names(row)) row$cluster_size else length(cdr3_seqs),
    stringsAsFactors = FALSE
  )
}) |> dplyr::bind_rows()

message(sprintf("    Parsed %d CDR3-to-cluster mappings", nrow(cluster_membership)))

# Map clusters to cell metadata via CDR3 sequence
cluster_meta <- cluster_membership |>
  dplyr::left_join(
    combined_trb |> dplyr::select(cell_id, CDR3b, sample, patient, tissue,
                                   diagnosis, source) |>
      dplyr::distinct(),
    by = "CDR3b"
  )

n_matched <- sum(!is.na(cluster_meta$diagnosis))
message(sprintf("    Matched %d / %d sequences to metadata", n_matched,
                nrow(cluster_meta)))

# --------------------------------------------------------------------------
# 3c. Disease enrichment analysis
# --------------------------------------------------------------------------
message("  [3c] Analyzing diagnosis enrichment in GLIPH clusters...")

total_by_diagnosis <- cluster_meta |>
  dplyr::filter(!is.na(diagnosis)) |>
  dplyr::count(diagnosis, name = "total")

all_dx_to_test <- setdiff(unique(cluster_meta$diagnosis[!is.na(cluster_meta$diagnosis)]), "CTRL")

diagnosis_enrichment <- lapply(
  unique(cluster_meta$cluster[!is.na(cluster_meta$cluster)]),
  function(cl) {
    cl_data <- cluster_meta |> dplyr::filter(cluster == cl, !is.na(diagnosis))
    if (nrow(cl_data) < 3) return(NULL)

    results <- lapply(all_dx_to_test, function(dx) {
      in_cluster_dx <- sum(cl_data$diagnosis == dx)
      in_cluster_ctrl <- sum(cl_data$diagnosis == "CTRL")

      total_dx <- total_by_diagnosis$total[total_by_diagnosis$diagnosis == dx]
      total_ctrl <- total_by_diagnosis$total[total_by_diagnosis$diagnosis == "CTRL"]

      if (length(total_dx) == 0 || length(total_ctrl) == 0 ||
          is.na(total_dx) || is.na(total_ctrl) ||
          total_dx == 0 || total_ctrl == 0) return(NULL)

      out_cluster_dx <- total_dx - in_cluster_dx
      out_cluster_ctrl <- total_ctrl - in_cluster_ctrl

      if (any(c(in_cluster_dx, in_cluster_ctrl, out_cluster_dx, out_cluster_ctrl) < 0))
        return(NULL)

      mat <- matrix(c(in_cluster_dx, in_cluster_ctrl,
                       out_cluster_dx, out_cluster_ctrl), nrow = 2)

      tryCatch({
        test <- fisher.test(mat)
        data.frame(
          cluster = cl, diagnosis = dx,
          n_dx_in_cluster = in_cluster_dx,
          n_ctrl_in_cluster = in_cluster_ctrl,
          cluster_size = nrow(cl_data),
          odds_ratio = test$estimate,
          p.value = test$p.value
        )
      }, error = function(e) NULL)
    }) |> dplyr::bind_rows()

    results
  }
) |> dplyr::bind_rows()

if (nrow(diagnosis_enrichment) > 0) {
  diagnosis_enrichment <- diagnosis_enrichment |>
    dplyr::mutate(
      p.adj = p.adjust(p.value, method = "BH"),
      enrichment_direction = ifelse(odds_ratio > 1, "Disease-enriched", "CTRL-enriched"),
      log2_OR = log2(odds_ratio + 0.01)
    ) |>
    dplyr::arrange(p.adj)

  sig_enrichments <- diagnosis_enrichment |> dplyr::filter(p.adj < 0.1, odds_ratio > 1)
  message(sprintf("    %d disease-enriched clusters (FDR < 0.1)", nrow(sig_enrichments)))

  for (i in seq_len(min(10, nrow(sig_enrichments)))) {
    row <- sig_enrichments[i, ]
    message(sprintf("      Cluster %s - %s: OR=%.2f, FDR=%.4f (n=%d vs %d CTRL)",
                    row$cluster, row$diagnosis, row$odds_ratio, row$p.adj,
                    row$n_dx_in_cluster, row$n_ctrl_in_cluster))
  }
}

# --------------------------------------------------------------------------
# 3d. Tissue compartmentalization
# --------------------------------------------------------------------------
message("  [3d] Analyzing tissue-specificity of GLIPH clusters...")

cluster_meta <- cluster_meta |>
  dplyr::mutate(
    tissue_standard = case_when(
      tissue == "CSF" ~ "CSF",
      tissue %in% c("PBMC", "blood") ~ "PBMC",
      tissue == "biopsy" ~ "Biopsy",
      str_detect(tissue, "blood") ~ "PBMC",
      str_detect(tissue, "CSF") ~ "CSF",
      TRUE ~ tissue
    )
  )

tissue_enrichment <- lapply(
  unique(cluster_meta$cluster[!is.na(cluster_meta$cluster)]),
  function(cl) {
    cl_data <- cluster_meta |> dplyr::filter(cluster == cl, !is.na(tissue_standard))
    if (nrow(cl_data) < 3) return(NULL)

    n_csf <- sum(cl_data$tissue_standard == "CSF")
    n_pbmc <- sum(cl_data$tissue_standard == "PBMC")
    if (n_csf + n_pbmc < 2) return(NULL)

    total_csf <- sum(cluster_meta$tissue_standard == "CSF", na.rm = TRUE)
    total_pbmc <- sum(cluster_meta$tissue_standard == "PBMC", na.rm = TRUE)

    mat <- matrix(c(n_csf, total_csf - n_csf,
                     n_pbmc, total_pbmc - n_pbmc), nrow = 2)
    ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
    if (is.null(ft)) return(NULL)

    data.frame(
      cluster = cl, n_csf = n_csf, n_pbmc = n_pbmc,
      csf_fraction = n_csf / (n_csf + n_pbmc),
      tissue_OR = ft$estimate, tissue_p = ft$p.value
    )
  }
) |> dplyr::bind_rows()

if (nrow(tissue_enrichment) > 0) {
  tissue_enrichment$tissue_padj <- p.adjust(tissue_enrichment$tissue_p, method = "BH")
  tissue_enrichment <- tissue_enrichment |>
    dplyr::mutate(tissue_bias = case_when(
      tissue_padj < 0.1 & tissue_OR > 1 ~ "CSF-enriched",
      tissue_padj < 0.1 & tissue_OR < 1 ~ "PBMC-enriched",
      TRUE ~ "No bias"
    ))

  # Combine tissue and diagnosis enrichment
  cluster_tissue_dx <- diagnosis_enrichment |>
    dplyr::filter(p.adj < 0.1, odds_ratio > 1) |>
    dplyr::left_join(
      tissue_enrichment |> dplyr::select(cluster, csf_fraction, tissue_OR,
                                          tissue_padj, tissue_bias),
      by = "cluster"
    )

  # Tissue bias summary per diagnosis
  tissue_bias_summary <- cluster_tissue_dx |>
    dplyr::count(diagnosis, tissue_bias) |>
    tidyr::pivot_wider(names_from = tissue_bias, values_from = n, values_fill = 0)

  message("    Tissue bias summary:")
  print(tissue_bias_summary)
}


# ==============================================================================
# Section 4: Myelin-reactive clone mapping
# ==============================================================================

message("  [4a] Myelin-reactive clone mapping...")

# Exact CDR3b matching: Heming data <-> supp table 2
heming_for_match <- tcr_metadata |>
  dplyr::select(cell_id, TRB_CDR3, sample, patient, tissue, diagnosis) |>
  dplyr::filter(!is.na(TRB_CDR3)) |>
  dplyr::distinct()

exact_matches <- heming_for_match |>
  dplyr::inner_join(
    myelin_lookup |> dplyr::select(CDR3b, specificity, pns_myelin_antigen,
                                    gliph_cluster_id, public_clonotype, pt),
    by = c("TRB_CDR3" = "CDR3b")
  )

n_exact <- nrow(exact_matches)
message(sprintf("    Exact CDR3b matches: %d", n_exact))

if (n_exact > 0) {
  message("    Matched clones:")
  for (i in seq_len(min(10, n_exact))) {
    row <- exact_matches[i, ]
    message(sprintf("      %s | %s | %s | specificity: %s",
                    row$TRB_CDR3, row$diagnosis, row$tissue, row$specificity))
  }
}

# GLIPH cluster-based matching
myelin_cdr3_set <- unique(myelin_lookup$CDR3b)

myelin_clusters <- cluster_meta |>
  dplyr::filter(CDR3b %in% myelin_cdr3_set) |>
  dplyr::distinct(cluster) |>
  dplyr::pull(cluster)

cluster_matches <- cluster_meta |>
  dplyr::filter(
    cluster %in% myelin_clusters,
    source == "Heming",
    !CDR3b %in% myelin_cdr3_set
  ) |>
  dplyr::select(cluster, CDR3b, patient, tissue, diagnosis) |>
  dplyr::distinct()

n_cluster_matches <- nrow(cluster_matches)
message(sprintf("    GLIPH cluster-based matches: %d Heming CDR3b in myelin-reactive clusters",
                n_cluster_matches))

# Fisher's exact: is any diagnosis enriched for myelin-reactive TCRs vs CTRL?
myelin_enrichment <- lapply(neuropathy_dx, function(dx) {
  dx_total <- sum(tcr_metadata$diagnosis == dx)
  ctrl_total <- sum(tcr_metadata$diagnosis == "CTRL")

  dx_match <- sum(
    exact_matches$diagnosis == dx,
    cluster_matches$diagnosis == dx
  )
  ctrl_match <- sum(
    exact_matches$diagnosis == "CTRL",
    cluster_matches$diagnosis == "CTRL"
  )

  mat <- matrix(c(dx_match, ctrl_match,
                   dx_total - dx_match, ctrl_total - ctrl_match), nrow = 2)
  ft <- tryCatch(fisher.test(mat), error = function(e) list(estimate = NA, p.value = NA))

  data.frame(
    diagnosis = dx,
    n_myelin_reactive = dx_match,
    n_total = dx_total,
    fraction = dx_match / dx_total,
    odds_ratio = ft$estimate,
    p.value = ft$p.value
  )
}) |> dplyr::bind_rows()

myelin_enrichment$p.adj <- p.adjust(myelin_enrichment$p.value, method = "BH")
message("    Myelin reactivity enrichment by diagnosis:")
print(myelin_enrichment |> dplyr::select(diagnosis, n_myelin_reactive, fraction, odds_ratio, p.adj))


# ==============================================================================
# Section 5: Figure generation
# ==============================================================================

message("\n============================================================")
message("Generating manuscript figures")
message("============================================================\n")

# --------------------------------------------------------------------------
# Figure 1: Enriched clusters overview bubble
# --------------------------------------------------------------------------
message("--- GLIPH overview bubble ---")

if (nrow(diagnosis_enrichment) > 0) {

  overview_data <- diagnosis_enrichment |>
    dplyr::filter(p.adj < 0.1, odds_ratio > 1) |>
    dplyr::mutate(
      neg_log10_fdr = -log10(pmax(p.adj, 1e-300)),
      diagnosis = factor(diagnosis, levels = intersect(c(neuropathy_dx, "GBS_Sukenikova"),
                                                        unique(diagnosis)))
    )

  if (nrow(overview_data) > 0) {
    if (nrow(tissue_enrichment) > 0) {
      overview_data <- overview_data |>
        dplyr::left_join(
          tissue_enrichment |> dplyr::select(cluster, tissue_bias),
          by = "cluster"
        ) |>
        dplyr::mutate(tissue_bias = tidyr::replace_na(tissue_bias, "No bias"))
    } else {
      overview_data$tissue_bias <- "No bias"
    }

    tissue_shapes <- c("CSF-enriched" = 24, "No bias" = 21, "PBMC-enriched" = 25)

    top_labels_s1 <- overview_data |>
      dplyr::group_by(diagnosis) |>
      dplyr::slice_max(order_by = cluster_size, n = 2, with_ties = FALSE) |>
      dplyr::ungroup()

    p_s1 <- ggplot(overview_data,
                   aes(x = diagnosis, y = cluster_size)) +
      geom_jitter(aes(fill = diagnosis, size = neg_log10_fdr, shape = tissue_bias),
                  width = 0.2, alpha = 0.75, stroke = 0.4) +
      scale_fill_manual(values = diagnosis_col, guide = "none") +
      scale_shape_manual(values = tissue_shapes, name = "Tissue compartment") +
      scale_size_continuous(range = c(2, 8), name = "-log10(FDR)") +
      geom_text_repel(
        data = top_labels_s1,
        aes(label = paste0("C", cluster)),
        size = 2.5, max.overlaps = 15, fontface = "italic",
        segment.color = "grey50", segment.size = 0.3, seed = 42
      ) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Disease-enriched GLIPH2 specificity groups overview",
        x = "",
        y = "Cluster size (# clonotypes)"
      )

    safe_save("fig_gliph_enriched_clusters_overview_bubble.pdf", p_s1, 9, 6.5)
  }
}

# --------------------------------------------------------------------------
# Figure 2: Tissue bias fraction by diagnosis (absolute counts)
# --------------------------------------------------------------------------
message("--- Tissue bias counts ---")

if (exists("tissue_bias_summary") && nrow(tissue_bias_summary) > 0) {

  tbs_long <- tissue_bias_summary |>
    tidyr::pivot_longer(
      cols = -diagnosis,
      names_to = "tissue_bias",
      values_to = "n_clusters"
    ) |>
    dplyr::mutate(
      tissue_bias = factor(tissue_bias, levels = c("CSF-enriched", "No bias", "PBMC-enriched")),
      diagnosis = factor(diagnosis, levels = intersect(c(dx_order, "GBS_Sukenikova"),
                                                        unique(diagnosis)))
    )

  tb_colors <- c("CSF-enriched" = "#E41A1C", "No bias" = "grey70", "PBMC-enriched" = "#377EB8")

  p_s3 <- ggplot(tbs_long, aes(x = diagnosis, y = n_clusters, fill = tissue_bias)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = tb_colors, name = "Tissue bias") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Tissue compartmentalization of disease-enriched clusters",
      x = "Diagnosis",
      y = "Number of enriched clusters"
    )

  safe_save("fig_gliph_tissue_bias_fraction_by_diagnosis.pdf", p_s3, 8, 5.5)
}

# --------------------------------------------------------------------------
# Figure 3: Myelin-reactive fraction by diagnosis
# --------------------------------------------------------------------------
message("--- Myelin enrichment bar ---")

if (nrow(myelin_enrichment) > 0) {

  me_plot <- myelin_enrichment |>
    dplyr::mutate(
      log2_OR = log2(pmax(odds_ratio, 1e-4)),
      sig_label = sapply(p.adj, format_pval),
      diagnosis = factor(diagnosis, levels = intersect(neuropathy_dx, diagnosis))
    ) |>
    dplyr::filter(!is.na(diagnosis))

  if (nrow(me_plot) > 0) {
    p_s5 <- ggplot(me_plot, aes(x = diagnosis, y = fraction * 100, fill = diagnosis)) +
      geom_col(width = 0.7, show.legend = FALSE) +
      geom_text(aes(label = paste0(sig_label, "\n", sprintf("OR=%.1f", odds_ratio))),
                vjust = -0.3, size = 2.8) +
      scale_fill_manual(values = diagnosis_col) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Myelin-reactive TCR fraction by diagnosis",
        subtitle = "Exact + cluster-level matches to Sukenikova myelin-reactive clones",
        x = "",
        y = "Myelin-reactive TCRs (%)"
      ) +
      expand_limits(y = max(me_plot$fraction * 100, na.rm = TRUE) * 1.4)

    safe_save("fig_gliph_myelin_reactive_fraction_by_diagnosis.pdf", p_s5, 7, 5.5)
  }
}

# --------------------------------------------------------------------------
# Figure 4: Sukenikova reactive table
# --------------------------------------------------------------------------
message("--- Sukenikova reactive table ---")

if (n_exact > 0) {
  exact_table_clean <- exact_matches |>
    dplyr::arrange(tissue, diagnosis, TRB_CDR3) |>
    dplyr::select(TRB_CDR3, patient, tissue, diagnosis, pt, specificity) |>
    dplyr::rename(
      "CDR3 beta" = TRB_CDR3, "Patient" = patient, "Tissue" = tissue,
      "Diagnosis" = diagnosis, "Sukenikova Patient" = pt,
      "Sukenikova Specificity" = specificity
    )

  table_plot <- gridExtra::tableGrob(
    exact_table_clean, rows = NULL,
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

  pdf(file.path(fig_dir, "sukenikova_reactive_table.pdf"),
      width = 12, height = 8)
  grid::grid.draw(table_plot)
  dev.off()
  message("  Saved: sukenikova_reactive_table.pdf (12 x 8)")
}

message("\nDone!")
