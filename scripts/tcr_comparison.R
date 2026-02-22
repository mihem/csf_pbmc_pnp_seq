# ==============================================================================
# TCR Repertoire Comparison with Sukenikova Integration
# ==============================================================================
#
# Consolidated analysis script for CSF/PBMC neuropathy TCR sequencing.
# Compares TCR repertoires across neuropathy diagnoses (GBS, CIDP, CIAP, MAG,
# MFS, PNC, CAN, PPN) vs healthy controls (CTRL).
#
# Sections:
#   0. Setup & Data Loading
#   1. Sukenikova Data Import (Adaptive immunoseq bulk TCRB)
#   2. V/J Gene Usage Analysis
#   3. GLIPH2 Combined Analysis (single-cell + Sukenikova)
#   4. Sukenikova Cross-Reference (myelin reactivity, cluster mirroring)
#   5. Save Cache for Figure Script
#
# ==============================================================================

# ==============================================================================
# Section 0: Setup & Data Loading
# ==============================================================================

if (nzchar(Sys.getenv("RENV_PROJECT"))) {
  Sys.setenv(RENV_PROJECT = "")
  message("Note: renv deactivated, using system packages")
}

set.seed(42)

library(qs)
library(Seurat)
library(tidyverse)
library(scRepertoire)
library(writexl)
library(readxl)
library(ggpubr)
library(ggsignif)
library(viridis)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(immApex)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(emmeans)
library(ggrepel)
library(gridExtra)
library(grid)
library(turboGliph)

# Constants
neuropathy_dx <- c("CIAP", "CIDP", "GBS", "MAG", "MFS", "PNC", "CAN", "PPN")
dx_order <- c("CTRL", neuropathy_dx)
aa_codes <- c("A","C","D","E","F","G","H","I","K","L",
              "M","N","P","Q","R","S","T","V","W","Y")

# Output directories
output_dirs <- c(
  "results/tcr_comparison",
  "results/tcr_comparison/gene_usage",
  "results/tcr_comparison/gliph",
  "results/tcr_comparison/sukenikova",
  "results/tcr_comparison/manuscript_figures",
  "results/tcr_comparison/tables",
  "data/cache"
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Load single-cell TCR data
message("Loading single-cell TCR data...")
sc_tcr <- qs::qread("data/sc_tcr.qs", nthreads = 6)

# Color palettes
diagnosis_col <- sc_tcr@misc$diagnosis_col
tissue_diagnosis_col <- sc_tcr@misc$tissue_diagnosis_col
cluster_col <- sc_tcr@misc$cluster_col
tissue_col <- c("CSF" = "#E41A1C", "PBMC" = "#377EB8")

# Add GBS_Sukenikova color
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

# Unique clonotypes per sample
tcr_clonotypes <- tcr_metadata |>
  dplyr::group_by(CTaa, sample, patient, tissue, group, diagnosis,
                  tissue_group, tissue_diagnosis) |>
  dplyr::summarize(
    n_cells = n(),
    TRA_CDR3 = first(TRA_CDR3),
    TRB_CDR3 = first(TRB_CDR3),
    TRA_length = first(TRA_length),
    TRB_length = first(TRB_length),
    .groups = "drop"
  )

sample_summary <- tcr_metadata |>
  dplyr::distinct(sample, patient, tissue, group, diagnosis)

message(sprintf("  Loaded %d cells with TCR, %d unique clonotypes, %d patients",
                nrow(tcr_metadata), nrow(tcr_clonotypes),
                length(unique(tcr_metadata$patient))))

# ==============================================================================
# Helper functions
# ==============================================================================

# Mixed model: diagnosis vs CTRL with patient random effect
run_diagnosis_stats <- function(data, value_col, include_tissue = TRUE) {
  data_clean <- data |>
    dplyr::filter(!is.na(.data[[value_col]]), !is.na(diagnosis), !is.na(patient)) |>
    dplyr::mutate(diagnosis = factor(diagnosis, levels = c("CTRL", neuropathy_dx)))

  if (nrow(data_clean) < 20 || length(unique(data_clean$patient)) < 5) {
    return(data.frame(
      diagnosis = neuropathy_dx, estimate = NA, SE = NA, t.ratio = NA,
      p.value = NA, p.adj = NA, method = "insufficient_data"
    ))
  }

  has_tissues <- length(unique(data_clean$tissue)) == 2

  tryCatch({
    if (include_tissue && has_tissues) {
      formula_str <- paste0(value_col, " ~ diagnosis * tissue + (1|patient)")
    } else {
      formula_str <- paste0(value_col, " ~ diagnosis + (1|patient)")
    }

    model <- lmerTest::lmer(as.formula(formula_str), data = data_clean)
    emm <- emmeans::emmeans(model, "diagnosis")
    contrasts <- emmeans::contrast(emm, method = "trt.vs.ctrl", ref = 1)

    result <- as.data.frame(contrasts) |>
      dplyr::mutate(
        diagnosis = gsub(" - CTRL", "", contrast),
        method = "emmeans_dunnett"
      ) |>
      dplyr::select(diagnosis, estimate, SE, t.ratio, p.value, method)

    means <- data_clean |>
      dplyr::group_by(diagnosis) |>
      dplyr::summarize(mean = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop")

    ctrl_mean <- means$mean[means$diagnosis == "CTRL"]
    result <- result |>
      dplyr::left_join(means, by = "diagnosis") |>
      dplyr::rename(mean_dx = mean) |>
      dplyr::mutate(mean_ctrl = ctrl_mean)

    result

  }, error = function(e) {
    ctrl_data <- data_clean |> dplyr::filter(diagnosis == "CTRL")

    results <- lapply(neuropathy_dx, function(dx) {
      dx_data <- data_clean |> dplyr::filter(diagnosis == dx)
      if (nrow(dx_data) < 3) {
        return(data.frame(diagnosis = dx, estimate = NA, SE = NA, t.ratio = NA,
                          p.value = NA, mean_dx = NA, mean_ctrl = NA, method = "insufficient"))
      }

      test <- tryCatch(
        wilcox.test(ctrl_data[[value_col]], dx_data[[value_col]]),
        error = function(e) list(p.value = NA)
      )

      data.frame(
        diagnosis = dx,
        estimate = mean(dx_data[[value_col]], na.rm = TRUE) - mean(ctrl_data[[value_col]], na.rm = TRUE),
        SE = NA, t.ratio = NA,
        p.value = test$p.value,
        mean_dx = mean(dx_data[[value_col]], na.rm = TRUE),
        mean_ctrl = mean(ctrl_data[[value_col]], na.rm = TRUE),
        method = "wilcoxon"
      )
    }) |> dplyr::bind_rows()

    results
  })
}

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
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

# Extract k-mer motifs from CDR3 sequences
extract_cdr3_motifs <- function(cdr3_seqs, k = 3) {
  motifs <- unlist(lapply(cdr3_seqs, function(seq) {
    if (is.na(seq) || nchar(seq) < k) return(character(0))
    sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
  }))
  table(motifs)
}


# ==============================================================================
# Section 1: Sukenikova Data Import
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Section 1: Sukenikova Data Import")
message(paste(rep("=", 60), collapse = ""))

# --------------------------------------------------------------------------
# 1a. Import bulk immunoseq TSVs
# --------------------------------------------------------------------------
message("  [1a] Importing Adaptive immunoseq TSV files...")

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
      # Use v_gene_name, fall back to v_gene_name_ties
      v_raw = dplyr::if_else(is.na(v_gene), v_gene_ties, v_gene),
      j_raw = dplyr::if_else(is.na(j_gene), j_gene_ties, j_gene),
      # Standardize V gene
      v = v_raw |>
        str_replace_all("TCRBV", "TRBV") |>
        str_extract("TRBV[0-9\\-]+"),
      # Standardize J gene
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
# 1b. Import supplementary table 2 (myelin-reactive clones)
# --------------------------------------------------------------------------
message("  [1b] Importing Sukenikova supp table 2 (myelin-reactive clones)...")

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
# 1c. Size normalization
# --------------------------------------------------------------------------
message("  [1c] Size normalization analysis...")

# Report size comparison
n_heming_trb <- n_distinct(tcr_metadata$TRB_CDR3)
n_suk_total <- n_distinct(sukenikova_bulk$amino_acid)

message(sprintf("    Heming single-cell unique TRB CDR3: %d", n_heming_trb))
message(sprintf("    Sukenikova bulk unique CDR3b: %d", n_suk_total))
message(sprintf("    Ratio: %.1fx", n_suk_total / n_heming_trb))

# Downsampling strategy: keep unique CDR3b per patient (deduplicated across samples),
# then cap at top N by total template count to match single-cell scale
sukenikova_for_gliph <- sukenikova_bulk |>
  # Deduplicate across samples: keep unique CDR3b per patient
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

# Cap per patient: keep top 1000 clonotypes by template count
max_per_patient <- 1000
sukenikova_for_gliph <- sukenikova_for_gliph |>
  dplyr::group_by(suk_patient) |>
  dplyr::slice_max(order_by = total_count, n = max_per_patient, with_ties = FALSE) |>
  dplyr::ungroup()

n_suk_downsampled <- n_distinct(sukenikova_for_gliph$amino_acid)
message(sprintf("    After top-%d/patient cap: %d unique CDR3b", max_per_patient,
                n_suk_downsampled))

# Save import summary
import_summary <- tibble::tibble(
  source = c("Heming (single-cell)", "Sukenikova (bulk, raw)",
             "Sukenikova (dedup)", "Sukenikova (downsampled)"),
  n_unique_cdr3b = c(n_heming_trb, n_suk_total, n_suk_dedup, n_suk_downsampled),
  description = c(
    "Unique TRB CDR3 from single-cell Seurat object",
    "All unique CDR3b across 22 Adaptive immunoseq files",
    "Unique CDR3b after per-patient deduplication",
    sprintf("Top %d clonotypes per patient by template count", max_per_patient)
  )
)

writexl::write_xlsx(
  list(
    size_comparison = import_summary,
    sukenikova_per_file = sukenikova_bulk |>
      dplyr::group_by(file_name, suk_patient, suk_timepoint, suk_tissue) |>
      dplyr::summarize(
        n_clonotypes = n(),
        n_unique_cdr3 = n_distinct(amino_acid),
        total_templates = sum(count),
        .groups = "drop"
      )
  ),
  "results/tcr_comparison/sukenikova/sukenikova_import_summary.xlsx"
)
message("    Saved import summary to sukenikova/sukenikova_import_summary.xlsx")


# ==============================================================================
# Section 2: V/J Gene Usage Analysis
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Section 2: V/J Gene Usage Analysis")
message(paste(rep("=", 60), collapse = ""))

gene_usage_stats <- list()
gene_usage_props <- list()

for (chain in c("TRA", "TRB")) {
  for (gene_type in c("V", "J")) {
    message(sprintf("  Analyzing %s %s genes", chain, gene_type))

    # Extract gene from CTgene column
    gene_data <- tcr_metadata |>
      dplyr::mutate(
        gene = case_when(
          chain == "TRA" & gene_type == "V" ~ gsub("\\..*", "", gsub("_.*", "", CTgene)),
          chain == "TRA" & gene_type == "J" ~ gsub(".*\\.", "", gsub("_.*", "", CTgene)),
          chain == "TRB" & gene_type == "V" ~ gsub("\\..*", "", gsub(".*_", "", CTgene)),
          chain == "TRB" & gene_type == "J" ~ gsub(".*\\.", "", gsub(".*_", "", CTgene))
        )
      ) |>
      dplyr::filter(!is.na(gene), gene != "", nchar(gene) > 2)

    # Top 25 genes by frequency
    top_genes <- gene_data |>
      dplyr::count(gene, sort = TRUE) |>
      dplyr::slice(1:25) |>
      dplyr::pull(gene)

    # Calculate proportions per sample
    sample_gene_props <- gene_data |>
      dplyr::filter(gene %in% top_genes) |>
      dplyr::count(sample, patient, tissue, diagnosis, group, gene) |>
      dplyr::group_by(sample) |>
      dplyr::mutate(prop = n / sum(n)) |>
      dplyr::ungroup()

    gene_usage_props[[paste(chain, gene_type, sep = "_")]] <- sample_gene_props

    # Run mixed model for each gene: prop ~ diagnosis * tissue + (1|patient)
    gene_results <- lapply(top_genes, function(g) {
      gene_df <- sample_gene_props |>
        dplyr::filter(gene == g) |>
        tidyr::complete(
          tidyr::nesting(sample, patient, tissue, diagnosis, group),
          fill = list(n = 0, prop = 0)
        )

      if (nrow(gene_df) < 10) return(NULL)

      result <- run_diagnosis_stats(gene_df, "prop")
      result$gene <- g
      result$chain <- chain
      result$gene_type <- gene_type
      result
    }) |> dplyr::bind_rows()

    if (nrow(gene_results) > 0) {
      gene_results <- gene_results |>
        dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))
      gene_usage_stats[[paste(chain, gene_type, sep = "_")]] <- gene_results
    }
  }
}

gene_usage_combined <- dplyr::bind_rows(gene_usage_stats, .id = "comparison")

# Volcano plots per diagnosis and chain/gene_type
message("  Creating volcano plots")
for (dx in neuropathy_dx) {
  for (comp in unique(gene_usage_combined$comparison)) {
    dx_data <- gene_usage_combined |>
      dplyr::filter(diagnosis == dx, comparison == comp,
                    !is.na(estimate), !is.na(p.value))

    if (nrow(dx_data) < 3) next

    comp_label <- gsub("_", " ", comp)

    p <- ggplot(dx_data, aes(x = estimate, y = -log10(p.value))) +
      geom_point(aes(color = p.adj < 0.1), size = 2, alpha = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggrepel::geom_text_repel(
        data = dx_data |> dplyr::filter(p.adj < 0.1),
        aes(label = gene), size = 3, max.overlaps = 20
      ) +
      scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red3")) +
      theme_pub() +
      labs(title = paste(dx, "vs CTRL:", comp_label, "gene usage"),
           x = "Difference in proportion", y = "-log10(p-value)") +
      theme(legend.position = "none")

    ggsave(paste0("results/tcr_comparison/gene_usage/volcano_", dx, "_", comp, ".pdf"),
           p, width = 6, height = 5)
  }
}

# Significant genes heatmap
sig_genes <- gene_usage_combined |>
  dplyr::filter(p.adj < 0.1) |>
  dplyr::distinct(gene) |>
  dplyr::pull(gene)

if (length(sig_genes) > 3) {
  heatmap_data <- gene_usage_combined |>
    dplyr::filter(gene %in% sig_genes) |>
    dplyr::select(comparison, gene, diagnosis, estimate) |>
    dplyr::mutate(gene_comp = paste(comparison, gene, sep = "_")) |>
    dplyr::select(gene_comp, diagnosis, estimate) |>
    tidyr::pivot_wider(names_from = diagnosis, values_from = estimate, values_fill = 0) |>
    tibble::column_to_rownames("gene_comp") |>
    as.matrix()

  pdf("results/tcr_comparison/gene_usage/significant_genes_heatmap.pdf",
      width = 12, height = 10)
  pheatmap::pheatmap(
    heatmap_data,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "Significant V/J genes (difference from CTRL)",
    cluster_cols = FALSE
  )
  dev.off()
}

writexl::write_xlsx(gene_usage_stats,
                    "results/tcr_comparison/gene_usage/gene_usage_statistics.xlsx")
message("  Saved gene usage statistics")


# ==============================================================================
# Section 3: GLIPH2 Combined Analysis
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Section 3: GLIPH2 Combined Analysis (Heming + Sukenikova)")
message(paste(rep("=", 60), collapse = ""))

# --------------------------------------------------------------------------
# 3a. Prepare combined GLIPH2 input
# --------------------------------------------------------------------------
message("  [3a] Preparing combined GLIPH2 input...")

# Extract TRB CDR3 from single-cell data
heming_trb <- tcr_metadata |>
  dplyr::select(cell_id, TRB_CDR3, sample, patient, tissue, diagnosis) |>
  dplyr::filter(!is.na(TRB_CDR3), TRB_CDR3 != "") |>
  dplyr::rename(CDR3b = TRB_CDR3) |>
  dplyr::mutate(source = "Heming")

# Extract V gene from CTgene for Heming data
heming_vgene <- tcr_metadata |>
  dplyr::mutate(TRBV = gsub("\\..*", "", gsub(".*_", "", CTgene))) |>
  dplyr::select(cell_id, TRBV)

heming_trb <- heming_trb |>
  dplyr::left_join(heming_vgene, by = "cell_id")

# Sukenikova data for GLIPH (downsampled)
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

# Combine
combined_trb <- dplyr::bind_rows(heming_trb, suk_trb)
message(sprintf("    Combined input: %d sequences (%d Heming, %d Sukenikova)",
                nrow(combined_trb), nrow(heming_trb), nrow(suk_trb)))

# Build turboGliph input data.frame
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

# Total cells per diagnosis
total_by_diagnosis <- cluster_meta |>
  dplyr::filter(!is.na(diagnosis)) |>
  dplyr::count(diagnosis, name = "total")

# Fisher's exact test: each diagnosis vs CTRL for each cluster
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

  # Report top hits
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

# Standardize Sukenikova tissue labels to match Heming
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

  # Export tissue-specificity results
  writexl::write_xlsx(
    list(
      cluster_tissue_dx = cluster_tissue_dx,
      tissue_enrichment = tissue_enrichment,
      tissue_bias_summary = tissue_bias_summary
    ),
    "results/tcr_comparison/gliph/tissue_specificity_analysis.xlsx"
  )
  message("    Saved tissue-specificity analysis")
}

# --------------------------------------------------------------------------
# 3e. CDR3 motif extraction from disease-enriched clusters
# --------------------------------------------------------------------------
message("  [3e] Extracting driving motifs from GLIPH clusters...")

motif_driving <- NULL
if (nrow(diagnosis_enrichment) > 0) {
  dx_with_enrichment <- unique(diagnosis_enrichment$diagnosis)

  motif_driving <- lapply(dx_with_enrichment, function(dx) {
    dx_enriched_cl <- diagnosis_enrichment |>
      dplyr::filter(diagnosis == dx, odds_ratio > 1, p.adj < 0.2) |>
      dplyr::pull(cluster)

    dx_cdr3 <- cluster_meta |>
      dplyr::filter(cluster %in% dx_enriched_cl, !is.na(CDR3b)) |>
      dplyr::pull(CDR3b) |> unique()

    if (length(dx_cdr3) < 3) return(NULL)

    # Extract 3-mer and 4-mer motifs
    lapply(c(3, 4), function(k) {
      dx_motifs <- extract_cdr3_motifs(dx_cdr3, k)
      all_cdr3 <- cluster_meta |> dplyr::filter(!is.na(CDR3b)) |>
        dplyr::pull(CDR3b) |> unique()
      bg_motifs <- extract_cdr3_motifs(all_cdr3, k)

      shared_motifs <- intersect(names(dx_motifs), names(bg_motifs))
      if (length(shared_motifs) == 0) return(NULL)

      motif_stats <- lapply(shared_motifs, function(m) {
        a <- as.integer(dx_motifs[m])
        b <- sum(dx_motifs) - a
        c_val <- as.integer(bg_motifs[m])
        d <- sum(bg_motifs) - c_val
        mat <- matrix(c(a, c_val, b, d), nrow = 2)
        ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
        if (is.null(ft)) return(NULL)
        data.frame(
          diagnosis = dx, motif = m, k = k,
          count_enriched = a, count_background = c_val,
          freq_enriched = a / sum(dx_motifs),
          freq_background = c_val / sum(bg_motifs),
          odds_ratio = ft$estimate, p.value = ft$p.value
        )
      }) |> dplyr::bind_rows()

      if (nrow(motif_stats) > 0) {
        motif_stats$p.adj <- p.adjust(motif_stats$p.value, method = "BH")
        motif_stats$log2_fc <- log2((motif_stats$freq_enriched + 1e-6) /
                                     (motif_stats$freq_background + 1e-6))
      }
      motif_stats
    }) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()

  if (!is.null(motif_driving) && nrow(motif_driving) > 0) {
    writexl::write_xlsx(
      motif_driving |> dplyr::arrange(p.adj),
      "results/tcr_comparison/gliph/motif_enrichment_by_diagnosis.xlsx"
    )
    message("    Saved motif enrichment by diagnosis")
  }
}

# Tissue-specific motifs (CSF vs PBMC within disease clusters)
message("  Analyzing tissue-specific motifs...")

if (nrow(tissue_enrichment) > 0) {
  csf_clusters <- tissue_enrichment$cluster[tissue_enrichment$tissue_bias == "CSF-enriched"]
  pbmc_clusters <- tissue_enrichment$cluster[tissue_enrichment$tissue_bias == "PBMC-enriched"]

  csf_cdr3 <- cluster_meta |>
    dplyr::filter(cluster %in% csf_clusters, !is.na(CDR3b)) |>
    dplyr::pull(CDR3b) |> unique()
  pbmc_cdr3 <- cluster_meta |>
    dplyr::filter(cluster %in% pbmc_clusters, !is.na(CDR3b)) |>
    dplyr::pull(CDR3b) |> unique()

  tissue_motif_comparison <- NULL
  if (length(csf_cdr3) >= 5 && length(pbmc_cdr3) >= 5) {
    tissue_motif_comparison <- lapply(c(3, 4), function(k) {
      csf_m <- extract_cdr3_motifs(csf_cdr3, k)
      pbmc_m <- extract_cdr3_motifs(pbmc_cdr3, k)
      shared <- intersect(names(csf_m), names(pbmc_m))
      if (length(shared) == 0) return(NULL)

      lapply(shared, function(m) {
        a <- as.integer(csf_m[m]); b <- sum(csf_m) - a
        c_val <- as.integer(pbmc_m[m]); d <- sum(pbmc_m) - c_val
        ft <- tryCatch(fisher.test(matrix(c(a, c_val, b, d), nrow = 2)),
                       error = function(e) NULL)
        if (is.null(ft)) return(NULL)
        data.frame(
          motif = m, k = k,
          count_csf = a, count_pbmc = c_val,
          freq_csf = a / sum(csf_m), freq_pbmc = c_val / sum(pbmc_m),
          odds_ratio = ft$estimate, p.value = ft$p.value
        )
      }) |> dplyr::bind_rows()
    }) |> dplyr::bind_rows()

    if (!is.null(tissue_motif_comparison) && nrow(tissue_motif_comparison) > 0) {
      tissue_motif_comparison$p.adj <- p.adjust(tissue_motif_comparison$p.value, method = "BH")
      tissue_motif_comparison$log2_fc <- log2((tissue_motif_comparison$freq_csf + 1e-6) /
                                               (tissue_motif_comparison$freq_pbmc + 1e-6))
      writexl::write_xlsx(
        tissue_motif_comparison |> dplyr::arrange(p.adj),
        "results/tcr_comparison/gliph/tissue_specific_motif_comparison.xlsx"
      )
      message("    Saved tissue-specific motif comparison")
    }
  }
}

# --------------------------------------------------------------------------
# 3f. Cross-cohort cluster analysis
# --------------------------------------------------------------------------
message("  [3f] Cross-cohort cluster analysis...")

# Identify clusters containing TCRs from both datasets
cross_cohort <- cluster_meta |>
  dplyr::filter(!is.na(source), !is.na(cluster)) |>
  dplyr::group_by(cluster) |>
  dplyr::summarize(
    has_heming = any(source == "Heming"),
    has_sukenikova = any(source == "Sukenikova"),
    n_heming = sum(source == "Heming"),
    n_sukenikova = sum(source == "Sukenikova"),
    n_total = n(),
    heming_diagnoses = paste(unique(diagnosis[source == "Heming"]), collapse = ","),
    .groups = "drop"
  ) |>
  dplyr::mutate(cross_cohort = has_heming & has_sukenikova)

n_cross <- sum(cross_cohort$cross_cohort)
message(sprintf("    Cross-cohort clusters (both datasets): %d / %d", n_cross,
                nrow(cross_cohort)))

# Which Heming diagnoses share clusters with Sukenikova?
if (n_cross > 0) {
  cross_dx <- cross_cohort |>
    dplyr::filter(cross_cohort) |>
    dplyr::select(cluster, heming_diagnoses, n_heming, n_sukenikova, n_total) |>
    tidyr::separate_rows(heming_diagnoses, sep = ",") |>
    dplyr::count(heming_diagnoses, name = "n_shared_clusters") |>
    dplyr::arrange(desc(n_shared_clusters))

  message("    Heming diagnoses sharing clusters with Sukenikova GBS:")
  for (i in seq_len(nrow(cross_dx))) {
    message(sprintf("      %s: %d clusters", cross_dx$heming_diagnoses[i],
                    cross_dx$n_shared_clusters[i]))
  }
}

# Export enriched CDR3 examples
top_enriched <- diagnosis_enrichment |>
  dplyr::filter(p.adj < 0.2, odds_ratio > 2) |>
  dplyr::group_by(diagnosis) |>
  dplyr::slice_max(order_by = odds_ratio, n = 5) |>
  dplyr::ungroup()

if (nrow(top_enriched) > 0) {
  enriched_cdr3 <- cluster_meta |>
    dplyr::filter(cluster %in% top_enriched$cluster) |>
    dplyr::left_join(top_enriched |> dplyr::select(cluster, diagnosis, odds_ratio, p.adj),
                     by = "cluster", suffix = c("", "_enriched")) |>
    dplyr::filter(!is.na(diagnosis_enriched)) |>
    dplyr::group_by(cluster, diagnosis_enriched) |>
    dplyr::summarize(
      n_sequences = n(),
      example_CDR3s = paste(head(unique(CDR3b), 5), collapse = ", "),
      odds_ratio = first(odds_ratio),
      p.adj = first(p.adj),
      .groups = "drop"
    )

  writexl::write_xlsx(
    enriched_cdr3,
    "results/tcr_comparison/gliph/disease_enriched_cdr3_motifs.xlsx"
  )
  message("    Saved disease-enriched CDR3 motifs")
}

writexl::write_xlsx(
  list(
    cross_cohort_summary = cross_cohort,
    cross_cohort_by_diagnosis = if (n_cross > 0) cross_dx else tibble::tibble(),
    diagnosis_enrichment = diagnosis_enrichment
  ),
  "results/tcr_comparison/gliph/cross_cohort_cluster_analysis.xlsx"
)
message("    Saved cross-cohort cluster analysis")

# --------------------------------------------------------------------------
# 3g. Cluster composition profiling
# --------------------------------------------------------------------------
message("  [3g] Cluster composition profiling...")

cluster_composition <- cluster_meta |>
  dplyr::filter(!is.na(diagnosis), !is.na(cluster)) |>
  dplyr::count(cluster, diagnosis) |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(
    fraction = n / sum(n),
    cluster_total = sum(n)
  ) |>
  dplyr::ungroup()

cluster_entropy <- cluster_composition |>
  dplyr::group_by(cluster) |>
  dplyr::summarize(
    cluster_total = first(cluster_total),
    n_diagnoses = n(),
    shannon_entropy = -sum(fraction * log2(fraction + 1e-10)),
    max_fraction = max(fraction),
    dominant_dx = diagnosis[which.max(fraction)],
    .groups = "drop"
  ) |>
  dplyr::mutate(
    specificity_class = dplyr::case_when(
      max_fraction >= 0.9 ~ "Private",
      max_fraction >= 0.6 ~ "Semi-private",
      TRUE ~ "Promiscuous"
    )
  )

enriched_cluster_ids <- diagnosis_enrichment |>
  dplyr::filter(p.adj < 0.1, odds_ratio > 1) |>
  dplyr::pull(cluster) |> unique()

enriched_cluster_entropy <- cluster_entropy |>
  dplyr::filter(cluster %in% enriched_cluster_ids)

specificity_summary <- enriched_cluster_entropy |>
  dplyr::left_join(
    diagnosis_enrichment |>
      dplyr::filter(p.adj < 0.1, odds_ratio > 1) |>
      dplyr::select(cluster, diagnosis) |>
      dplyr::rename(enriched_for = diagnosis),
    by = "cluster"
  ) |>
  dplyr::count(enriched_for, specificity_class) |>
  tidyr::pivot_wider(names_from = specificity_class, values_from = n, values_fill = 0)

message(sprintf("    Enriched clusters: %d private, %d semi-private, %d promiscuous",
                sum(enriched_cluster_entropy$specificity_class == "Private"),
                sum(enriched_cluster_entropy$specificity_class == "Semi-private"),
                sum(enriched_cluster_entropy$specificity_class == "Promiscuous")))

writexl::write_xlsx(
  list(
    cluster_composition = cluster_composition,
    cluster_entropy = cluster_entropy,
    enriched_cluster_entropy = enriched_cluster_entropy,
    specificity_summary = specificity_summary
  ),
  "results/tcr_comparison/gliph/cluster_composition_profiling.xlsx"
)
message("    Saved cluster composition profiling")

# --------------------------------------------------------------------------
# 3h. CDR3 length distribution in enriched clusters
# --------------------------------------------------------------------------
message("  [3h] CDR3 length distributions in enriched clusters...")

bg_cdr3_lengths <- cluster_meta |>
  dplyr::filter(!is.na(CDR3b)) |>
  dplyr::mutate(cdr3_length = nchar(CDR3b))

enriched_cdr3_lengths <- cluster_meta |>
  dplyr::filter(cluster %in% enriched_cluster_ids, !is.na(CDR3b)) |>
  dplyr::left_join(
    diagnosis_enrichment |>
      dplyr::filter(p.adj < 0.1, odds_ratio > 1) |>
      dplyr::select(cluster, diagnosis) |>
      dplyr::rename(enriched_for = diagnosis),
    by = "cluster"
  ) |>
  dplyr::filter(!is.na(enriched_for)) |>
  dplyr::mutate(cdr3_length = nchar(CDR3b))

cdr3_length_tests <- lapply(unique(enriched_cdr3_lengths$enriched_for), function(dx) {
  dx_lengths <- enriched_cdr3_lengths |>
    dplyr::filter(enriched_for == dx) |>
    dplyr::pull(cdr3_length)
  bg_lengths <- bg_cdr3_lengths$cdr3_length

  if (length(dx_lengths) < 5) return(NULL)

  ks <- ks.test(dx_lengths, bg_lengths)
  data.frame(
    diagnosis = dx,
    n_enriched = length(dx_lengths),
    mean_enriched = mean(dx_lengths),
    sd_enriched = sd(dx_lengths),
    mean_background = mean(bg_lengths),
    sd_background = sd(bg_lengths),
    ks_statistic = ks$statistic,
    ks_pvalue = ks$p.value
  )
}) |> dplyr::bind_rows()

if (nrow(cdr3_length_tests) > 0) {
  cdr3_length_tests$ks_padj <- p.adjust(cdr3_length_tests$ks_pvalue, method = "BH")
  message("    CDR3 length constraint results:")
  print(cdr3_length_tests |> dplyr::select(diagnosis, mean_enriched, mean_background, ks_padj))
}

writexl::write_xlsx(
  list(
    cdr3_length_tests = cdr3_length_tests,
    enriched_cdr3_lengths = enriched_cdr3_lengths |>
      dplyr::count(enriched_for, cdr3_length)
  ),
  "results/tcr_comparison/gliph/cdr3_length_analysis.xlsx"
)
message("    Saved CDR3 length analysis")

# --------------------------------------------------------------------------
# 3i. V-gene usage within GLIPH clusters
# --------------------------------------------------------------------------
message("  [3i] V-gene usage within GLIPH clusters...")

cluster_vgene <- cluster_meta |>
  dplyr::filter(cluster %in% enriched_cluster_ids, !is.na(CDR3b)) |>
  dplyr::left_join(
    combined_trb |> dplyr::select(CDR3b, TRBV) |> dplyr::distinct(),
    by = "CDR3b"
  ) |>
  dplyr::filter(!is.na(TRBV), TRBV != "") |>
  dplyr::left_join(
    diagnosis_enrichment |>
      dplyr::filter(p.adj < 0.1, odds_ratio > 1) |>
      dplyr::select(cluster, diagnosis) |>
      dplyr::rename(enriched_for = diagnosis),
    by = "cluster"
  ) |>
  dplyr::filter(!is.na(enriched_for))

bg_vgene_freq <- combined_trb |>
  dplyr::filter(!is.na(TRBV), TRBV != "") |>
  dplyr::count(TRBV) |>
  dplyr::mutate(bg_freq = n / sum(n)) |>
  dplyr::rename(bg_n = n)

vgene_in_clusters <- cluster_vgene |>
  dplyr::count(enriched_for, TRBV) |>
  dplyr::group_by(enriched_for) |>
  dplyr::mutate(cluster_freq = n / sum(n), cluster_total = sum(n)) |>
  dplyr::ungroup() |>
  dplyr::left_join(bg_vgene_freq, by = "TRBV") |>
  dplyr::mutate(
    fold_enrichment = (cluster_freq + 1e-6) / (bg_freq + 1e-6),
    log2_fe = log2(fold_enrichment)
  )

vgene_enrichment_tests <- lapply(
  seq_len(nrow(vgene_in_clusters)),
  function(i) {
    row <- vgene_in_clusters[i, ]
    a <- row$n
    b <- row$cluster_total - a
    c_val <- row$bg_n
    d <- sum(bg_vgene_freq$bg_n) - c_val
    if (any(c(a, b, c_val, d) < 0, na.rm = TRUE)) return(NULL)
    ft <- tryCatch(
      fisher.test(matrix(c(a, c_val, b, d), nrow = 2)),
      error = function(e) NULL
    )
    if (is.null(ft)) return(NULL)
    data.frame(
      enriched_for = row$enriched_for, TRBV = row$TRBV,
      n_in_clusters = a, cluster_freq = row$cluster_freq,
      bg_freq = row$bg_freq, fold_enrichment = row$fold_enrichment,
      OR = ft$estimate, p.value = ft$p.value
    )
  }
) |> dplyr::bind_rows()

if (nrow(vgene_enrichment_tests) > 0) {
  vgene_enrichment_tests$p.adj <- p.adjust(vgene_enrichment_tests$p.value, method = "BH")
  sig_vgenes <- vgene_enrichment_tests |> dplyr::filter(p.adj < 0.1)
  message(sprintf("    %d significant V-gene enrichments within GLIPH clusters", nrow(sig_vgenes)))
}

writexl::write_xlsx(
  list(
    vgene_in_clusters = vgene_in_clusters,
    vgene_enrichment_tests = vgene_enrichment_tests
  ),
  "results/tcr_comparison/gliph/vgene_within_clusters.xlsx"
)
message("    Saved V-gene within clusters analysis")

# --------------------------------------------------------------------------
# 3j. Public clone analysis
# --------------------------------------------------------------------------
message("  [3j] Public clone analysis...")

public_clones <- combined_trb |>
  dplyr::filter(!is.na(CDR3b), CDR3b != "", !is.na(diagnosis)) |>
  dplyr::distinct(CDR3b, patient, diagnosis) |>
  dplyr::group_by(CDR3b, diagnosis) |>
  dplyr::summarize(
    n_patients = dplyr::n_distinct(patient),
    patients = paste(unique(patient), collapse = ","),
    .groups = "drop"
  ) |>
  dplyr::filter(n_patients >= 2)

total_unique_cdr3_per_dx <- combined_trb |>
  dplyr::filter(!is.na(CDR3b), !is.na(diagnosis)) |>
  dplyr::distinct(CDR3b, diagnosis) |>
  dplyr::count(diagnosis, name = "total_unique")

public_clone_rates <- public_clones |>
  dplyr::count(diagnosis, name = "n_public") |>
  dplyr::left_join(total_unique_cdr3_per_dx, by = "diagnosis") |>
  dplyr::mutate(public_rate = n_public / total_unique)

message("    Public clone rates:")
print(public_clone_rates)

public_in_gliph <- public_clones |>
  dplyr::inner_join(
    cluster_meta |>
      dplyr::filter(cluster %in% enriched_cluster_ids) |>
      dplyr::select(CDR3b, cluster) |>
      dplyr::distinct(),
    by = "CDR3b"
  )

message(sprintf("    %d public clones found in disease-enriched GLIPH clusters",
                dplyr::n_distinct(public_in_gliph$CDR3b)))

suk_public <- myelin_lookup |>
  dplyr::filter(public_clonotype == "yes") |>
  dplyr::pull(CDR3b) |> unique()

public_suk_overlap <- public_clones |>
  dplyr::filter(CDR3b %in% suk_public)

message(sprintf("    %d public clones overlap with Sukenikova public clonotypes",
                nrow(public_suk_overlap)))

writexl::write_xlsx(
  list(
    public_clones = public_clones,
    public_clone_rates = public_clone_rates,
    public_in_gliph = public_in_gliph,
    public_suk_overlap = public_suk_overlap
  ),
  "results/tcr_comparison/gliph/public_clone_analysis.xlsx"
)
message("    Saved public clone analysis")

# ==============================================================================
# Section 4: Sukenikova Cross-Reference Analysis
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Section 4: Sukenikova Cross-Reference Analysis")
message(paste(rep("=", 60), collapse = ""))

# --------------------------------------------------------------------------
# 4a. Myelin-reactive clone mapping
# --------------------------------------------------------------------------
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

# GLIPH cluster-based matching: Heming CDR3b in same cluster as myelin-reactive clone
myelin_cdr3_set <- unique(myelin_lookup$CDR3b)

# Find clusters containing any myelin-reactive clone
myelin_clusters <- cluster_meta |>
  dplyr::filter(CDR3b %in% myelin_cdr3_set) |>
  dplyr::distinct(cluster) |>
  dplyr::pull(cluster)

# Find Heming cells in those clusters (but not exact matches)
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

# Quantify myelin-reactive % by diagnosis
myelin_by_dx <- dplyr::bind_rows(
  exact_matches |>
    dplyr::mutate(match_type = "exact") |>
    dplyr::select(diagnosis, tissue, match_type),
  cluster_matches |>
    dplyr::mutate(match_type = "cluster") |>
    dplyr::select(diagnosis, tissue, match_type)
) |>
  dplyr::count(diagnosis, match_type) |>
  dplyr::left_join(
    tcr_metadata |> dplyr::count(diagnosis, name = "total_cells"),
    by = "diagnosis"
  ) |>
  dplyr::mutate(fraction = n / total_cells)

# Fisher's exact: is any diagnosis enriched for myelin-reactive TCRs vs CTRL?
all_myelin_matches <- unique(c(exact_matches$cell_id, cluster_matches$CDR3b))

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

# --------------------------------------------------------------------------
# 4b. Mirror Sukenikova GLIPH analysis
# --------------------------------------------------------------------------
message("  [4b] Mirroring Sukenikova GLIPH cluster assignments...")

# Check if Heming CDR3b sequences fall into their published GLIPH clusters
# Their clusters are identified by gliph_cluster_id in supp table 2
suk_cluster_members <- myelin_lookup |>
  dplyr::filter(!is.na(gliph_cluster_id)) |>
  dplyr::group_by(gliph_cluster_id) |>
  dplyr::summarize(
    cluster_cdr3s = list(unique(CDR3b)),
    n_members = n(),
    specificities = paste(unique(specificity), collapse = ","),
    .groups = "drop"
  )

# For each Sukenikova cluster, check if any of our CDR3b are exact matches
# or fall into the same de novo GLIPH2 cluster
heming_cdr3_set <- unique(heming_trb$CDR3b)

suk_cluster_overlap <- suk_cluster_members |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_heming_exact = sum(cluster_cdr3s %in% heming_cdr3_set),
    heming_exact_cdr3 = paste(intersect(unlist(cluster_cdr3s), heming_cdr3_set),
                              collapse = ",")
  ) |>
  dplyr::ungroup() |>
  dplyr::select(-cluster_cdr3s)

n_overlap_clusters <- sum(suk_cluster_overlap$n_heming_exact > 0)
message(sprintf("    %d / %d Sukenikova GLIPH clusters have exact Heming CDR3b matches",
                n_overlap_clusters, nrow(suk_cluster_overlap)))

# Compare de novo vs published cluster assignments
# Find CDR3b that are in both our GLIPH2 and their published clusters
shared_cdr3_denovo <- cluster_meta |>
  dplyr::filter(CDR3b %in% myelin_cdr3_set, !is.na(cluster)) |>
  dplyr::select(CDR3b, cluster) |>
  dplyr::rename(denovo_cluster = cluster) |>
  dplyr::left_join(
    myelin_lookup |> dplyr::select(CDR3b, gliph_cluster_id),
    by = "CDR3b"
  ) |>
  dplyr::rename(published_cluster = gliph_cluster_id) |>
  dplyr::distinct()

message(sprintf("    %d CDR3b sequences have both de novo and published cluster assignments",
                nrow(shared_cdr3_denovo)))

# --------------------------------------------------------------------------
# 4c. Antigen specificity enrichment
# --------------------------------------------------------------------------
message("  [4c] Antigen specificity enrichment...")

# For all matches (exact + cluster), what antigens are they reactive to?
all_matched_cdr3 <- unique(c(exact_matches$TRB_CDR3, cluster_matches$CDR3b))

antigen_mapping <- myelin_lookup |>
  dplyr::filter(CDR3b %in% myelin_cdr3_set) |>
  dplyr::select(CDR3b, specificity, pns_myelin_antigen)

# Map through exact + cluster matches to diagnoses
antigen_by_dx <- dplyr::bind_rows(
  exact_matches |>
    dplyr::select(TRB_CDR3, diagnosis, tissue, specificity) |>
    dplyr::rename(CDR3b = TRB_CDR3) |>
    dplyr::mutate(match_type = "exact"),
  cluster_matches |>
    dplyr::left_join(
      cluster_meta |>
        dplyr::filter(CDR3b %in% myelin_cdr3_set) |>
        dplyr::select(cluster, CDR3b) |>
        dplyr::left_join(myelin_lookup |> dplyr::select(CDR3b, specificity), by = "CDR3b") |>
        dplyr::select(cluster, specificity) |>
        dplyr::distinct(),
      by = "cluster"
    ) |>
    dplyr::mutate(match_type = "cluster")
)

if (nrow(antigen_by_dx) > 0) {
  antigen_summary <- antigen_by_dx |>
    dplyr::filter(!is.na(specificity)) |>
    dplyr::count(diagnosis, specificity, match_type) |>
    dplyr::arrange(diagnosis, desc(n))

  message("    Antigen specificity by diagnosis:")
  print(antigen_summary)
}

# --------------------------------------------------------------------------
# 4d. Hamming distance near-matches
# --------------------------------------------------------------------------
message("  [4d] Hamming distance near-matches to myelin-reactive CDR3b...")

hamming_distance <- function(s1, s2) {
  if (nchar(s1) != nchar(s2)) return(Inf)
  chars1 <- strsplit(s1, "")[[1]]
  chars2 <- strsplit(s2, "")[[1]]
  sum(chars1 != chars2)
}

heming_cdr3_unique <- unique(heming_for_match$TRB_CDR3)

myelin_by_length <- split(myelin_cdr3_set, nchar(myelin_cdr3_set))
heming_by_length <- split(heming_cdr3_unique, nchar(heming_cdr3_unique))

near_matches <- lapply(names(myelin_by_length), function(len) {
  myelin_seqs <- myelin_by_length[[len]]
  heming_seqs <- heming_by_length[[len]]
  if (is.null(heming_seqs) || length(heming_seqs) == 0) return(NULL)

  results <- lapply(myelin_seqs, function(mseq) {
    distances <- sapply(heming_seqs, function(hseq) hamming_distance(mseq, hseq))
    near_idx <- which(distances == 1)
    if (length(near_idx) == 0) return(NULL)
    data.frame(
      myelin_CDR3b = mseq,
      heming_CDR3b = heming_seqs[near_idx],
      hamming_distance = 1,
      cdr3_length = as.integer(len)
    )
  }) |> dplyr::bind_rows()

  results
}) |> dplyr::bind_rows()

n_near <- if (!is.null(near_matches) && nrow(near_matches) > 0) nrow(near_matches) else 0
message(sprintf("    Hamming<=1 near-matches: %d", n_near))

if (n_near > 0) {
  near_matches_annotated <- near_matches |>
    dplyr::left_join(
      heming_for_match |>
        dplyr::select(TRB_CDR3, patient, tissue, diagnosis) |>
        dplyr::distinct(),
      by = c("heming_CDR3b" = "TRB_CDR3")
    ) |>
    dplyr::left_join(
      myelin_lookup |> dplyr::select(CDR3b, specificity, pns_myelin_antigen) |> dplyr::distinct(),
      by = c("myelin_CDR3b" = "CDR3b")
    )

  near_match_enrichment <- lapply(neuropathy_dx, function(dx) {
    dx_total <- sum(tcr_metadata$diagnosis == dx)
    ctrl_total <- sum(tcr_metadata$diagnosis == "CTRL")
    dx_match <- sum(near_matches_annotated$diagnosis == dx, na.rm = TRUE) +
      sum(exact_matches$diagnosis == dx) +
      sum(cluster_matches$diagnosis == dx)
    ctrl_match <- sum(near_matches_annotated$diagnosis == "CTRL", na.rm = TRUE) +
      sum(exact_matches$diagnosis == "CTRL") +
      sum(cluster_matches$diagnosis == "CTRL")

    mat <- matrix(c(dx_match, ctrl_match,
                     dx_total - dx_match, ctrl_total - ctrl_match), nrow = 2)
    ft <- tryCatch(fisher.test(mat), error = function(e) list(estimate = NA, p.value = NA))

    data.frame(
      diagnosis = dx,
      n_all_matches = dx_match,
      n_total = dx_total,
      fraction = dx_match / dx_total,
      odds_ratio = ft$estimate,
      p.value = ft$p.value
    )
  }) |> dplyr::bind_rows()
  near_match_enrichment$p.adj <- p.adjust(near_match_enrichment$p.value, method = "BH")

  message("    Enrichment including near-matches:")
  print(near_match_enrichment |> dplyr::select(diagnosis, n_all_matches, fraction, odds_ratio, p.adj))
} else {
  near_matches_annotated <- tibble::tibble()
  near_match_enrichment <- tibble::tibble()
}

# --------------------------------------------------------------------------
# 4e. Per-antigen cluster enrichment
# --------------------------------------------------------------------------
message("  [4e] Per-antigen cluster enrichment...")

antigens <- unique(na.omit(myelin_lookup$specificity))

antigen_cluster_enrichment <- lapply(antigens, function(ag) {
  ag_cdr3 <- myelin_lookup |>
    dplyr::filter(specificity == ag) |>
    dplyr::pull(CDR3b) |> unique()

  ag_clusters <- cluster_meta |>
    dplyr::filter(CDR3b %in% ag_cdr3, !is.na(cluster)) |>
    dplyr::distinct(cluster) |>
    dplyr::pull(cluster)

  if (length(ag_clusters) == 0) return(NULL)

  lapply(ag_clusters, function(cl) {
    cl_data <- cluster_meta |> dplyr::filter(cluster == cl, !is.na(CDR3b))
    n_ag_in_cl <- sum(cl_data$CDR3b %in% ag_cdr3)
    n_total_cl <- nrow(cl_data)

    data.frame(
      antigen = ag,
      cluster = cl,
      n_antigen_in_cluster = n_ag_in_cl,
      cluster_size = n_total_cl,
      antigen_fraction = n_ag_in_cl / n_total_cl
    )
  }) |> dplyr::bind_rows()
}) |> dplyr::bind_rows()

if (nrow(antigen_cluster_enrichment) > 0) {
  antigen_cluster_enrichment <- antigen_cluster_enrichment |>
    dplyr::left_join(
      diagnosis_enrichment |>
        dplyr::filter(p.adj < 0.1, odds_ratio > 1) |>
        dplyr::select(cluster, diagnosis, odds_ratio, p.adj) |>
        dplyr::rename(dx_enriched = diagnosis, dx_OR = odds_ratio, dx_padj = p.adj),
      by = "cluster"
    )

  message(sprintf("    %d antigen-cluster associations found", nrow(antigen_cluster_enrichment)))
}

writexl::write_xlsx(
  antigen_cluster_enrichment,
  "results/tcr_comparison/sukenikova/antigen_cluster_enrichment.xlsx"
)
message("    Saved antigen cluster enrichment")

# --------------------------------------------------------------------------
# 4f. Sukenikova timepoint analysis
# --------------------------------------------------------------------------
message("  [4f] Sukenikova timepoint analysis...")

suk_by_timepoint <- sukenikova_bulk |>
  dplyr::filter(!is.na(suk_timepoint)) |>
  dplyr::distinct(amino_acid, suk_timepoint) |>
  dplyr::rename(CDR3b = amino_acid)

cluster_timepoint <- cluster_meta |>
  dplyr::filter(source == "Sukenikova", !is.na(cluster)) |>
  dplyr::left_join(suk_by_timepoint, by = "CDR3b") |>
  dplyr::filter(!is.na(suk_timepoint)) |>
  dplyr::count(cluster, suk_timepoint)

timepoint_by_dx <- cluster_meta |>
  dplyr::filter(cluster %in% enriched_cluster_ids, source == "Heming", !is.na(diagnosis)) |>
  dplyr::select(cluster, diagnosis) |>
  dplyr::distinct() |>
  dplyr::left_join(
    cluster_timepoint |>
      tidyr::pivot_wider(names_from = suk_timepoint, values_from = n, values_fill = 0),
    by = "cluster"
  )

if (!"AC" %in% names(timepoint_by_dx)) timepoint_by_dx$AC <- 0L
if (!"REC" %in% names(timepoint_by_dx)) timepoint_by_dx$REC <- 0L

timepoint_summary <- timepoint_by_dx |>
  dplyr::group_by(diagnosis) |>
  dplyr::summarize(
    n_clusters_with_suk = sum(AC > 0 | REC > 0),
    total_AC = sum(AC),
    total_REC = sum(REC),
    ac_fraction = total_AC / (total_AC + total_REC + 1e-10),
    .groups = "drop"
  )

message("    Timepoint analysis by diagnosis:")
print(timepoint_summary)

writexl::write_xlsx(
  list(
    cluster_timepoint = cluster_timepoint,
    timepoint_by_dx = timepoint_by_dx,
    timepoint_summary = timepoint_summary
  ),
  "results/tcr_comparison/sukenikova/timepoint_analysis.xlsx"
)
message("    Saved timepoint analysis")

# Export all cross-reference results
writexl::write_xlsx(
  list(
    exact_matches = exact_matches,
    cluster_matches = cluster_matches,
    near_matches = near_matches_annotated,
    myelin_by_diagnosis = myelin_by_dx,
    myelin_enrichment = myelin_enrichment,
    near_match_enrichment = near_match_enrichment,
    suk_cluster_overlap = suk_cluster_overlap,
    denovo_vs_published = shared_cdr3_denovo,
    antigen_by_diagnosis = if (nrow(antigen_by_dx) > 0) antigen_summary else tibble::tibble(),
    antigen_cluster_enrichment = antigen_cluster_enrichment
  ),
  "results/tcr_comparison/sukenikova/cross_reference_analysis.xlsx"
)
message("    Saved cross-reference analysis")

# Table visualization of exact matches
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

  pdf("results/tcr_comparison/sukenikova/sukenikova_reactive_table.pdf",
      width = 12, height = 8)
  grid::grid.draw(table_plot)
  dev.off()
  message("    Saved sukenikova_reactive_table.pdf")
}


# ==============================================================================
# Section 5: Save Cache for Figure Script
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Section 5: Saving analysis cache")
message(paste(rep("=", 60), collapse = ""))

all_dx <- c(neuropathy_dx, "GBS_Sukenikova")

gliph_cache <- list(
  # Metadata and constants
  diagnosis_col = diagnosis_col,
  tissue_col = tissue_col,
  neuropathy_dx = neuropathy_dx,
  all_dx = all_dx,
  dx_order = dx_order,

  # Section 2: V/J gene usage
  gene_usage_combined = gene_usage_combined,
  gene_usage_props = gene_usage_props,

  # Section 3: GLIPH core
  combined_trb = combined_trb,
  gliph_result = gliph_result,
  cluster_meta = cluster_meta,
  cluster_membership = cluster_membership,

  # Section 3c: Disease enrichment
  diagnosis_enrichment = diagnosis_enrichment,
  total_by_diagnosis = total_by_diagnosis,

  # Section 3d: Tissue
  tissue_enrichment = tissue_enrichment,
  cluster_tissue_dx = if (exists("cluster_tissue_dx")) cluster_tissue_dx else tibble::tibble(),
  tissue_bias_summary = if (exists("tissue_bias_summary")) tissue_bias_summary else tibble::tibble(),

  # Section 3e: Motifs
  motif_driving = motif_driving,
  tissue_motif_comparison = if (exists("tissue_motif_comparison")) tissue_motif_comparison else NULL,

  # Section 3f: Cross-cohort
  cross_cohort = cross_cohort,
  cross_dx = if (exists("cross_dx")) cross_dx else tibble::tibble(),
  enriched_cdr3 = if (exists("enriched_cdr3")) enriched_cdr3 else tibble::tibble(),

  # Section 3g: Cluster composition
  cluster_composition = cluster_composition,
  cluster_entropy = cluster_entropy,
  enriched_cluster_entropy = enriched_cluster_entropy,
  specificity_summary = specificity_summary,
  enriched_cluster_ids = enriched_cluster_ids,

  # Section 3h: CDR3 length
  enriched_cdr3_lengths = enriched_cdr3_lengths,
  bg_cdr3_lengths = bg_cdr3_lengths,
  cdr3_length_tests = cdr3_length_tests,

  # Section 3i: V-gene within clusters
  vgene_in_clusters = vgene_in_clusters,
  vgene_enrichment_tests = vgene_enrichment_tests,

  # Section 3j: Public clones
  public_clones = public_clones,
  public_clone_rates = public_clone_rates,
  public_in_gliph = public_in_gliph,
  public_suk_overlap = public_suk_overlap,

  # Section 4: Sukenikova cross-reference
  exact_matches = exact_matches,
  cluster_matches = cluster_matches,
  myelin_enrichment = myelin_enrichment,
  myelin_by_dx = myelin_by_dx,
  myelin_lookup = myelin_lookup,
  myelin_cdr3_set = myelin_cdr3_set,
  myelin_clusters = myelin_clusters,
  suk_cluster_overlap = suk_cluster_overlap,
  shared_cdr3_denovo = shared_cdr3_denovo,
  antigen_by_dx = antigen_by_dx,
  antigen_summary = if (exists("antigen_summary")) antigen_summary else tibble::tibble(),

  # Section 4d: Hamming near-matches
  near_matches_annotated = near_matches_annotated,
  near_match_enrichment = near_match_enrichment,

  # Section 4e: Per-antigen cluster enrichment
  antigen_cluster_enrichment = antigen_cluster_enrichment,

  # Section 4f: Timepoint
  timepoint_summary = timepoint_summary,
  timepoint_by_dx = timepoint_by_dx,
  cluster_timepoint = cluster_timepoint,

  # Sukenikova raw
  sukenikova_bulk = sukenikova_bulk,

  # Timestamp
  cache_timestamp = Sys.time()
)

qs::qsave(gliph_cache, "data/cache/gliph_results.qs", nthreads = 6)
message(sprintf("  Cache saved: data/cache/gliph_results.qs (%d objects)",
                length(gliph_cache)))

# ==============================================================================
# Done
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("Analysis complete!")
message(paste(rep("=", 60), collapse = ""))
message(sprintf("  Output directory: results/tcr_comparison/"))
message(sprintf("  Cache: data/cache/gliph_results.qs"))
message(sprintf("  Run scripts/tcr_gliph_figures.R for manuscript figures"))
