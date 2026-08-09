targets_preprocess <- list(
  tar_target(
    sc_list,
    create_sample_objects(
      rna_manifest,
      lookup,
      donor_assignments,
      rna_creation_config,
      rna_files
    )
  ),
  tar_target(
    sc_filtered,
    filter_sample_objects(sc_list, qc_thresholds_file, rna_filter_config)
  ),
  tar_target(sc_merge_pre, merge_sample_objects(sc_filtered, lookup)),
  tar_target(qc_count_cells, summarize_cell_counts(sc_list, sc_merge_pre)),
  tar_target(qc_count_genes, summarize_gene_counts(sc_merge_pre)),
  tar_target(qc_cellranger_metrics, summarize_cellranger_metrics(rna_manifest)),
  tar_target(
    qc_count_cells_file,
    write_csv_result(qc_count_cells, "results/targets/qc/count_cells.csv"),
    format = "file"
  ),
  tar_target(
    qc_count_genes_file,
    write_csv_result(qc_count_genes, "results/targets/qc/count_genes.csv"),
    format = "file"
  ),
  tar_target(
    qc_cellranger_metrics_file,
    write_csv_result(
      qc_cellranger_metrics,
      "results/targets/qc/cellranger_metrics.csv"
    ),
    format = "file"
  ),
  tar_target(
    sc_preprocessed, normalize_and_reduce(sc_merge_pre, rna_normalization_config)
  )
)
