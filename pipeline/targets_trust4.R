targets_trust4 <- list(
  tar_target(
    sural_metadata_umap_file,
    analysis_config$paths$sural_metadata_umap,
    format = "file"
  ),
  tar_target(
    sural_trust4_manifest,
    discover_sural_trust4_inputs(
      analysis_config$paths$sural_trust4,
      analysis_config$trust4$sample_map
    )
  ),
  tar_target(
    sural_trust4_files,
    sural_trust4_input_files(sural_trust4_manifest),
    format = "file"
  ),
  tar_target(
    sural_metadata_umap,
    read_sural_metadata_umap(
      sural_metadata_umap_file,
      analysis_config$trust4$reduction
    )
  ),
  tar_target(
    sural_trust4_mapped,
    map_sural_trust4_barcodes(
      sural_trust4_manifest,
      sural_trust4_files,
      sural_metadata_umap,
      analysis_config$trust4$cluster_column
    )
  ),
  tar_target(
    sural_trust4_table_file,
    write_sural_trust4_table(sural_trust4_mapped),
    format = "file"
  ),
  tar_target(
    sural_trust4_umap_file,
    write_sural_trust4_umap(
      sural_metadata_umap,
      sural_trust4_mapped
    ),
    format = "file"
  ),
  tar_target(
    sural_trust4_cluster_summary_file,
    write_sural_trust4_cluster_summary(sural_trust4_mapped),
    format = "file"
  ),
  tar_target(
    sural_ic_metadata_umap_file,
    analysis_config$paths$sural_ic_metadata_umap,
    format = "file"
  ),
  tar_target(
    sural_ic_metadata_umap,
    read_sural_metadata_umap(
      sural_ic_metadata_umap_file,
      analysis_config$trust4$ic_reduction
    )
  ),
  tar_target(
    sural_ic_trust4_mapped,
    map_sural_trust4_barcodes(
      sural_trust4_manifest,
      sural_trust4_files,
      sural_ic_metadata_umap,
      analysis_config$trust4$ic_cluster_column
    )
  ),
  tar_target(
    sural_ic_trust4_table_file,
    write_sural_trust4_table(
      sural_ic_trust4_mapped,
      "trust4_immune_cells_with_subclusters.xlsx"
    ),
    format = "file"
  ),
  tar_target(
    sural_ic_trust4_umap_file,
    write_sural_trust4_umap(
      sural_ic_metadata_umap,
      sural_ic_trust4_mapped,
      "trust4_immune_cells_umap.png",
      "Immune cells with TRUST4 receptor calls"
    ),
    format = "file"
  ),
  tar_target(
    sural_ic_trust4_cluster_summary_file,
    write_sural_trust4_cluster_summary(
      sural_ic_trust4_mapped,
      "trust4_receptor_cells_by_immune_subcluster.pdf",
      "TRUST4 receptor-positive cells by immune subcluster"
    ),
    format = "file"
  ),
  tar_target(
    sural_tcr_cross_tissue_comparison,
    prepare_sural_tcr_comparison(
      tcr_contigs,
      sural_trust4_mapped$records,
      analysis_config$trust4$patient_map
    )
  ),
  tar_target(
    sural_tcr_cross_tissue_workbook_file,
    write_sural_tcr_comparison_workbook(
      sural_tcr_cross_tissue_comparison
    ),
    format = "file"
  ),
  tar_target(
    sural_tcr_intersection_plot_file,
    write_sural_tcr_intersection_plot(
      sural_tcr_cross_tissue_comparison
    ),
    format = "file"
  ),
  tar_target(
    sural_tcr_tracking_plot_file,
    write_sural_tcr_tracking_plot(
      sural_tcr_cross_tissue_comparison
    ),
    format = "file"
  ),
  tar_target(
    sural_tcr_chain_comparisons,
    prepare_sural_tcr_chain_comparisons(
      tcr_contigs,
      combined_tcr,
      sural_trust4_mapped$records,
      analysis_config$trust4$patient_map
    )
  ),
  tar_target(
    sural_tcr_chain_comparison_workbook_file,
    write_sural_chain_comparison_workbook(sural_tcr_chain_comparisons),
    format = "file"
  ),
  tar_target(
    sural_tra_intersection_plot_file,
    write_sural_chain_intersection_plot(
      sural_tcr_chain_comparisons$alpha,
      "TRA",
      "tra_shared_tissue_intersections.pdf"
    ),
    format = "file"
  ),
  tar_target(
    sural_tra_tracking_plot_file,
    write_sural_chain_tracking_plot(
      sural_tcr_chain_comparisons$alpha,
      "TRA",
      "sural_shared_tra_clonotype_tracking.pdf"
    ),
    format = "file"
  ),
  tar_target(
    sural_paired_tcr_intersection_plot_file,
    write_sural_chain_intersection_plot(
      sural_tcr_chain_comparisons$paired,
      "paired TRA + TRB",
      "paired_tra_trb_shared_tissue_intersections.pdf"
    ),
    format = "file"
  ),
  tar_target(
    sural_paired_tcr_tracking_plot_file,
    write_sural_chain_tracking_plot(
      sural_tcr_chain_comparisons$paired,
      "paired TRA + TRB",
      "sural_shared_paired_tra_trb_tracking.pdf"
    ),
    format = "file"
  ),
  tar_target(
    sural_tcr_expansion,
    prepare_sural_tcr_expansion(
      sural_tcr_cross_tissue_comparison,
      sural_tcr_chain_comparisons
    )
  ),
  tar_target(
    sural_tcr_expansion_workbook_file,
    write_sural_tcr_expansion_workbook(sural_tcr_expansion),
    format = "file"
  ),
  tar_target(
    sural_tcr_expansion_plot_file,
    write_sural_tcr_expansion_plot(sural_tcr_expansion),
    format = "file"
  )
)
