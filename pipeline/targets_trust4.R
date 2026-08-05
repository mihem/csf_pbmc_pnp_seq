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
  )
)
