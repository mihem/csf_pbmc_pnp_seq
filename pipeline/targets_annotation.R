targets_annotation <- list(
  tar_target(
    annotation_checkpoint_file,
    "references/checkpoints/sc_merge_batch.qs",
    format = "file"
  ),
  tar_target(
    sc_annotation_input,
    apply_annotation_checkpoint(sc_integrated, annotation_checkpoint_file)
  ),
  tar_target(
    sc_clustered,
    cluster_integrated_object(sc_annotation_input, analysis_config)
  ),
  tar_target(
    manual_overrides,
    read_manual_overrides(analysis_config, manual_files)
  ),
  tar_target(sc_annotated, apply_annotations(sc_clustered, manual_overrides, analysis_config)),
  tar_target(
    azimuth_level2_stacas_map_file,
    write_azimuth_level2_map(
      sc_annotated,
      "results/targets/map/map_pbmcref_level2_stacas.ss.all.png",
      seed = 123L
    ),
    format = "file"
  ),
  tar_target(
    integrated_tissue_umap_file,
    write_integrated_tissue_umap(
      sc_annotated,
      "results/targets/umap/stacas.ss.all_umap_tissue.png",
      seed = 123L
    ),
    format = "file"
  ),
  tar_target(
    annotation_marker_file,
    analysis_config$paths$markers,
    format = "file"
  ),
  tar_target(
    annotation_dotplot_cellmarkers_file,
    write_annotation_dotplot(
      sc_annotated,
      annotation_marker_file,
      "cellmarkers_seed",
      "results/targets/dotplot/dp_sc_merge_cellmarkers_seed.pdf",
      10,
      7
    ),
    format = "file"
  ),
  tar_target(
    annotation_dotplot_cd8tem3_file,
    write_annotation_dotplot(
      sc_annotated,
      annotation_marker_file,
      "cd8tem3",
      "results/targets/dotplot/dp_sc_merge_cd8tem3.pdf",
      5.5,
      7
    ),
    format = "file"
  ),
  tar_target(
    annotated_umap,
    write_annotation_umap(sc_annotated, "results/targets/umap/umap_annotated.pdf"),
    format = "file"
  )
)
