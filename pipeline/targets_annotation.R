targets_annotation_input <- if (use_legacy_validation()) {
  list(
    tar_target(
      sc_annotation_input,
      stabilize_annotation_input(
        sc_integrated,
        legacy_sc_integrated_file,
        validate_sc_integrated_result
      )
    )
  )
} else {
  list(tar_target(sc_annotation_input, sc_integrated))
}

targets_annotation <- c(
  targets_annotation_input,
  list(
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
    cerebro_checkpoint,
    export_cerebro_checkpoint(sc_clustered, "objects/targets/sc_cerebro.crb"),
    format = "file"
  ),
  tar_target(
    annotated_umap,
    write_annotation_umap(sc_annotated, "results/targets/umap/umap_annotated.pdf"),
    format = "file"
  )
  )
)
