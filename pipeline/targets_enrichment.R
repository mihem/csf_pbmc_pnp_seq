targets_enrichment <- list(
  tar_target(
    enrichment_deg_results,
    collect_enrichment_deg_results(
      combined_results = list(
        cidp_ctrl_csf = deg_combined_result_cidp_ctrl_csf,
        gbs_ctrl_csf = deg_combined_result_gbs_ctrl_csf,
        cidp_ctrl_pbmc = deg_combined_result_cidp_ctrl_pbmc,
        gbs_ctrl_pbmc = deg_combined_result_gbs_ctrl_pbmc
      ),
      cluster_results = list(
        cidp_ctrl_csf = deg_cluster_results_cidp_ctrl_csf,
        gbs_ctrl_csf = deg_cluster_results_gbs_ctrl_csf,
        cidp_ctrl_pbmc = deg_cluster_results_cidp_ctrl_pbmc,
        gbs_ctrl_pbmc = deg_cluster_results_gbs_ctrl_pbmc
      )
    )
  ),
  tar_target(
    enrichment_background,
    enrichment_background_genes(sc_annotated)
  ),
  tar_target(
    cd8tem3_markers,
    find_cd8tem3_markers(sc_annotated)
  ),
  tar_target(
    cd8tem3_marker_ora,
    run_cd8tem3_marker_ora(cd8tem3_markers, enrichment_background)
  ),
  tar_target(
    cd8tem3_marker_ora_files,
    write_cd8tem3_marker_ora(cd8tem3_marker_ora),
    format = "file"
  ),
  tar_target(
    enrichment_ora_results,
    run_enrichment_ora(enrichment_deg_results, enrichment_background)
  ),
  tar_target(
    enrichment_gsea_results,
    run_enrichment_gsea(enrichment_deg_results, seed = 123L)
  ),
  tar_target(
    enrichment_ora_files,
    write_enrichment_ora_outputs(enrichment_ora_results),
    format = "file"
  ),
  tar_target(
    enrichment_gsea_files,
    write_enrichment_gsea_outputs(enrichment_gsea_results),
    format = "file"
  )
)
