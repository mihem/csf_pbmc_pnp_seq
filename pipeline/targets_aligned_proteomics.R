targets_aligned_proteomics <- c(list(
  targets::tar_target(
    aligned_bh_parameters, aligned_bh_settings(analysis_config)
  ),
  targets::tar_target(
    proteomics_luminex_marker_inventory_file,
    write_aligned_bh_luminex_inventory(luminex_llod_analysis),
    format = "file"
  )
  ), tarchetypes::tar_map(
    values = luminex_contrasts(), names = comparison, unlist = TRUE,
    targets::tar_target(
      aligned_bh_olink_result,
      analyze_aligned_bh_olink(
        olink_quant_diagnosis, olink_quant_input, olink_metadata,
        comparison, group1, group2, aligned_bh_parameters$fdr
      )
    ),
    targets::tar_target(
      aligned_bh_olink_files,
      write_aligned_bh_artifacts(aligned_bh_olink_result, "olink_quant"),
      format = "file"
    ),
    targets::tar_target(
      aligned_bh_luminex_result,
      analyze_aligned_bh_luminex(
        luminex_llod_analysis, comparison, group1, group2,
        aligned_bh_parameters$fdr
      )
    ),
    targets::tar_target(
      aligned_bh_luminex_files,
      write_aligned_bh_artifacts(aligned_bh_luminex_result, "luminex"),
      format = "file"
    )
  ), list(
    targets::tar_target(
      aligned_bh_shared_markers,
      dplyr::bind_rows(
        compare_aligned_bh_markers(
          aligned_bh_olink_result_GBS_vs_CTRL,
          aligned_bh_luminex_result_GBS_vs_CTRL,
          c("CCL2", "CCL3", "CXCL8")
        ),
        compare_aligned_bh_markers(
          aligned_bh_olink_result_CIDP_vs_CTRL,
          aligned_bh_luminex_result_CIDP_vs_CTRL,
          c("CCL2", "CCL3", "CXCL8")
        ),
        compare_aligned_bh_markers(
          aligned_bh_olink_result_CIDP_vs_GBS,
          aligned_bh_luminex_result_CIDP_vs_GBS,
          c("CCL2", "CCL3", "CXCL8")
        )
      )
    ),
    targets::tar_target(
      aligned_bh_shared_marker_file,
      write_aligned_bh_comparison(aligned_bh_shared_markers), format = "file"
    )
  ))
