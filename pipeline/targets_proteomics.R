targets_proteomics <- c(list(
  targets::tar_target(
    proteomics_parameters, proteomics_settings(analysis_config)
  ),
  targets::tar_target(
    proteomics_luminex_marker_inventory_file,
    write_proteomics_luminex_inventory(luminex_llod_analysis),
    format = "file"
  )
  ), tarchetypes::tar_map(
    values = luminex_contrasts(), names = comparison, unlist = TRUE,
    targets::tar_target(
      proteomics_olink_result,
      analyze_proteomics_olink(
        olink_quant_data, olink_quant_input, olink_metadata,
        comparison, group1, group2, proteomics_parameters$fdr
      )
    ),
    targets::tar_target(
      proteomics_olink_files,
      write_proteomics_artifacts(proteomics_olink_result, "olink_quant"),
      format = "file"
    ),
    targets::tar_target(
      proteomics_luminex_result,
      analyze_proteomics_luminex(
        luminex_llod_analysis, comparison, group1, group2,
        proteomics_parameters$fdr
      )
    ),
    targets::tar_target(
      proteomics_luminex_files,
      write_proteomics_artifacts(proteomics_luminex_result, "luminex"),
      format = "file"
    )
  ), list(
    targets::tar_target(
      proteomics_shared_markers,
      dplyr::bind_rows(
        compare_proteomics_markers(
          proteomics_olink_result_GBS_vs_CTRL,
          proteomics_luminex_result_GBS_vs_CTRL,
          c("CCL2", "CCL3", "CXCL8")
        ),
        compare_proteomics_markers(
          proteomics_olink_result_CIDP_vs_CTRL,
          proteomics_luminex_result_CIDP_vs_CTRL,
          c("CCL2", "CCL3", "CXCL8")
        ),
        compare_proteomics_markers(
          proteomics_olink_result_CIDP_vs_GBS,
          proteomics_luminex_result_CIDP_vs_GBS,
          c("CCL2", "CCL3", "CXCL8")
        )
      )
    ),
    targets::tar_target(
      proteomics_shared_marker_file,
      write_proteomics_comparison(proteomics_shared_markers), format = "file"
    )
  ))
