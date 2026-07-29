targets_integration <- list(
  tar_target(
    stacas_ss_all,
    compute_stacas_reduction(sc_azimuth, lookup, analysis_config)
  ),
  tar_target(
    sc_integrated,
    assemble_stacas_integration(
      sc_azimuth,
      lookup,
      stacas_ss_all,
      analysis_config
    )
  ),
  tar_target(
    azimuth_level2_stacas_map_file,
    write_azimuth_level2_map(
      sc_integrated,
      "results/targets/map/map_pbmcref_level2_stacas.ss.all.png",
      seed = 123L
    ),
    format = "file"
  )
)
