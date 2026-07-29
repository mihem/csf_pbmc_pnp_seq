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
  )
)
