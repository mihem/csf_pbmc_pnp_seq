targets_luminex_llod <- list(
  targets::tar_target(
    luminex_llod_analysis, analyze_luminex_llod(luminex_input)
  )
)
