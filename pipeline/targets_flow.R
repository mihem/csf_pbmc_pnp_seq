if (!exists("prepare_flow_data", mode = "function")) {
  source(file.path("R", "downstream_flow.R"))
}

targets_flow <- list(
  targets::tar_target(
    flow_raw_file,
    file.path("raw", "flow", "flowbasic_v7.xlsx"),
    format = "file"
  ),
  targets::tar_target(
    flow_frontiers_lookup_file,
    file.path("raw", "flow", "cidp_stats_safe.xls"),
    format = "file"
  ),
  targets::tar_target(
    flow_seed_lookup_file,
    file.path("lookup", "SEED_lookup_v11_anonymized.xlsx"),
    format = "file"
  ),
  targets::tar_target(
    flow_prepared,
    prepare_flow_data(
      flow_raw_file, flow_frontiers_lookup_file, flow_seed_lookup_file
    )
  ),
  targets::tar_target(flow_data, flow_prepared$flow),
  targets::tar_target(flow_match_summary, flow_prepared$match_summary),
  targets::tar_target(
    flow_object_file,
    write_flow_object(
      flow_data, file.path(flow_result_dir(), "flow_pre.qs")
    ),
    format = "file"
  ),
  targets::tar_target(
    flow_gating_comparison_file,
    write_flow_table(
      flow_prepared$gating_comparison,
      file.path(flow_result_dir(), "gating_comparison_full.xlsx")
    ),
    format = "file"
  ),
  targets::tar_target(flow_comparison_config, flow_comparisons()),
  targets::tar_target(
    flow_volcano_results,
    run_flow_volcanoes(
      flow_data,
      flow_comparison_config,
      seed = 123L,
      n_permutations = 1000L
    )
  ),
  targets::tar_target(
    flow_volcano_tables_file,
    write_flow_volcano_tables(
      flow_volcano_results,
      file.path(flow_result_dir(), "volcano_results_all.xlsx")
    ),
    format = "file"
  ),
  targets::tar_target(
    flow_boxplot_files,
    write_flow_boxplots(
      flow_data, file.path(flow_result_dir(), "plots"), seed = 123L
    ),
    format = "file"
  ),
  targets::tar_target(
    flow_volcano_plot_files,
    write_flow_volcano_plots(
      flow_volcano_results,
      file.path(flow_result_dir(), "plots"),
      seed = 123L
    ),
    format = "file"
  ),
  targets::tar_target(
    flow_batch_plot_files,
    write_flow_batch_plots(
      flow_data, file.path(flow_result_dir(), "plots", "batch")
    ),
    format = "file"
  )
)

targets_flow
